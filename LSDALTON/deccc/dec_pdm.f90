!> @file
!> Subroutines used in scalapack related context within dec
!> the scalapack within dec is thought to adopt most of the
!> properties already defined in matop_scalapack.f90
!> \author Patrick Ettenhuber
!> \date May 2012
module dec_pdm_module


  ! Outside DEC directory
  use precision
  use ptr_assoc_module, only: ass_D1to4
  use dec_typedef_module
  use matrix_operations_scalapack!, only:SLGrid,BLOCK_SIZE
  use matrix_operations!, only:matrix_type,mtype_scalapack
  use memory_handling!, only: mem_alloc, mem_dealloc
#ifdef VAR_LSMPI
  use infpar_module
  use lsmpi_type
#endif

  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
  use manual_utils_module
  use array_memory_manager
  use dec_fragment_utils


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
    type(array) :: a(n_arrays)
    !> current address on node
    integer :: curr_addr_on_node=1
    !> counter for how many arrays were allocated in the persisten array
    integer :: arrays_allocated = 0
    !> counter for how many arrays were deallocated
    integer :: arrays_deallocated = 0
    !> conter for the arrays currently in use
    integer :: arrays_in_use = 0
    !> offset for the first tile allocation to get a better load distribution
    integer :: new_offset = 0
    !> list of n logicals as indicator wheter an adress is free to allocate a
    !new array
    logical :: free_addr_on_node(n_arrays)=.true.
  endtype persistent_array

  save
  
  ! job parameters for pdm jobs
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
  integer,parameter :: JOB_CHANGE_INIT_TYPE    = 22

  !> definition of the persistent array 
  type(persistent_array) :: p_arr


  contains

  !> \brief main subroutine for the communication of nodes on grid handling arr structures
  !> \author Patrick Ettenhuber
  !> \date May 2012
  subroutine pdm_array_sync(job,a,b,c,d)
    implicit none
    !> job is input for master and output for slaves, the arguments have to be
    !in the job paramenters list in top of this file
    integer :: job
    !the array(s) to be passed to the slaves for which the operation is
    !performed
    type(array),optional :: a,b,c,d
    integer,pointer ::   TMPI(:), dims(:)
    !character :: TMPC(12)
    integer :: i, j, context,modes(3),counter, stat,ierr,basic
    integer(kind=ls_mpik) :: sendctr
    modes=0
#ifdef VAR_LSMPI

    basic = 12


    IF(infpar%lg_mynum.eq.infpar%master) then
      !**************************************************************************************
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!code for MASTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !**************************************************************************************
      !Wake up slaves
      call ls_mpibcast(PDMA4SLV,infpar%master,infpar%lg_comm)
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
      call ls_mpibcast(counter,infpar%master,infpar%lg_comm)
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
      if(infpar%lg_nodtot>1.and.present(A).and..not.associated(A%addr_p_arr))&
      &call lsquit("ERROR(pdm_arr_sync):addr_p_arr for array A not associated",DECinfo%output)
      if(infpar%lg_nodtot>1.and.present(B).and..not.associated(B%addr_p_arr))&
      &call lsquit("ERROR(pdm_arr_sync):addr_p_arr for array B not associated",DECinfo%output)
      if(infpar%lg_nodtot>1.and.present(C).and..not.associated(C%addr_p_arr))&
      &call lsquit("ERROR(pdm_arr_sync):addr_p_arr for array C not associated",DECinfo%output)
      if(infpar%lg_nodtot>1.and.present(D).and..not.associated(D%addr_p_arr))&
      &call lsquit("ERROR(pdm_arr_sync):addr_p_arr for array D not associated",DECinfo%output)
      
      !print *,"master",TMPI      

      do sendctr=1,infpar%lg_nodtot-1
        IF (PRESENT(A)) TMPI(2)  = A%addr_p_arr(sendctr)
        IF (PRESENT(B)) TMPI(3)  = B%addr_p_arr(sendctr)
        IF (PRESENT(C)) TMPI(4)  = C%addr_p_arr(sendctr)
        IF (PRESENT(D)) TMPI(5)  = D%addr_p_arr(sendctr)
        call ls_mpisendrecv(TMPI,counter,infpar%lg_comm,infpar%master,sendctr)
      enddo
      call mem_dealloc(TMPI)


    else  

      !**************************************************************************************
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!code for SLAVES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !**************************************************************************************
      call ls_mpibcast(counter,infpar%master,infpar%lg_comm)
      call mem_alloc(TMPI,counter)
      call ls_mpisendrecv(TMPI,counter,infpar%lg_comm,infpar%master,infpar%lg_mynum)

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

#ifdef VAR_LSMPI
    !Get the slaves to this routine
    if(infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(JOB_GET_FRAG_CC_ENERGY,t1,t2,gmo)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(occ_num,infpar%master)
      call ls_mpi_buffer(occ_idx,occ_num,infpar%master)
      call ls_mpi_buffer(virt_num,infpar%master)
      call ls_mpi_buffer(virt_idx,virt_num,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    endif
    call memory_allocate_array_dense(gmo)
    call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)

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

    call memory_deallocate_array_dense(gmo)
    
    call lsmpi_local_reduction(Eocc,infpar%master)
    call lsmpi_local_reduction(Evirt,infpar%master)

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
    integer :: lt,i,j,a,b,o(t2%mode)

#ifdef VAR_LSMPI
    !Get the slaves to this routine
    if(infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(JOB_GET_CC_ENERGY,t1,t2,gmo)
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

      !count over local indices
      do j=1,t2%ti(lt)%d(4)
        do i=1,t2%ti(lt)%d(3)
          do b=1,t2%ti(lt)%d(2)
            do a=1,t2%ti(lt)%d(1)
     
              E2 = E2 + t(a,b,i,j)*&
              & (2.0E0_realk*  gmo%elm4(i+o(3),a+o(1),j+o(4),b+o(2))-gmo%elm4(i+o(3),b+o(2),j+o(4),a+o(1)))
              E1 = E1 + ( t1%elm2(a+o(1),i+o(3))*t1%elm2(b+o(2),j+o(4)) ) * &
                   (2.0E0_realk*gmo%elm4(i+o(3),a+o(1),j+o(4),b+o(2))-gmo%elm4(i+o(3),b+o(2),j+o(4),a+o(1)))
   
            enddo 
          enddo
        enddo
      enddo
      nullify(t)
    enddo

    call memory_deallocate_array_dense(gmo)
    
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
    integer :: t(4)

#ifdef VAR_LSMPI

    !CHECK if the distributions are the same, if it becomes necessary, that they
    !are not, then this routine has to be rewritten
    if(omega2%tdim(1)/=prec%tdim(1).or.omega2%tdim(2)/=prec%tdim(2).or.&
    &omega2%tdim(3)/=prec%tdim(3).or.omega2%tdim(4)/=prec%tdim(4))then
      call lsquit("ERROR(precondition_doubles_parallel):omega2 and prec have&
      &different distributions", DECinfo%output)
    endif
    !Get the slaves to this routine
    if(infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(JOB_PREC_DOUBLES_PAR,omega2,ppfock,qqfock,prec)
    endif

    dims=prec%dims

    
    !do a loop over the local tiles of the preconditioned matrix and get the
    !corresponding tiles of the residual to form the preconditioned residual
    do lt=1,prec%nlti
      call array_get_tile(omega2,prec%ti(lt)%gt,prec%ti(lt)%t,prec%ti(lt)%e)
      call ass_D1to4(prec%ti(lt)%t,om,prec%ti(lt)%d)
      
      !get offset for global indices
      call get_midx(prec%ti(lt)%gt,dims,prec%ntpm,prec%mode)
      do j=1,prec%mode
        dims(j)=(dims(j)-1)*prec%tdim(j)
      enddo

      !count over local indices
      do j=1,prec%ti(lt)%d(4)
        do i=1,prec%ti(lt)%d(3)
          do b=1,prec%ti(lt)%d(2)
            do a=1,prec%ti(lt)%d(1)
     
              om(a,b,i,j) = om(a,b,i,j) / &
                 ( ppfock%elm2(i+dims(3),i+dims(3)) - qqfock%elm2(a+dims(1),a+dims(1)) + &
                   ppfock%elm2(j+dims(4),j+dims(4)) - qqfock%elm2(b+dims(2),b+dims(2)) )
   
            enddo 
          enddo
        enddo
      enddo
      nullify(om)
    enddo
    
    !crucial barrier, wait for all slaves to finish their jobs
    call lsmpi_barrier(infpar%lg_comm)
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

#ifdef VAR_LSMPI
    !check if the init-types are the same
    if(arr1%init_type/=arr2%init_type)then
      call lsquit("ERROR(array_ddot_par):different init types of the&
      & arrays is not possible",DECinfo%output)
    endif

    !check if the destination to collet the resut makes sense in connection with
    !the init_type
    if(arr1%init_type==MASTER_INIT.and.dest/=0)then
      call lsquit("ERROR(array_ddot_par): the choice of destnation is&
      & useless",DECinfo%output)
    endif

    !get the slaves to this routine
    if(arr1%init_type==MASTER_INIT.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(JOB_DDOT_PAR,arr1,arr2)
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
      call lsmpi_local_allreduce(res)
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
#ifdef VAR_LSMPI

    !check if the init_types are the same
    if(x%init_type/=y%init_type)then
      call lsquit("ERROR(array_add_par):different init types&
      & impossible",DECinfo%output)
    endif

    !IF NOT MASTER_INIT all processes should know b on call-time, else b is
    !broadcasted here
    if(x%init_type==MASTER_INIT.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(JOB_ADD_PAR,x,y)
      call ls_mpibcast(b,infpar%master,infpar%lg_comm)
    elseif(x%init_type==MASTER_INIT.and.infpar%lg_mynum/=infpar%master)then
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
#ifdef VAR_LSMPI

    !check for the same init_types
    if(from%init_type/=to_ar%init_type)then
      call lsquit("ERROR(array_cp_tiled):different init types&
      & impossible",DECinfo%output)
    endif

    !get the slaves
    if(from%init_type==MASTER_INIT.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(JOB_CP_ARR,from,to_ar)
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
#ifdef VAR_LSMPI
    !get the slaves here
    if(a%init_type==MASTER_INIT.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(JOB_ARRAY_ZERO,a)
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
    !> integer specifying the init_type of the array
    integer,intent(in) :: pdm
    integer(kind=long) :: i,j
    integer :: rnk,stat,help,addr,tdimdummy(nmodes)
    integer :: dflt(nmodes),nelms
    integer(kind=ls_mpik) :: nlocalnodes,ierr
    integer, pointer :: buf(:)
    logical :: master

    !allocate all pdm in p_arr therefore get free address and associate it with
    !the array, and increment the array counter
    p_arr%curr_addr_on_node=get_free_address(.true.)
    addr=p_arr%curr_addr_on_node
    p_arr%arrays_in_use = p_arr%arrays_in_use + 1

    !set the initial values and overwrite them later
    nlocalnodes=1
    master=.true.
    p_arr%a(addr)%init_type=pdm

#ifdef VAR_LSMPI
    !assign if master and the number of nodes in the local group
    if(.not.infpar%lg_mynum==infpar%master)master=.false.
    nlocalnodes=infpar%lg_nodtot
#endif

    !SET MODE
    p_arr%a(addr)%mode=nmodes
    !SET DIMS
    call arr_set_dims(p_arr%a(addr),dims,nmodes)
    !SET ARRAY TYPE
    p_arr%a(addr)%atype=REPLICATED
    !SET NELMS
    nelms=1
    do i=1,nmodes
      nelms=nelms*dims(i)
    enddo
    p_arr%a(addr)%nelms=nelms

    !put 0 in tdim, since for the replicated array it is not important
    tdimdummy=0
    call arr_set_tdims(p_arr%a(addr),tdimdummy,nmodes)

    !In the initialization the addess has to be set, since pdm_array_sync
    !depends on the  adresses, but setting them correctly is done later
    if(p_arr%a(addr)%init_type>=1)then
      call mem_alloc(buf,nlocalnodes)
      buf = 0
    endif
    
    !if master init only master has to init the addresses addresses before
    !pdm syncronization
    if(master .and. p_arr%a(addr)%init_type==MASTER_INIT)then
      call arr_set_addr(p_arr%a(addr),buf,nlocalnodes)
#ifdef VAR_LSMPI
      call pdm_array_sync(JOB_INIT_ARR_REPLICATED,p_arr%a(addr))
#endif
    endif

    !if all_init all have to have the addresses allocated
    if(p_arr%a(addr)%init_type==ALL_INIT)call arr_set_addr(p_arr%a(addr),buf,nlocalnodes)

#ifdef VAR_LSMPI
    !SET THE ADDRESSES ON ALL NODES     
    buf(infpar%lg_mynum+1)=addr 
    call lsmpi_local_allreduce(buf,nlocalnodes)
    call arr_set_addr(p_arr%a(addr),buf,nlocalnodes)
#endif
  
    !ALLOCATE STORAGE SPACE FOR THE ARRAY
    call memory_allocate_array_dense(p_arr%a(addr))

    !RETURN THE CURRENLY ALLOCATE ARRAY
    arr=p_arr%a(addr)

    if(arr%init_type>=1) call mem_dealloc(buf)
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
#ifdef VAR_LSMPI

    !get the slaves
    if(infpar%lg_mynum==infpar%master.and.arr%init_type==MASTER_INIT)then
      call pdm_array_sync(JOB_GET_NORM_REPLICATED,arr)
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
#ifdef VAR_LSMPI

    !give meaningful quit statement for useless input
    if(present(fromnode).and.arr%init_type==MASTER_INIT)then
      call lsquit("ERROR(array_sync_replicated): This combintion of input&
      &elements does not give sense",DECinfo%output)
      ! why would you want to collect the data on a node you cannot direcly
      ! access, or if you can access the data in the calling subroutine on the
      ! specified node, why is the init_tyep MASTER_INIT?
    endif

    ! get slaves
    if(infpar%lg_mynum==infpar%master.and.arr%init_type==MASTER_INIT)then
      call pdm_array_sync(JOB_SYNC_REPLICATED,arr)
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
#ifdef VAR_LSMPI
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
  function array_init_tiled(dims,nmodes,pdm,tdims,zeros_in_tiles)result(arr)
    implicit none
    type(array) :: arr
    integer,intent(in) :: nmodes,dims(nmodes)
    integer,optional :: tdims(nmodes),pdm
    logical, optional :: zeros_in_tiles
    integer(kind=long) :: i,j
    integer :: rnk,stat,ierr,help,addr,pdmt,k,div
    integer :: dflt(nmodes),cdims
    integer, pointer :: buf(:)
    integer(kind=ls_mpik) :: nlocalnodes
    logical :: master,defdims
    !allocate all tiled arrays in p_arr, get free
    p_arr%curr_addr_on_node=get_free_address(.true.)

    addr=p_arr%curr_addr_on_node

    p_arr%arrays_in_use = p_arr%arrays_in_use + 1

    nlocalnodes=1
    master=.true.
    !call array_free_basic(arr)

    if(present(pdm))then
      p_arr%a(addr)%init_type=pdm
    else
      p_arr%a(addr)%init_type=NO_PDM
    endif
#ifdef VAR_LSMPI
    if(.not.infpar%lg_mynum==infpar%master)master=.false.
    !if(p_arr%a(addr)%init_type==ALL_INIT)call lsmpi_barrier(infpar%lg_comm)
    nlocalnodes=infpar%lg_nodtot
#endif

    !INITIALIZE TILE STRUCTURE, if master from basics, if slave most is already
    !there
    defdims=.false.
    if(present(zeros_in_tiles)) p_arr%a(addr)%zeros=zeros_in_tiles

    p_arr%a(addr)%mode=nmodes

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
    if(p_arr%a(addr)%init_type>=1) p_arr%a(addr)%atype=TILED_DIST
    if(p_arr%a(addr)%init_type==0) p_arr%a(addr)%atype=TILED
    
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

    p_arr%a(addr)%ntiles=1
    p_arr%a(addr)%tsize=1
    do i=1,p_arr%a(addr)%mode
      p_arr%a(addr)%ntiles = p_arr%a(addr)%ntiles * p_arr%a(addr)%ntpm(i)
      p_arr%a(addr)%tsize  = p_arr%a(addr)%tsize  * min(p_arr%a(addr)%tdim(i),p_arr%a(addr)%dims(i))
    enddo

    !In the initialization the addess has to be set, since pdm_array_sync
    !depends on the  adresses, but setting them correctly is done later
    if(p_arr%a(addr)%init_type>=1)then
      call mem_alloc(buf,2*nlocalnodes)
      buf = 0
    endif
    
    !if master init only master has to get addresses
    if(master .and. p_arr%a(addr)%init_type==MASTER_INIT)then
      call arr_set_addr(p_arr%a(addr),buf,nlocalnodes)
      call pdm_array_sync(JOB_INIT_ARR_TILED,p_arr%a(addr))
    endif

    !if all_init only all have to know the addresses
    if(p_arr%a(addr)%init_type==ALL_INIT)call arr_set_addr(p_arr%a(addr),buf,nlocalnodes)
    !print *,infpar%mynum,"set things now entering tile alloc",arr%init_type


    if(p_arr%a(addr)%init_type>=1)then
#ifdef VAR_LSMPI
      call get_distribution_info(p_arr%a(addr))
      buf(infpar%lg_mynum+1)=addr 
      buf(nlocalnodes+infpar%lg_mynum+1)=p_arr%a(addr)%offset
      call lsmpi_local_allreduce(buf,2*nlocalnodes)
      call arr_set_addr(p_arr%a(addr),buf,nlocalnodes)
      do i=1,nlocalnodes
        if(buf(nlocalnodes+i)/=p_arr%a(addr)%offset)then
          print * ,infpar%lg_mynum,"found",buf(nlocalnodes+i),p_arr%a(addr)%offset,i
          call lsquit("ERROR(array_init_tiled):offset &
          &is not the same on all nodes",DECinfo%output)
        endif
      enddo
#endif
    else
      p_arr%a(addr)%nlti=p_arr%a(addr)%ntiles
    endif

    call memory_allocate_tiles(p_arr%a(addr))

    arr=p_arr%a(addr)
    !print *,infpar%lg_mynum,associated(arr%wi),"peristent",associated(p_arr%a(addr)%wi)

    if(arr%init_type>=1) call mem_dealloc(buf)
    !print *,infpar%lg_mynum,"init done returning"
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
    integer, intent(in) :: nelms
    !> logical specifying whether the tiles are in pdm
    logical, intent(in) :: pdm
    !> reorder if the array is reorder with respect to the fortran array
    integer, intent(in), optional :: order(arr%mode)
    integer :: i,j,k,tmdidx(arr%mode)
    integer :: l,nelintile,tdim(arr%mode)
    real(realk), pointer :: tmp(:)

    !check nelms
    if(nelms/=arr%nelms)call lsquit("ERROR(cp_tileddate2fort):array&
        &dimensions are not the same",DECinfo%output)

    !allocate space 
    if(pdm)then
      call mem_alloc(tmp,arr%tsize)
    endif

    do i=1,arr%ntiles
      call get_midx(i,tmdidx,arr%ntpm,arr%mode)
      call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
      if(pdm)then
        call array_get_tile(arr,i,tmp,nelintile)
      else
        tmp => arr%ti(i)%t
      endif
      if(present(order))call add_tile_to_fort(tmp,i,arr%tdim,fort,arr%dims,arr%mode,b,order)
      if(.not.present(order))call add_tile_to_fort(tmp,i,arr%tdim,fort,arr%dims,arr%mode,b)
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
    integer, intent(in) :: nelms
    logical, intent(in) :: pdm
    integer, intent(in), optional :: order(arr%mode)
    integer :: i,j,k,tmdidx(arr%mode),minimode(arr%mode),o(arr%mode)
    integer :: glbmodeidx(arr%mode),glbidx,l,nelintile,tdim(arr%mode)
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

    do i=1,arr%ntiles
      call get_midx(i,tmdidx,arr%ntpm,arr%mode)
      call get_tile_dim(l,i,arr%dims,arr%tdim,arr%mode,2)
      call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
      if(pdm)then
        call array_get_tile(arr,i,tmp,nelintile)
      else
        tmp => arr%ti(i)%t
      endif
      call put_tile_in_fort(tmp,i,arr%tdim,fort,arr%dims,arr%mode,o)
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
    integer,intent(in) :: nelms
    integer(kind=ls_mpik),intent(in) :: nod
    integer,intent(in),optional :: optorder(arr%mode)
    real(realk),pointer :: buf(:)
    integer :: nelmsit,i, order(arr%mode)
    integer :: ltidx
    integer(kind=ls_mpik) :: nnod,dest,me
#ifdef VAR_LSMPI

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE
    do i=1,arr%mode
      order(i)=i
    enddo
    if(present(optorder))order=optorder
   

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
        call extract_tile_from_fort(A,arr%mode,i,arr%dims,arr%tdim,buf,order)
        call daxpy(nelmsit,sc,buf,1,arr%ti(ltidx)%t,1)
      elseif(nod==me)then
        call extract_tile_from_fort(A,arr%mode,i,arr%dims,arr%tdim,buf,order)
        call lsmpi_send(buf,nelmsit,infpar%lg_comm,dest)
      elseif(dest==me)then
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



  subroutine array_gatheradd_tilestofort(arr,sc,fort,nelms,nod,optorder)
    implicit none
    type(array),intent(in) :: arr
    real(realk),intent(in) :: sc
    real(realk),intent(inout) :: fort(*)
    integer, intent(in) :: nelms
    integer(kind=ls_mpik) :: nod
    integer, intent(in), optional :: optorder(arr%mode)
    integer(kind=ls_mpik) :: src,me,nnod
    integer :: i,ltidx,order(arr%mode)
    integer :: nelintile
    real(realk), pointer :: tmp(:)
#ifdef VAR_LSMPI

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

    call mem_alloc(tmp,arr%tsize)

    do i=1,arr%ntiles
      src=get_residence_of_tile(i,arr)
      if(src==me.or.nod==me)then
        call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
        if(src==me.and.nod==me)then
          ltidx = (i - 1) /nnod + 1
          call add_tile_to_fort(arr%ti(ltidx)%t,i,arr%tdim,fort,arr%dims,arr%mode,sc,order)
        elseif(src==me)then
          ltidx = (i - 1) /nnod + 1
          call lsmpi_send(arr%ti(ltidx)%t,nelintile,infpar%lg_comm,nod)
        elseif(nod==me)then
          call lsmpi_recv(tmp,nelintile,infpar%lg_comm,src)
          call add_tile_to_fort(tmp,i,arr%tdim,fort,arr%dims,arr%mode,sc,order)
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
    integer, intent(in) :: nelms
    integer(kind=ls_mpik) :: nod
    integer, intent(in), optional :: optorder(arr%mode)
    integer(kind=ls_mpik) :: src,me,nnod
    integer :: i,j,k,ltidx
    integer :: nelintile,order(arr%mode)
    real(realk), pointer :: tmp(:)
#ifdef VAR_LSMPI
    do i = 1, arr%mode
      order(i) = i
    enddo
    if(present(optorder))order=optorder
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
          call put_tile_in_fort(arr%ti(ltidx)%t,i,arr%tdim,fort,arr%dims,arr%mode,order)
        elseif(src==me)then
          ltidx = (i - 1) /nnod + 1
          call lsmpi_send(arr%ti(ltidx)%t,nelintile,infpar%lg_comm,nod)
        elseif(nod==me)then
          call lsmpi_recv(tmp,nelintile,infpar%lg_comm,src)
          call put_tile_in_fort(tmp,i,arr%tdim,fort,arr%dims,arr%mode,order)
        endif
      endif
    enddo

    call mem_dealloc(tmp)
#else
    call lsquit("ERROR(array_gather_tilesinfort):this routine is MPI only",-1)
#endif
  end subroutine array_gather_tilesinfort


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
    pdm=(arr%atype==TILED_DIST)

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE

    me = 0
    nnod=1
#ifdef VAR_LSMPI
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
      !if(comp_ti>arr%ntiles)then
      !  call lsquit("tiles out of range",DECinfo%output)
      !endif
      !check where the current tile resides and jump the following steps if not
      !master where the full matrix resides or the destination slave
      !dest = mod(comp_ti-1+arr%offset,nnod) 
      if(pdm)dest = get_residence_of_tile(comp_ti,arr) 
      !print *,comp_ti,arr%offset,dest,infpar%lg_nodtot,i,nccblocks
      !if(dest>infpar%lg_nodtot)then
      !  call lsquit("destination out of range",DECinfo%output)
      !endif
    
      !get the dimensions of the remote tile
      call get_tile_dim(remote_td,comp_ti,arr%dims,arr%tdim,arr%mode)
      !do j=1, arr%mode
      !  if(((arr%dims(j)-(in_tile_mode(j)-1)*arr%tdim(j))/arr%tdim(j))>=1)then
      !    remote_td(j)=arr%tdim(j)
      !  elseif(((arr%dims(j)-(in_tile_mode(j)-1)*arr%tdim(j))/arr%tdim(j))<1 )then
      !    remote_td(j)=mod(arr%dims(j)-(in_tile_mode(j)-1)*arr%tdim(j),arr%tdim(j))
      !  endif
      !enddo



      !now get position of the first element of the batch in the current tile
      do j=1,arr%mode
        elm_in_tile(j) = mod(orig_addr(j)-1,arr%tdim(j)) + 1
      enddo

      !get the one index element number for the remote tile
      comp_el=get_cidx(elm_in_tile,remote_td,arr%mode)

      !copy data to the identified places
      if(pdm)then
#ifdef VAR_LSMPI
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
#ifdef VAR_LSMPI
    master = .true.
    if(infpar%lg_mynum/=infpar%master)master=.false.

    if(.not.present(infoonmaster))then
      if(master.and..not.allaccs)then
        call pdm_array_sync(JOB_PRINT_MEM_INFO1)
      endif
      do i=1,infpar%lg_nodtot
        if(infpar%lg_mynum+1==i)then
          write(output,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
          write(output,'("Printing memory information for rank",I3)'),infpar%lg_mynum
          write(output,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
          call array_print_memory_currents(output)
          write(output,'("currently",I4," arrays allocated")')p_arr%arrays_in_use
        endif
        call lsmpi_barrier(infpar%lg_comm)
      enddo
    else
      if(master.and..not.allaccs)then
        call pdm_array_sync(JOB_PRINT_MEM_INFO2)
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
#endif
  end subroutine print_mem_per_node

  subroutine add_data2tiled_intiles(arr,mult,A,dims,mode,order)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*),mult
    integer,intent(in) :: mode, dims(mode)
    integer,intent(in),optional :: order(mode)
    real(realk),pointer :: buf(:)
    integer ::nnod,fib,lt,ce,j,me,dest,step,act_step,mod_step,iter,nccblocks,ierr,st
    integer :: nelmsit,loc_ti,comp_ti,comp_el,i,nelms,fe_in_block,o(mode)
    integer, pointer :: elm_in_tile(:),in_tile_mode(:),orig_addr(:),remote_td(:)
    do i=1,mode
      o(i)=i
    enddo
    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE

    me = 0
    nnod=1
#ifdef VAR_LSMPI
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
      if(.not.present(order))call extract_tile_from_fort(A,arr%mode,i,arr%dims,arr%tdim,buf,o)
      if(present(order))call extract_tile_from_fort(A,arr%mode,i,arr%dims,arr%tdim,buf,order)
      call get_tile_dim(nelmsit,arr,i)
      if(mult/=1.0E0_realk)call dscal(nelmsit,mult,buf,1)
#ifdef VAR_LSMPI
      call array_accumulate_tile(arr,i,buf,nelmsit)
#endif
    enddo
    call mem_dealloc(buf)
  end subroutine add_data2tiled_intiles

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

    pdm=(arr%atype==TILED_DIST)

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE

    me = 0
    nnod=1
#ifdef VAR_LSMPI
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
      !if(comp_ti>arr%ntiles)then
      !  call lsquit("tiles out of range",DECinfo%output)
      !endif
      !check where the current tile resides and jump the following steps if not
      !master where the full matrix resides or the destination slave
      !dest = mod(comp_ti-1+arr%offset,nnod) 
      if(pdm)dest = get_residence_of_tile(comp_ti,arr) 
      !print *,comp_ti,arr%offset,dest,infpar%lg_nodtot,i,nccblocks
      !if(dest>infpar%lg_nodtot)then
      !  call lsquit("destination out of range",DECinfo%output)
      !endif
    
      !get the dimensions of the remote tile
      call get_tile_dim(remote_td,comp_ti,arr%dims,arr%tdim,arr%mode)
      !do j=1, arr%mode
      !  if(((arr%dims(j)-(in_tile_mode(j)-1)*arr%tdim(j))/arr%tdim(j))>=1)then
      !    remote_td(j)=arr%tdim(j)
      !  elseif(((arr%dims(j)-(in_tile_mode(j)-1)*arr%tdim(j))/arr%tdim(j))<1 )then
      !    remote_td(j)=mod(arr%dims(j)-(in_tile_mode(j)-1)*arr%tdim(j),arr%tdim(j))
      !  endif
      !enddo



      !now get position of the first element of the batch in the current tile
      do j=1,arr%mode
        elm_in_tile(j) = mod(orig_addr(j)-1,arr%tdim(j)) + 1
      enddo

      !get the one index element number for the remote tile
      comp_el=get_cidx(elm_in_tile,remote_td,arr%mode)

      !copy data to the identified places
      if(pdm)then
#ifdef VAR_LSMPI
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

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE
    do i=1,mode
      order(i)=i
    enddo
    if(present(optorder))order=optorder
   

    me = 0
    nnod=1
#ifdef VAR_LSMPI
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
      call extract_tile_from_fort(A,arr%mode,i,arr%dims,arr%tdim,buf,order)
      call get_tile_dim(nelmsit,arr,i)
      !call pn(buf,nelmsit)
      !copy data to the identified places
#ifdef VAR_LSMPI
      !print *,"copying",i,arr%ntiles
      call array_put_tile(arr,i,buf,nelmsit)
#endif
    enddo
    call mem_dealloc(buf)
  end subroutine cp_data2tiled_intiles


  subroutine array_scatter_densetotiled(arr,A,nelms,nod,optorder)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*)
    integer,intent(in) :: nelms
    integer(kind=ls_mpik),intent(in) :: nod
    integer,intent(in),optional :: optorder(arr%mode)
    real(realk),pointer :: buf(:)
    integer :: nelmsit,i, order(arr%mode)
    integer :: ltidx
    integer(kind=ls_mpik) :: nnod,dest,me
#ifdef VAR_LSMPI

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE
    do i=1,arr%mode
      order(i)=i
    enddo
    if(present(optorder))order=optorder
   

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
        call extract_tile_from_fort(A,arr%mode,i,arr%dims,arr%tdim,arr%ti(ltidx)%t,order)
      elseif(nod==me)then
        call extract_tile_from_fort(A,arr%mode,i,arr%dims,arr%tdim,buf,order)
        call lsmpi_send(buf,nelmsit,infpar%lg_comm,dest)
      elseif(dest==me)then
        ltidx = (i - 1) /nnod + 1
        call lsmpi_recv(arr%ti(ltidx)%t,nelmsit,infpar%lg_comm,nod)
      endif
    enddo
    call mem_dealloc(buf)
#else
    call lsquit("ERROR(array_scatter_densetotiled):this routine is MPI only",-1)
#endif

  end subroutine array_scatter_densetotiled



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


  !> \autor Patrick Ettenhuber
  !> \brief Subroutine to put a specific tile in a general matrix
  !> \date January 2013
  subroutine add_tile_to_fort(tilein,tnr,tdims,fort,dims,mode,sc,optorder)
    implicit none
    !> input array with data in arbtirary order
    real(realk), intent(inout) :: fort(*)
    !> mode infortmation about how to interpret data
    integer, intent(in) :: mode
    !> the tile number in column major ordering of the tiles
    integer, intent(in) :: tnr
    !> dimension infortmation about how to interpret data
    integer, intent(in) :: dims(mode)
    !> batch information for the tiles
    integer, intent(in) :: tdims(mode)
    !> logical specifying wheter to copy or to add
    real(realk),intent(in) :: sc
    !> specify the order if the tile stems from a matrix with a different order
    integer,intent(in),optional :: optorder(mode)
    !> tile output
    real(realk), intent(in) :: tilein(*)
    integer :: i,nccblcks,nelms,k
    integer :: tmodeidx(mode),o(mode)
    integer :: idxintile(mode),glbidx,rorder(mode),olddims(mode)
    integer :: ccels,ntimes,el,acttdim(mode),nels,fels(mode)
    integer :: pos1,pos2,ntpm(mode),glbmodeidx(mode)
    integer :: simpleorder,bs
    bs=int(((8000.0*1000.0)/(8.0*2.0))**(1.0/float(mode)))
    !bs=5
    nels=1
    do i=1,mode
      o(i)=i
      nels=nels*dims(i)
    enddo
    simpleorder=0

    if(present(optorder))o=optorder
    
    do i=1,mode
      rorder(o(i))=i
      if(o(i)/=i)simpleorder=-1
      ntpm(i) = dims(i)/tdims(i)
      if(mod(dims(i),tdims(i))>0)ntpm(i) = ntpm(i) + 1
    enddo
    if(mode==4)then
      if(o(1)==2.and.o(2)==1.and.o(3)==4.and.o(4)==3)simpleorder=1
    endif
    call get_midx(tnr,tmodeidx,ntpm,mode)

    ntimes=1
    do i=1,mode
      fels(i) = (tmodeidx(i)-1) * tdims(i) + 1
      olddims(i) = dims(rorder(i))
      if(tmodeidx(i)*tdims(i)>dims(i))then
        acttdim(i)=mod(dims(i),tdims(i))
      else
        acttdim(i)=tdims(i)
      endif
      if(i>1)ntimes=ntimes*acttdim(i)
    enddo
    

    select case(simpleorder)
      case(0)
        ccels=acttdim(1)
        do i=1,ntimes
          call get_midx(i,idxintile(2:mode),acttdim(2:mode),mode-1)
          idxintile(1)=1
          do k=1,mode
            glbmodeidx(k)=idxintile(k) +(tmodeidx(k)-1) *tdims(k)
          enddo
          pos1=get_cidx(glbmodeidx,dims,mode)
          call daxpy(ccels,sc,tilein(1+(i-1)*ccels),1,fort(pos1),1)
        enddo
      case(1)
        call manual_2143_reordering_tile2full(bs,tdims,dims,fels,sc,tilein,1.0E0_realk,fort)
      case default
        print *,"default part reorder add",o
        !count elements in the current tile for loop over elements
        !identify their original position and put them in tile
        nelms=1
        do i=1,mode
          nelms = nelms * acttdim(i)
        enddo
        do i = 1,nelms
          !get mode index of element in tile
          call get_midx(i,idxintile,acttdim,mode)
          do k=1,mode
            glbmodeidx(o(k))=idxintile(k) + (tmodeidx(k)-1)*tdims(k)
          enddo
          pos1=get_cidx(glbmodeidx,olddims,mode)
          fort(pos1)=fort(pos1)+sc*tilein(i)
        enddo
    end select
  end subroutine add_tile_to_fort

  !> \autor Patrick Ettenhuber
  !> \brief Subroutine to put a specific tile in a general matrix
  !> \date January 2013
  subroutine put_tile_in_fort(tilein,tnr,tdims,fort,dims,mode,o)
    implicit none
    !> input array with data in arbtirary order
    real(realk), intent(inout) :: fort(*)
    !> mode infortmation about how to interpret data
    integer, intent(in) :: mode
    !> the tile number in column major ordering of the tiles
    integer, intent(in) :: tnr
    !> dimension infortmation about how to interpret data
    integer, intent(in) :: dims(mode)
    !> batch information for the tiles
    integer, intent(in) :: tdims(mode)
    !> specify the order if the tile stems from a matrix with a different order
    integer,intent(in) :: o(mode)
    !> tile output
    real(realk), intent(in) :: tilein(*)
    real(realk),pointer :: dummy(:)
    integer :: i,nccblcks,nelms,k
    integer :: tmodeidx(mode),rtd(mode)
    integer :: idxintile(mode),glbidx,ro(mode),olddims(mode)
    integer :: ccels,ntimes,el,acttdim(mode),nels,fels(mode)
    integer :: pos1,pos2,ntpm(mode),glbmodeidx(mode),dummyidx(mode)
    integer :: simpleorder,bs
    bs=int(((8000.0*1000.0)/(8.0*2.0))**(1.0/float(mode)))
    !bs=5
    nels=1
    simpleorder=0
    do i=1,mode
      if(o(i)/=i)simpleorder=-1
      nels=nels*dims(i)
      ro(o(i))=i
    enddo

    
    do i=1,mode
      olddims(i) = dims(ro(i))
      ntpm(i) = dims(i)/tdims(i)
      if(mod(dims(i),tdims(i))>0)ntpm(i) = ntpm(i) + 1
    enddo
    if(mode==4)then
      if(o(1)==2.and.o(2)==1.and.o(3)==4.and.o(4)==3)simpleorder=1
      if(o(1)==1.and.o(2)==4.and.o(3)==2.and.o(4)==3)simpleorder=2
    endif
    call get_midx(tnr,tmodeidx,ntpm,mode)

    ntimes=1
    !call get_tile_dim(tdim,tnr,arr%dims,arr%tdim,arr%mode)
    do i=1,mode
      fels(i) = (tmodeidx(o(i))-1) * tdims(o(i)) + 1
      if(tmodeidx(i)*tdims(i)>dims(i))then
        acttdim(i)=mod(dims(i),tdims(i))
      else
        acttdim(i)=tdims(i)
      endif
      if(i>1)ntimes=ntimes*acttdim(i)
      rtd(o(i)) = acttdim(i)
    enddo
    

    select case(simpleorder)
      case(0)
        ccels=acttdim(1)
        do i=1,ntimes
          call get_midx(i,idxintile(2:mode),acttdim(2:mode),mode-1)
          idxintile(1)=1
          do k=1,mode
            glbmodeidx(k)=idxintile(k) +(tmodeidx(k)-1) *tdims(k)
          enddo
          pos1=get_cidx(glbmodeidx,dims,mode)
          call dcopy(ccels,tilein(1+(i-1)*ccels),1,fort(pos1),1)
        enddo
      case(1)
        call manual_2143_reordering_tile2full(bs,acttdim,dims,fels,1.0E0_realk,tilein,0.0E0_realk,fort)
      !case(2)
      !  print *,dims
      !  print *,tdims
      !  print *,o
      !  print *,olddims
      !  print *,rtd
      !  print *,ro
      !  call mem_alloc(dummy,dims(1)*dims(2)*dims(3)*dims(4))
      !  do i=1,dims(1)*dims(2)*dims(3)*dims(4)
      !    dummy(i) = fort(i)
      !  enddo
      !  call manual_1423_reordering_tile2full(bs,acttdim,dims,fels,1.0E0_realk,tilein,0.0E0_realk,dummy)
      !  print *,infpar%lg_mynum,norm2(dummy)
      !  nelms=1
      !  do i=1,mode
      !    nelms = nelms * acttdim(i)
      !  enddo
      !  do i = 1,nelms
      !    !get mode index of element in tile
      !    call get_midx(i,idxintile,acttdim,mode)
      !    if(i==nelms)print *,idxintile,tmodeidx
      !    do k=1,mode
      !      glbmodeidx(o(k))=idxintile(k) + (tmodeidx(k)-1)*tdims(k)
      !    enddo
      !    if(i==nelms)print *,glbmodeidx,olddims
      !    pos1=get_cidx(glbmodeidx,olddims,mode)
      !    fort(pos1)=tilein(i)
      !  enddo
      !  print *,infpar%lg_mynum,norm2(fort(1:dims(1)*dims(2)*dims(3)*dims(4)))
      !  do i=1,dims(1)*dims(2)*dims(3)*dims(4)
      !    if(abs(dummy(i)-fort(i)) > 1.0E-12)then
      !      call get_midx(i,dummyidx,dims,4)
      !      print *, "firs diff in ",i,"ie",dummyidx
      !      print *,"orig",fort(i)
      !      print *,"new ",dummy(i)
      !      exit
      !    endif
      !  enddo
      !  call mem_dealloc(dummy)
    case default
      print *,"default part reorder put",o
      !count elements in the current tile for loop over elements
      !identify their original position and put them in tile
      nelms=1
      do i=1,mode
        nelms = nelms * acttdim(i)
      enddo
      do i = 1,nelms
        !get mode index of element in tile
        call get_midx(i,idxintile,acttdim,mode)
        do k=1,mode
          glbmodeidx(o(k))=idxintile(k) + (tmodeidx(k)-1)*tdims(k)
        enddo
        pos1=get_cidx(glbmodeidx,olddims,mode)
        fort(pos1)=tilein(i)
      enddo
    end select
  end subroutine put_tile_in_fort

  !> \autor Patrick Ettenhuber
  !> \brief Subroutine to extract a specific tile from a general matrix
  !> \date November 2012
  subroutine extract_tile_from_fort_add(pre,fort,mode,tnr,dims,tdims,tileout,optorder)
    implicit none
    !> input array with data in arbtirary order
    real(realk), intent(in) :: fort(*),pre
    !> mode infortmation about how to interpret data
    integer, intent(in) :: mode
    !> the tile number in column major ordering of the tiles
    integer, intent(in) :: tnr
    !> dimension infortmation about how to interpret data
    integer, intent(in) :: dims(mode)
    !> batch information for the tiles
    integer, intent(in) :: tdims(mode)
    !> reorder information for the array with respect to the original array,
    !> if optorder is given, then the dimensions, tdims and tnr are with
    !reference to the tile to calculate
    integer, intent(in),optional :: optorder(mode)
    !> tile output
    real(realk), intent(out) :: tileout(*)
    integer :: i,nccblcks,nels,k
    integer :: tmodeidx(mode)
    integer :: idxintile(mode),glbidx,order(mode)
    integer :: ccels,ntimes,el,acttdim(mode),olddims(mode),nelms
    integer :: pos1,pos2,ntpm(mode),glbmodeidx(mode),rorder(mode)
    logical :: simpleorder
    do i=1,mode
      order(i)=i
    enddo
    simpleorder=.true.

    if(present(optorder))order=optorder

    do i=1,mode
      if(order(i)/=i)simpleorder=.false.
      !get the reverse order information
      rorder(order(i))=i
    enddo

    !calculate number of tiles per mode
    nels=1
    do i=1,mode
      olddims(i) = dims(rorder(i))
      ntpm(i) = dims(i)/tdims(i)
      if(mod(dims(i),tdims(i))>0)ntpm(i) = ntpm(i) + 1
      nels=nels*dims(i)
    enddo
    !get global tile index as mode-index
    call get_midx(tnr,tmodeidx,ntpm,mode)
    !get number of consecutive elements in tile form array and the actual
    !dimension of the tile
    ntimes=1
    do i=1,mode
      if(tmodeidx(i)*tdims(i)>dims(i))then
        acttdim(i)=mod(dims(i),tdims(i))
      else
        acttdim(i)=tdims(i)
      endif
      if(i>1)ntimes=ntimes*acttdim(i)
    enddo

    if(simpleorder)then

      ccels=acttdim(1)
      !loop over the remaining not-consecutive dimensions
      do i=1,ntimes
        !get the mode-index in the remaining dimensions
        call get_midx(i,idxintile(2:mode),acttdim(2:mode),mode-1)
        !get the position of the first element in the consecutive stretch
        idxintile(1)=1
        do k=1,mode
          glbmodeidx(k)=idxintile(k) +(tmodeidx(k)-1) *tdims(k)
        enddo
        pos1=get_cidx(glbmodeidx,dims,mode)
        call daxpy(ccels,pre,fort(pos1),1,tileout(1+(i-1)*ccels),1)
      enddo
    else
      print *,"default part reorder extract",order
      
      !count elements in the current tile for loop over elements
      !identify their original position and put them in tile
      nelms=1
      do i=1,mode
        nelms = nelms * acttdim(i)
      enddo
      do i = 1,nelms
        !get mode index of element in tile
        call get_midx(i,idxintile,acttdim,mode)
        !get global index of element, example order = 2 3 1 4 of new array with
        !respect to old --> element 54 3 27 8 of old goes to 3 27 54 8 of new -
        ! old with respect to new 3 1 2 4 
        do k=1,mode
          glbmodeidx(order(k))=idxintile(k) + (tmodeidx(k)-1)*tdims(k)
        enddo
        pos1=get_cidx(glbmodeidx,olddims,mode)
        tileout(i)=tileout(i)+pre*fort(pos1)
      enddo
    endif
    
  end subroutine extract_tile_from_fort_add


  !> \autor Patrick Ettenhuber
  !> \brief Subroutine to extract a specific tile from a general matrix
  !> \date November 2012
  subroutine extract_tile_from_fort(fort,mode,tnr,dims,tdims,tileout,o)
    implicit none
    !> input array with data in arbtirary order
    real(realk), intent(in) :: fort(*)
    !> mode infortmation about how to interpret data
    integer, intent(in) :: mode
    !> the tile number in column major ordering of the tiles
    integer, intent(in) :: tnr
    !> dimension infortmation about how to interpret data
    integer, intent(in) :: dims(mode)
    !> batch information for the tiles
    integer, intent(in) :: tdims(mode)
    !> reorder information for the array with respect to the original array,
    !> if optorder is given, then the dimensions, tdims and tnr are with
    !reference to the tile to calculate
    integer, intent(in) :: o(mode)
    !> tile output
    real(realk), intent(out) :: tileout(*)
    integer :: i,nccblcks,nels,k
    integer :: tmodeidx(mode)
    integer :: idxintile(mode),glbidx
    integer :: ccels,ntimes,el,acttdim(mode),olddims(mode),nelms
    integer :: pos1,pos2,ntpm(mode),glbmodeidx(mode),ro(mode),rtd(mode),fels(mode)
    integer :: simpleorder,bs,a,b,c,d
    real(realk),pointer :: dummy(:)
    bs=int(((8000.0*1000.0)/(8.0*2.0))**(1.0/float(mode)))
    !bs=5
    simpleorder=0
    do i=1,mode
      if(o(i)/=i)simpleorder=-1
    enddo


    do i=1,mode
      !get the reverse order information
      ro(o(i))=i
    enddo

    !calculate number of tiles per mode
    nels=1
    do i=1,mode
      olddims(i) = dims(ro(i))
      ntpm(i) = dims(i)/tdims(i)
      if(mod(dims(i),tdims(i))>0)ntpm(i) = ntpm(i) + 1
      nels=nels*dims(i)
    enddo
    call mem_alloc(dummy,nels)
    if(mode==4)then
      if(o(1)==2.and.o(2)==1.and.o(3)==4.and.o(4)==3)simpleorder=1
      if(o(1)==1.and.o(2)==3.and.o(3)==2.and.o(4)==4)simpleorder=2
      if(o(1)==1.and.o(2)==3.and.o(3)==4.and.o(4)==2)simpleorder=3
      if(o(1)==2.and.o(2)==1.and.o(3)==3.and.o(4)==4)simpleorder=4
      if(o(1)==1.and.o(2)==4.and.o(3)==2.and.o(4)==3)simpleorder=5
    endif
    call get_midx(tnr,tmodeidx,ntpm,mode)
    ntimes=1
    do i=1,mode
      fels(o(i)) = (tmodeidx(i)-1) * tdims(i) + 1
      if(tmodeidx(i)*tdims(i)>dims(i))then
        acttdim(i)=mod(dims(i),tdims(i))
      else
        acttdim(i)=tdims(i)
      endif
      if(i>1)ntimes=ntimes*acttdim(i)
      rtd(o(i))     = acttdim(i)
    enddo

    select case(simpleorder)
      case(0)
        ccels=acttdim(1)
        !loop over the remaining not-consecutive dimensions
        do i=1,ntimes
          !get the mode-index in the remaining dimensions
          call get_midx(i,idxintile(2:mode),acttdim(2:mode),mode-1)
          !get the position of the first element in the consecutive stretch
          idxintile(1)=1
          do k=1,mode
            glbmodeidx(k)=idxintile(k) +(tmodeidx(k)-1) *tdims(k)
          enddo
          pos1=get_cidx(glbmodeidx,dims,mode)
          call dcopy(ccels,fort(pos1),1,tileout(1+(i-1)*ccels),1)
        enddo
      case(1)
        call manual_2143_reordering_full2tile(bs,rtd,olddims,fels,1.0E0_realk,fort,0.0E0_realk,tileout)
      case(2)
        call manual_1324_reordering_full2tile(bs,rtd,olddims,fels,1.0E0_realk,fort,0.0E0_realk,tileout)
      case(3)
        call manual_1342_reordering_full2tile(bs,rtd,olddims,fels,1.0E0_realk,fort,0.0E0_realk,tileout)
      case(4)
        call manual_2134_reordering_full2tile(bs,rtd,olddims,fels,1.0E0_realk,fort,0.0E0_realk,tileout)
      case(5)
        call manual_1423_reordering_full2tile(bs,rtd,olddims,fels,1.0E0_realk,fort,0.0E0_realk,tileout)
      case default
        print *,"default part reorder extract2",o
        !count elements in the current tile for loop over elements
        !identify their original position and put them in tile
        nelms=1
        do i=1,mode
          nelms = nelms * acttdim(i)
        enddo
        do i = 1,nelms
          !get mode index of element in tile
          call get_midx(i,idxintile,acttdim,mode)
          !get global index of element, example order = 2 3 1 4 of new array with
          !respect to old --> element 54 3 27 8 of old goes to 3 27 54 8 of new -
          ! old with respect to new 3 1 2 4 
          do k=1,mode
            glbmodeidx(o(k))=idxintile(k) + (tmodeidx(k)-1)*tdims(k)
          enddo
          pos1=get_cidx(glbmodeidx,olddims,mode)
          tileout(i)=fort(pos1)
        enddo
    end select

    call mem_dealloc(dummy)

  end subroutine extract_tile_from_fort

  subroutine array_free_pdm(arr)
    implicit none
    type(array) :: arr
#ifdef VAR_LSMPI
    if(arr%init_type==MASTER_INIT.and.&
    &infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(JOB_FREE_ARR_PDM,arr)
    endif
    p_arr%free_addr_on_node(arr%addr_p_arr(infpar%lg_mynum+1))=.true.
    p_arr%arrays_in_use = p_arr%arrays_in_use - 1 
    call array_free_basic(p_arr%a(arr%addr_p_arr(infpar%lg_mynum+1))) 
    call array_nullify_pointers(arr)
#endif
  end subroutine array_free_pdm

  subroutine get_distribution_info(arr)
    implicit none
    type(array),intent(inout) :: arr
    integer :: i,ntiles2dis
#ifdef VAR_LSMPI
    arr%offset=p_arr%new_offset
    p_arr%new_offset=mod(p_arr%new_offset+arr%ntiles,infpar%lg_nodtot)
    arr%nlti=arr%ntiles/infpar%lg_nodtot
    if(mod(arr%ntiles,infpar%lg_nodtot)>&
    &mod(infpar%lg_mynum+infpar%lg_nodtot-arr%offset,infpar%lg_nodtot))arr%nlti=arr%nlti+1
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
#ifdef VAR_LSMPI
    gtnr=globtinr
    if(arr%init_type==MASTER_INIT.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(JOB_PRINT_TI_NRM,arr)
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
    elseif(infpar%lg_mynum==0)then
      write(DECinfo%output,'("LOCAL TILE NORM ON",I3,f20.15)'),dest,sqrt(norm)
    endif
#endif
  end subroutine array_tiled_pdm_print_ti_nrm

  function array_tiled_pdm_get_nrm2(arr) result(nrm)
    implicit none
    type(array), intent(in) :: arr
    real(realk) :: nrm
    integer :: i,j,should
#ifdef VAR_LSMPI
    if(infpar%lg_mynum==infpar%master.and.arr%init_type==MASTER_INIT) then
      call pdm_array_sync(JOB_GET_NRM2_TILED,arr)
    endif
    !print *,infpar%lg_mynum,"calcing nrm",arr%addr_p_arr(infpar%lg_mynum+1),arr%nlti
    !call lsmpi_barrier(infpar%lg_comm)
    nrm=0.0E0_realk
    do i=1,arr%nlti
      !call get_tile_dim(should,arr,arr%ti(i)%gt)
      !print *,infpar%lg_mynum,i,arr%ti(i)%e,arr%ti(i)%gt,should
      do j=1,arr%ti(i)%e
        nrm = nrm +(arr%ti(i)%t(j) * arr%ti(i)%t(j))
      enddo
    enddo
    !print *,infpar%lg_mynum,"do reduction",nrm
    !call lsmpi_barrier(infpar%lg_comm)
    if(arr%init_type==MASTER_INIT)call lsmpi_local_reduction(nrm,infpar%master)
    if(arr%init_type==ALL_INIT)call lsmpi_local_allreduce(nrm)
    !if(infpar%lg_mynum==infpar%master)print *,"on master",nrm
#else
    nrm = 0.0E0_realk
#endif
  end function array_tiled_pdm_get_nrm2

   subroutine change_init_type_td(arr,totype)
     implicit none
     type(array),intent(inout) :: arr
     integer,intent(in) :: totype
#ifdef VAR_LSMPI
     if(totype/=REPLICATED.and.totype/=DENSE.and.totype/=TILED_DIST.and.totype/=TILED)then
       call lsquit("ERROR(change_init_type_td): wrong type given",-1)
     endif
     if(infpar%lg_mynum==infpar%master.and.arr%init_type==MASTER_INIT) then
       call pdm_array_sync(JOB_CHANGE_INIT_TYPE,arr)
     endif
     arr%init_type=totype
#endif
   end subroutine change_init_type_td


  !> \brief direct communication routine for the accumulation of arrays,
  !> interface to the combined index routine
  !> \author Patrick Ettenhuber
  subroutine array_accumulate_tile_modeidx(arr,modidx,fort,nelms)
    implicit none
    !> input array for which a tile should be accumulated
    type(array),intent(in) ::arr
    !> input, the index of the tile in modular form and the number of elements
    integer,intent(in) :: modidx(arr%mode),nelms
    !> input the fortan array which should be transferred to the tile
    real(realk),intent(inout) :: fort(nelms)
    integer :: cidx
    cidx=get_cidx(modidx,arr%ntpm,arr%mode)
    call array_accumulate_tile(arr,cidx,fort,nelms)
  end subroutine array_accumulate_tile_modeidx
  subroutine array_accumulate_tile_combidx4(arr,globtilenr,fort,nelms)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=4) :: nelms
    !> input the fortan array which should be transferred to the tile
    real(realk),intent(inout) :: fort(nelms)
    integer(kind=ls_mpik) :: assert,dest, ierr,n
#ifdef VAR_LSMPI
    integer(kind=MPI_ADDRESS_KIND) ::offset
    assert=0
    offset=0
    n=nelms
    dest=get_residence_of_tile(globtilenr,arr)
    call lsmpi_win_lock(dest,arr%wi(globtilenr),'e')
    call lsmpi_acc(fort,nelms,1,dest,arr%wi(globtilenr))
    CALL lsmpi_win_unlock(dest, arr%wi(globtilenr))
#endif
  end subroutine array_accumulate_tile_combidx4
  subroutine array_accumulate_tile_combidx8(arr,globtilenr,fort,nelms)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=8) :: nelms
    !> input the fortan array which should be transferred to the tile
    real(realk),intent(inout) :: fort(nelms)
    integer(kind=ls_mpik) :: assert,dest, ierr,n
#ifdef VAR_LSMPI
    integer(kind=MPI_ADDRESS_KIND) ::offset
    assert=0
    offset=0
    n=nelms
    dest=get_residence_of_tile(globtilenr,arr)
    call lsmpi_win_lock(dest,arr%wi(globtilenr),'e')
    call lsmpi_acc(fort,nelms,1,dest,arr%wi(globtilenr))
    call lsmpi_win_unlock(dest,arr%wi(globtilenr))
#endif
  end subroutine array_accumulate_tile_combidx8
  subroutine array_puttile_modeidx(arr,modidx,fort,nelms)
    implicit none
    type(array),intent(in) ::arr
    integer,intent(in) :: modidx(arr%mode),nelms
    real(realk),intent(inout) :: fort(nelms)
    integer :: cidx
    cidx=get_cidx(modidx,arr%ntpm,arr%mode)
    call array_put_tile(arr,cidx,fort,nelms)
  end subroutine array_puttile_modeidx

  subroutine array_puttile_combidx8(arr,globtilenr,fort,nelms)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=8),intent(in) :: nelms
    real(realk),intent(inout) :: fort(nelms)
    integer(kind=ls_mpik) :: assert,dest, ierr,n
#ifdef VAR_LSMPI
    integer(kind=MPI_ADDRESS_KIND) ::offset
    assert=0
    offset=0
    n=nelms
    dest=get_residence_of_tile(globtilenr,arr)
    call lsmpi_win_lock(dest,arr%wi(globtilenr),'e')
    call lsmpi_put(fort,nelms,1,dest,arr%wi(globtilenr))
    call lsmpi_win_unlock(dest,arr%wi(globtilenr))
#endif
  end subroutine array_puttile_combidx8
  subroutine array_puttile_combidx4(arr,globtilenr,fort,nelms)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=4),intent(in) :: nelms
    real(realk),intent(inout) :: fort(nelms)
    integer(kind=ls_mpik) :: assert,dest, ierr,n
#ifdef VAR_LSMPI
    integer(kind=MPI_ADDRESS_KIND) ::offset
    assert=0
    offset=0
    n=nelms
    dest=get_residence_of_tile(globtilenr,arr)
    call lsmpi_win_lock(dest,arr%wi(globtilenr),'e')
    call lsmpi_put(fort,nelms,1,dest,arr%wi(globtilenr))
    call lsmpi_win_unlock(dest,arr%wi(globtilenr))
#endif
  end subroutine array_puttile_combidx4

  !interface to the array_gettile_combidx
  subroutine array_gettile_modeidx(arr,modidx,fort,nelms)
    implicit none
    type(array),intent(in) ::arr
    integer,intent(in) :: modidx(arr%mode),nelms
    real(realk),intent(inout) :: fort(nelms)
    integer :: cidx
    cidx=get_cidx(modidx,arr%ntpm,arr%mode)
    call array_get_tile(arr,cidx,fort,nelms)
  end subroutine array_gettile_modeidx

  subroutine array_gettile_combidx8(arr,globtilenr,fort,nelms)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=8),intent(in) :: nelms
    real(realk),intent(inout) :: fort(nelms)
    integer(kind=ls_mpik) :: assert,source, ierr,n
#ifdef VAR_LSMPI
    integer(kind=MPI_ADDRESS_KIND) ::offset
    assert=0
    offset=0
    n=nelms
    source=get_residence_of_tile(globtilenr,arr)
    call lsmpi_win_lock(source,arr%wi(globtilenr),'e')
    call lsmpi_get(fort,nelms,1,source,arr%wi(globtilenr))
    call lsmpi_win_unlock(source,arr%wi(globtilenr))
#endif
  end subroutine array_gettile_combidx8
  subroutine array_gettile_combidx4(arr,globtilenr,fort,nelms)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=4),intent(in) :: nelms
    real(realk),intent(inout) :: fort(nelms)
    integer(kind=ls_mpik) :: assert,source, ierr,n
#ifdef VAR_LSMPI
    integer(kind=MPI_ADDRESS_KIND) ::offset
    assert=0
    offset=0
    n=nelms
    source=get_residence_of_tile(globtilenr,arr)
    call lsmpi_win_lock(source,arr%wi(globtilenr),'s')
    call lsmpi_get(fort,nelms,1,source,arr%wi(globtilenr))
    call lsmpi_win_unlock(source,arr%wi(globtilenr))
#endif
  end subroutine array_gettile_combidx4



end module dec_pdm_module




subroutine pdm_array_slave()
  use precision
  use matrix_operations_scalapack, only: BLOCK_SIZE, SLGrid, DLEN_
  use memory_handling, only: mem_alloc,mem_dealloc
  use matrix_operations, only: mtype_scalapack, matrix_type
  use dec_typedef_module
#ifdef VAR_LSMPI
  use infpar_module
  use lsmpi_type
#endif

  ! DEC DEPENDENCIES (within deccc directory)
  ! *****************************************
  use dec_pdm_module
  use array_memory_manager, only: MASTER_INIT,&
       & arr_free_aux

   IMPLICIT NONE
   TYPE(array) :: A, B, C, D, AUX
   CHARACTER    :: T(2)
   INTEGER      :: JOB, i, j
   INTEGER      :: DESC_A(DLEN_), DESC_B(DLEN_), DESC_C(DLEN_), DESC_AF(DLEN_),dist(2)
   real(REALK)  :: AF, AB(2)
   real(REALK),pointer :: dummy(:)
   integer, pointer    :: idiag(:)
   logical :: logi
   integer, external :: numroc
   integer, pointer :: dims(:),dims2(:)
#ifdef VAR_LSMPI
   CALL PDM_ARRAY_SYNC(JOB,A,B,C,D) !Job is output
   !print *,"slave in pdm-> job is",JOB
   !print *,"array1_dims",A%dims
   !print *,"array2_dims",B%dims
   !print *,"array3_dims",C%dims

   SELECT CASE(JOB)
     CASE(JOB_INIT_ARR_TILED)
       call mem_alloc(idiag,A%mode)
       call mem_alloc(dims,A%mode)
       idiag=A%tdim
       dims =A%dims
       call arr_free_aux(A)
       A=array_init_tiled(dims,A%mode,MASTER_INIT,idiag,A%zeros) 
       call mem_dealloc(idiag)
       call mem_dealloc(dims)
     CASE(JOB_FREE_ARR_PDM)
       call array_free_pdm(A) 
     CASE(JOB_INIT_ARR_REPLICATED)
       call mem_alloc(dims,A%mode)
       dims =A%dims
       call arr_free_aux(A)
       A=array_init_replicated(dims,A%mode,MASTER_INIT) 
       call mem_dealloc(dims)
     CASE(JOB_PRINT_MEM_INFO1)
       call print_mem_per_node(DECinfo%output,.false.)
     CASE(JOB_PRINT_MEM_INFO2)
       call mem_alloc(dummy,1)
       call print_mem_per_node(DECinfo%output,.false.,dummy)
       call mem_dealloc(dummy)
     CASE(JOB_GET_NRM2_TILED)
       AF=array_tiled_pdm_get_nrm2(A)
     CASE(JOB_DATA2TILED_DIST)
       !dummy has to be allocated, otherwise seg faults might occur with 
       !some compilers
       call mem_alloc(dummy,1)
       !call cp_data2tiled(A,dummy,A%dims,A%mode,.true.)
       print *,"not necessary"
       call mem_dealloc(dummy)
     CASE(JOB_GET_TILE_SEND)
       !i,dummy and j are just dummy arguments
       call array_get_tile(A,i,dummy,j)
     CASE(JOB_PRINT_TI_NRM)
       call array_tiled_pdm_print_ti_nrm(A,0)
     CASE(JOB_SYNC_REPLICATED)
       call array_sync_replicated(A)
     CASE(JOB_GET_NORM_REPLICATED)
       AF=array_print_norm_repl(A)
     CASE(JOB_PREC_DOUBLES_PAR)
       call precondition_doubles_parallel(A,B,C,D)
     CASE(JOB_DDOT_PAR)
       AF=array_ddot_par(A,B,0)
     CASE(JOB_ADD_PAR)
       call array_add_par(A,AF,B)
     CASE(JOB_CP_ARR)
       call array_cp_tiled(A,B)
     CASE(JOB_ARRAY_ZERO)
       call array_zero_tiled_dist(A)
     CASE(JOB_GET_CC_ENERGY)
       AF = get_cc_energy_parallel(A,B,C)
     CASE(JOB_GET_FRAG_CC_ENERGY)
       !the counterpart to this buffer is in get_fragment_cc_energy
       call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
       call ls_mpi_buffer(i,infpar%master)
       call mem_alloc(dims,i)
       call ls_mpi_buffer(dims,i,infpar%master)
       call ls_mpi_buffer(j,infpar%master)
       call mem_alloc(dims2,j)
       call ls_mpi_buffer(dims2,j,infpar%master)
       call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

       AF = get_fragment_cc_energy_parallel(A,B,C,i,j,dims,dims2)

       call mem_dealloc(dims)
       call mem_dealloc(dims2)
     CASE(JOB_CHANGE_INIT_TYPE)
       call change_init_type_td(A,i)
   END SELECT
#endif
end subroutine pdm_array_slave
