!> @file
!> Operations for general arrays
!> \author Patrick Ettenhuber

module tensor_interface_module


  ! Outside DEC directory
  use memory_handling
  use precision
  use files!,only: lsopen,lsclose
  use LSTIMING!,only:lstimer
  use reorder_frontend_module
  use lspdm_tensor_operations_module
  use matrix_module


  !> Number of created arrays
  integer(kind=long) :: ArraysCreated=0
  !> Number of destroyed arrays
  integer(kind=long) :: ArraysDestroyed=0
  !> Number of created arrays
  integer(kind=long) :: CreatedPDMArrays=0
  integer(kind=long) :: DestroyedPDMArrays=0


  !> Overloaded operator for adding arrays


!> convert arrays, the idea is for a general conversion only the interface
!should be called
  interface array_convert
     module procedure array_convert_fort2arr_wrapper1,&
     &array_convert_fort2arr_wrapper2,array_convert_fort2arr_wrapper3,&
     &array_convert_fort2arr_wrapper4,array_convert_array22array,&
     &array_convert_arr2fort_wrapper1,array_convert_arr2fort_wrapper2,&
     &array_convert_arr2fort_wrapper3,array_convert_arr2fort_wrapper4
  end interface array_convert

!> print norms of array, array2 array3, array4 and fortran arrays
  interface print_norm
    module procedure print_norm_fort_wrapper1_nrm,&
                    &print_norm_fort_wrapper2_nrm,&
                    &print_norm_fort_wrapper3_nrm,&
                    &print_norm_fort_wrapper4_nrm,&
                    &array_print_norm_nrm,&
                    &array2_print_norm_nrm,&
                    &array4_print_norm_nrm,&
                    &matrix_print_norm_nrm,&
                    &print_norm_fort_wrapper1_customprint,&
                    &print_norm_fort_wrapper2_customprint,&
                    &print_norm_fort_wrapper3_customprint,&
                    &print_norm_fort_wrapper4_customprint,&
                    &array_print_norm_customprint,&
                    &array2_print_norm_customprint,&
                    &array4_print_norm_customprint
  end interface print_norm


  interface array_add
    module procedure array_add_normal, array_add_arr2fullfort,array_add_fullfort2arr
  end interface array_add

  !interface array_contract
  !  module procedure array_contract_pref
  !end interface array_contract


contains
  subroutine copy_array(arr_in,arr_out)
    implicit none
    type(array), intent(in) :: arr_in
    type(array), intent(inout) :: arr_out
    integer :: i
    arr_out%mode = arr_in%mode
    arr_out%nlti = arr_in%nlti
    arr_out%tsize = arr_in%tsize
    arr_out%itype = arr_in%itype
    arr_out%nelms = arr_in%nelms
    arr_out%ntiles = arr_in%ntiles
    arr_out%access_type = arr_in%access_type
    if(associated(arr_in%dims))call arr_set_dims(arr_out,arr_in%dims,arr_out%mode)
    if(associated(arr_in%ntpm))call arr_set_ntpm(arr_out,arr_in%ntpm,arr_out%mode)
    if(associated(arr_in%tdim))call arr_set_tdims(arr_out,arr_in%tdim,arr_out%mode)
    if(associated(arr_in%addr_p_arr))call arr_set_addr(arr_out,arr_in%addr_p_arr,size(arr_in%addr_p_arr,kind=ls_mpik))
    !arr_out%dims = arr_in%dims
    !arr_out%tdim = arr_in%tdim
    !arr_out%ntpm = arr_in%ntpm
    if(associated(arr_in%elm1))then
      call memory_allocate_array_dense(arr_out)
      arr_out%elm1=arr_in%elm1
    endif
    if(associated(arr_in%ti))then
      call memory_allocate_tiles(arr_out)
      do i=1,arr_in%nlti
        arr_out%ti(i)%t=arr_in%ti(i)%t
      enddo
    endif
  end subroutine copy_array


  ! x = x + b * y
  !> \brief add a scaled array to another array. The data may have different
  !distributions in the two arrays to be added
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_add_normal(x,b,y)
    implicit none
    !> array input, this is the result array with overwritten data
    type(array),intent(inout) :: x
    !> array to add
    type(array),intent(in) :: y
    !> scaling factor for array y
    real(realk),intent(in) :: b
    real(realk),pointer :: buffer(:)
    integer :: ti,i,nel,o(x%mode)
    do i=1,x%mode
      o(i) = i
    enddo

    select case(x%itype)
    case(DENSE)
      select case(y%itype)
      case(DENSE)
        call daxpy(int(x%nelms),b,y%elm1,1,x%elm1,1)
      case(REPLICATED)
        call daxpy(int(x%nelms),b,y%elm1,1,x%elm1,1)
        call array_sync_replicated(x)
      case(TILED_DIST)
        call mem_alloc(buffer,y%tsize)
        do ti=1,y%ntiles
          call get_tile_dim(nel,y,ti)
          call array_get_tile(y,ti,buffer,nel)
          call tile_in_fort(b,buffer,ti,y%tdim,1.0E0_realk,x%elm1,x%dims,x%mode,o)
        enddo
        call mem_dealloc(buffer)
      case default
        print *,x%itype,y%itype
        call lsquit("ERROR(array_add):not yet implemented",DECinfo%output)
      end select
    case(REPLICATED)
      select case(y%itype)
      case(DENSE)
        call daxpy(int(x%nelms),b,y%elm1,1,x%elm1,1)
        call array_sync_replicated(x)
      case(REPLICATED)
        call daxpy(int(x%nelms),b,y%elm1,1,x%elm1,1)
        call array_sync_replicated(x)
      case(TILED_DIST)
        call mem_alloc(buffer,y%tsize)
        do ti=1,y%ntiles
          call get_tile_dim(nel,y,ti)
          call array_get_tile(y,ti,buffer,nel)
          call tile_in_fort(b,buffer,ti,y%tdim,1.0E0_realk,x%elm1,x%dims,x%mode,o)
        enddo
        call mem_dealloc(buffer)
        call array_sync_replicated(x)
      case default
        print *,x%itype,y%itype
        call lsquit("ERROR(array_add):not yet implemented",DECinfo%output)
      end select
    case(TILED_DIST)
      select case(y%itype)
      case(TILED_DIST)
        call array_add_par(x,b,y)
      case default
        print *,x%itype,y%itype
        call lsquit("ERROR(array_add):not yet implemented",DECinfo%output)
      end select
    case default
      print *,x%itype,y%itype
      call lsquit("ERROR(array_add):not yet implemented",DECinfo%output)
    end select
  end subroutine array_add_normal

  ! x = x + b * y
  !> \brief add a scaled fortran-array to an array. The data may have arbitrary
  !distribution for the array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_add_fullfort2arr(arrx,b,fortarry,order,wrk,iwrk)
    implicit none
    !> full fortan arra´y, this corresponds to y
    real(realk), intent(in) :: fortarry(*)
    !> scaling factor for fortran array
    real(realk), intent(in) :: b
    !> array which is overwritten
    type(array), intent(inout) :: arrx
    !> order of the fortran array with respect to the array
    integer, intent(in),optional :: order(arrx%mode)
    !> optinally workspace can be passed, the size is defined as iwrk
    integer(kind=8), intent(in),optional :: iwrk
    real(realk), intent(inout),optional :: wrk(*)
    integer :: o(arrx%mode)
    !> check if there is enough memory to send a full tile, this will die out
    integer :: i
    real(realk) :: MemFree,tilemem
    do i=1,arrx%mode
      o(i) = i
    enddo
    if(present(order))o=order
    select case(arrx%itype)
      case(DENSE)
        if(.not.present(order))then
          call daxpy(int(arrx%nelms),b,fortarry,1,arrx%elm1,1)
        else
          call lsquit("ERROR(array_add_fullfort2arr1):not implemented",-1)
        endif
      case(TILED)
        call lsquit("ERROR(array_add_fullfort2arr):not implemented",-1)
      case(TILED_DIST)
        if(present(wrk).and.present(iwrk))then
          call add_data2tiled_intiles_explicitbuffer(arrx,b,fortarry,arrx%dims,arrx%mode,o,wrk,iwrk)
        else
          if(b==1.0E0_realk)then
            call add_data2tiled_intiles_nobuffer(arrx,fortarry,arrx%dims,arrx%mode,o)
          else
            call add_data2tiled_intiles_stackbuffer(arrx,b,fortarry,arrx%dims,arrx%mode,o)
          endif
        endif
    end select
  end subroutine array_add_fullfort2arr

  ! x = x + b * y
  !> \brief add a scaled array to a fortran- array. The data may have arbitrary
  !distribution for the array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_add_arr2fullfort(fortarrx,b,arry,order)
    implicit none
    !> full fortan arra´y, this corresponds to x and is overwritten
    real(realk), intent(inout) :: fortarrx(*)
    !> scaling factor for the array y
    real(realk), intent(in) :: b
    !> array to add --> y
    type(array), intent(in) :: arry
    !> order of the fortran array with respect to the array
    integer, intent(in),optional :: order(arry%mode)
    !> check if there is enough memory to send a full tile, this will die out
    integer :: i
    real(realk) :: MemFree,tilemem
    select case(arry%itype)
      case(DENSE)
        call daxpy(int(arry%nelms),b,arry%elm1,1,fortarrx,1)
      case(TILED)
        call lsquit("ERROR(array_add_fullfort2arr):not implemented",-1)
      case(TILED_DIST)
        if(present(order))call add_tileddata2fort(arry,b,fortarrx,arry%nelms,.true.,order)
        if(.not.present(order))call add_tileddata2fort(arry,b,fortarrx,arry%nelms,.true.)
    end select
  end subroutine array_add_arr2fullfort


  !> \brief perform a contraction of the outer indices of two arrays, r means
  !the right-most index for the first array, l means the left-most index for the
  !sectond array, p1 * left * right + p2 * res = res 
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_contract_outer_indices_rl(p1,left,right,p2,res)
    implicit none
    !> prefactors for the left*right part and the res part
    real(realk),intent(in) :: p1,p2
    !> the resulting array
    type(array),intent(inout) :: res
    !> the arrays to be multiplied
    type(array),intent(in) :: left,right
    integer :: i,m,n,k
    if(left%dims(left%mode)/=right%dims(1))then
      call lsquit("ERROR(array_contract_outer_indices_rl):wrong&
      &contraction dimensions!!",DECinfo%output)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Introduce a lot of dimensionality checks here!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,left%mode-1
      if(left%dims(i)/=res%dims(i))then
        call lsquit("ERROR(array_contract_outer_indices_rl):wrong&
        &result (left) dimensions!!",DECinfo%output)
      endif
    enddo
    do i=right%mode,2,-1
      if(right%dims(i)/=res%dims(left%mode+i-2))then
        call lsquit("ERROR(array_contract_outer_indices_rl):wrong&
        &result (right) dimensions!!",DECinfo%output)
      endif
    enddo
    if(left%itype==TILED.or.left%itype==TILED_DIST.or.right%itype==TILED&
    &.or.right%itype==TILED_DIST)call lsquit("ERROR(array_contract_outer_&
    &indices_rl):not yet implemented for tiled/PDM",DECinfo%output)
   

    !do that only if type dense --> the other routines are to come
    select case(res%itype)
      case(DENSE)
        m=1
        do i=1,left%mode-1
          m=m*left%dims(i)
        enddo
        n=1
        do i=right%mode,2,-1
          n=n*right%dims(i)
        enddo
        k=right%dims(1) 
        call dgemm('n','n',m,n,k,p1,left%elm1,m,right%elm1,k,p2,res%elm1,m)
      case(REPLICATED)
        m=1
        do i=1,left%mode-1
          m=m*left%dims(i)
        enddo
        n=1
        do i=right%mode,2,-1
          n=n*right%dims(i)
        enddo
        k=right%dims(1) 
        call dgemm('n','n',m,n,k,p1,left%elm1,m,right%elm1,k,p2,res%elm1,m)
        if(res%access_type==MASTER_ACCESS)call array_sync_replicated(res)
      case(TILED)
              call lsquit("ERROR(array_contract_outer_indices_rl):not yet implemented for tiled",DECinfo%output)
      case(TILED_DIST)
              call lsquit("ERROR(array_contract_outer_indices_rl):not yet implemented for PDM",DECinfo%output)
      case default
              call lsquit("operation for your choice of arrays not yet implemented",DECinfo%output)
    end select
    
  end subroutine array_contract_outer_indices_rl

  !> \brief perform a contraction of the outer indices of two arrays, the first l means
  !the leftt-most index for the first array, l means the left-most index for the
  !sectond array, p1 * left^T * right + p2 * res = res 
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_contract_outer_indices_ll(p1,left,right,p2,res)
    implicit none
    !> prefactors for the left*right part and the res part
    real(realk),intent(in) :: p1,p2
    !> the resulting array
    type(array),intent(inout) :: res
    !> the arrays to be multiplied
    type(array),intent(in) :: left,right
    integer :: i,m,n,k
    if(left%dims(1)/=right%dims(1))then
      call lsquit("ERROR(array_contract_outer_indices_ll):wrong&
      &dimensions!!",DECinfo%output)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Introduce a lot of dimensionality checks here!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=2,left%mode
      if(left%dims(i)/=res%dims(i-1))then
        call lsquit("ERROR(array_contract_outer_indices_ll):wrong&
        &result (left) dimensions!!",DECinfo%output)
      endif
    enddo
    do i=right%mode,2,-1
      if(right%dims(i)/=res%dims(left%mode+i-2))then
        call lsquit("ERROR(array_contract_outer_indices_ll):wrong&
        &result (right) dimensions!!",DECinfo%output)
      endif
    enddo
    if(left%itype==TILED.or.left%itype==TILED_DIST.or.right%itype==TILED&
    &.or.right%itype==TILED_DIST)call lsquit("ERROR(array_contract_outer_&
    &indices_ll):not yet implemented for tiled/PDM",DECinfo%output)

    !do that only if type dense --> the other routines are to come
    select case(res%itype)
      case(DENSE)
        m=1
        do i=left%mode,2,-1
          m=m*left%dims(i)
        enddo
        n=1
        do i=right%mode,2,-1
          n=n*right%dims(i)
        enddo
        k=right%dims(1) 
        call dgemm('t','n',m,n,k,p1,left%elm1,k,right%elm1,k,p2,res%elm1,m)
      case(REPLICATED)
        m=1
        do i=left%mode,2,-1
          m=m*left%dims(i)
        enddo
        n=1
        do i=right%mode,2,-1
          n=n*right%dims(i)
        enddo
        k=right%dims(1) 
        call dgemm('t','n',m,n,k,p1,left%elm1,k,right%elm1,k,p2,res%elm1,m)
        if(res%access_type==MASTER_ACCESS)call array_sync_replicated(res)
      case(TILED)
              call lsquit("ERROR(array_contract_outer_indices_ll):not yet implemented for tiled",DECinfo%output)
      case(TILED_DIST)
              call lsquit("ERROR(array_contract_outer_indices_ll):not yet implemented for PDM",DECinfo%output)
      case default
              call lsquit("operation for your choice of arrays not yet implemented",DECinfo%output)
    end select
    
  end subroutine array_contract_outer_indices_ll

  !> \brief perform a contraction of the outer indices of two arrays, the first l means
  !the leftt-most index for the first array, r means the r-most index for the
  !sectond array, p1 * left^T * right^T + p2 * res = res 
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_contract_outer_indices_lr(p1,left,right,p2,res)
    implicit none
    !> prefactors for the left*right part and the res part
    real(realk),intent(in) :: p1,p2
    !> the resulting array
    type(array),intent(inout) :: res
    !> the arrays to be multiplied
    type(array),intent(in) :: left,right
    integer :: i,m,n,k
    if(left%dims(1)/=right%dims(right%mode))then
      call lsquit("ERROR(array_contract_outer_indices_lr):wrong&
      &contraction dimensions!!",DECinfo%output)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Introduce a lot of dimensionality checks here!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=2,left%mode
      if(left%dims(i)/=res%dims(i-1))then
        call lsquit("ERROR(array_contract_outer_indices_lr):wrong&
        &result (left) dimensions!!",DECinfo%output)
      endif
    enddo
    do i=right%mode-1,1,-1
      if(right%dims(i)/=res%dims(left%mode+i-1))then
        call lsquit("ERROR(array_contract_outer_indices_lr):wrong&
        &result (right) dimensions!!",DECinfo%output)
      endif
    enddo
    if(left%itype==TILED.or.left%itype==TILED_DIST.or.right%itype==TILED&
    &.or.right%itype==TILED_DIST)call lsquit("ERROR(array_contract_outer_&
    &indices_lr):not yet implemented for tiled/PDM",DECinfo%output)

    !do that only if type dense --> the other routines are to come
    select case(res%itype)
      case(DENSE)
        m=1
        do i=left%mode,2,-1
          m=m*left%dims(i)
        enddo
        n=1
        do i=1,right%mode-1
          n=n*right%dims(i)
        enddo
        k=left%dims(1) 
        print *,m,n,k
        print *,left%dims
        print *,right%dims
        print *,res%dims
        call dgemm('t','t',m,n,k,p1,left%elm1,k,right%elm1,n,p2,res%elm1,m)
      case(REPLICATED)
        m=1
        do i=left%mode,2,-1
          m=m*left%dims(i)
        enddo
        n=1
        do i=1,right%mode-1
          n=n*right%dims(i)
        enddo
        k=left%dims(1) 
        print *,m,n,k
        print *,left%dims
        print *,right%dims
        print *,res%dims
        call dgemm('t','t',m,n,k,p1,left%elm1,k,right%elm1,n,p2,res%elm1,m)
        if(res%access_type==MASTER_ACCESS)call array_sync_replicated(res)
      case(TILED)
              call lsquit("ERROR(array_contract_outer_indices_lr):not yet implemented for tiled",DECinfo%output)
      case(TILED_DIST)
              call lsquit("ERROR(array_contract_outer_indices_lr):not yet implemented for PDM",DECinfo%output)
      case default
              call lsquit("operation for your choice of arrays not yet implemented",DECinfo%output)
    end select
    
  end subroutine array_contract_outer_indices_lr

  !> \brief perform a contraction of the outer indices of two arrays, the first r means
  !the right-most index for the first array, r means the r-most index for the
  !sectond array, p1 * left * right^T + p2 * res = res 
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_contract_outer_indices_rr(p1,left,right,p2,res)
    implicit none
    !> prefactors for the left*right part and the res part
    real(realk),intent(in) :: p1,p2
    !> the resulting array
    type(array),intent(inout) :: res
    !> the arrays to be multiplied
    type(array),intent(in) :: left,right
    integer :: i,m,n,k
    if(left%dims(left%mode)/=right%dims(right%mode))then
      call lsquit("ERROR(array_contract_outer_indices_rr):wrong&
      &dimensions!!",DECinfo%output)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Introduce a lot of dimensionality checks here!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,left%mode-1
      if(left%dims(i)/=res%dims(i))then
        call lsquit("ERROR(array_contract_outer_indices_rr):wrong&
        &result (left) dimensions!!",DECinfo%output)
      endif
    enddo
    do i=right%mode-1,1,-1
      if(right%dims(i)/=res%dims(left%mode+i-1))then
        call lsquit("ERROR(array_contract_outer_indices_rr):wrong&
        &result (right) dimensions!!",DECinfo%output)
      endif
    enddo
    if(left%itype==TILED.or.left%itype==TILED_DIST.or.right%itype==TILED&
    &.or.right%itype==TILED_DIST)call lsquit("ERROR(array_contract_outer_&
    &indices_rr):not yet implemented for tiled/PDM",DECinfo%output)

    select case(res%itype)
      case(DENSE)
        m=1
        do i=1,left%mode-1
          m=m*left%dims(i)
        enddo
        n=1
        do i=1,right%mode-1
          n=n*right%dims(i)
        enddo
        k=left%dims(left%mode) 
        call dgemm('n','t',m,n,k,p1,left%elm1,m,right%elm1,n,p2,res%elm1,m)
      case(REPLICATED)
        m=1
        do i=1,left%mode-1
          m=m*left%dims(i)
        enddo
        n=1
        do i=1,right%mode-1
          n=n*right%dims(i)
        enddo
        k=left%dims(left%mode) 
        call dgemm('n','t',m,n,k,p1,left%elm1,m,right%elm1,n,p2,res%elm1,m)
        if(res%access_type==MASTER_ACCESS)call array_sync_replicated(res)
      case(TILED)
              call lsquit("ERROR(array_contract_outer_indices_rr):not yet implemented for tiled",DECinfo%output)
      case(TILED_DIST)
              call lsquit("ERROR(array_contract_outer_indices_rr):not yet implemented for PDM",DECinfo%output)
      case default
              call lsquit("operation for your choice of arrays not yet implemented",DECinfo%output)
    end select
    
  end subroutine array_contract_outer_indices_rr

  !> \brief array dotprduct, the arrays may have different distributions
  !> \author Patrick Ettenhuber
  !> \date some time at the end 2012
  function array_ddot(arr1,arr2,opt_par) result(res)
    implicit none
    !> the two arrays to calculate the dotproduct from
    type(array),intent(in) :: arr1,arr2
    !> optional integer specifying on which node the result should be stored
    integer,optional,intent(in) :: opt_par
    real(realk) :: res
    integer :: dest
    real(realk), external :: ddot
    
    if(arr1%nelms/=arr2%nelms)then
      call lsquit("ERROR(array_ddot):operation not defined for arrays with&
      & different nelms",DECinfo%output)
    endif

    !get the destination of the contraction
    dest = -1
    if(arr1%access_type==MASTER_ACCESS)dest=0
    if(present(opt_par))dest=opt_par

    select case(arr1%itype)
    case(DENSE)
      select case(arr2%itype)
      case(DENSE)
        res=ddot(int(arr1%nelms),arr1%elm1,1,arr2%elm1,1)
      case default
        call lsquit("ERROR(array_ddot):operation not yet&
        & implemented",DECinfo%output)
      end select
            
    case(TILED_DIST)
      select case(arr2%itype)
      case(TILED_DIST)
        res=array_ddot_par(arr1,arr2,dest)
      case default
        call lsquit("ERROR(array_ddot):operation not yet&
        & implemented",DECinfo%output)
      end select
    case default
      call lsquit("ERROR(array_ddot):operation not yet&
      & implemented",DECinfo%output)
    end select

  end function array_ddot


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY de-/init ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> \author Patrick Ettenhuber
  !> \date January 2013
  function array_minit( dims, nmodes, local, atype, tdims) result(arr)
    !> the output array
    type(array) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer, intent(in)              :: nmodes, dims(nmodes)
    integer, intent(in),optional     :: tdims(nmodes)
    logical, intent(in),optional     :: local
    character(4),intent(in),optional :: atype
    character(4)  :: at
    integer       :: it
    logical :: loc

    ! Sanity check
    !if(arr%initialized)call lsquit("ERROR(array_minit):array already initialized",-1) 
    do i=1, nmodes
      if (dims(i) == 0) call lsquit("ERROR(array_minit): 0 dimendion not allowed",-1)
    end do

    !set defaults
    loc = .true.
    at  = 'LDAR'

    if(present(local)) loc = local
    if(present(atype)) at  = atype


#ifdef VAR_MPI
    if(loc) then
      select case(at)
      case('LDAR','REAR','REPD','TDAR','TDPD','RTAR')
        arr=array_init_standard(dims,nmodes,pdm=NO_PDM_ACCESS)
        arr%atype='LDAR'
      !case('TDAR','TDPD')
      !  arr=array_init_tiled(dims,nmodes,pdm=NO_PDM_ACCESS)
      !  arr%atype='LTAR'
      case default
        call lsquit("ERROR(array_minit): atype not known",-1)
      end select
    else
      select case(at)
      case('LDAR')
        !INITIALIZE a Local Dense ARray
        arr              = array_init_standard(dims,nmodes,pdm=MASTER_ACCESS)
        arr%atype        = 'LDAR'
      case('TDAR')
        !INITIALIZE a Tiled Distributed ARray
        it               = TILED_DIST
        if(present(tdims))then
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=MASTER_ACCESS,tdims=tdims)
        else
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=MASTER_ACCESS)
        endif
        CreatedPDMArrays = CreatedPDMArrays+1
      case('RTAR')
        !INITIALIZE a Replicated Tiled ARray (all nodes have all tiles)
        it               = TILED
        if(present(tdims))then
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=MASTER_ACCESS,tdims=tdims)
        else
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=MASTER_ACCESS)
        endif
        CreatedPDMArrays = CreatedPDMArrays+1
      case('REAR')
        !INITIALIZE a REplicated ARray
        arr              = array_init_replicated(dims,nmodes,pdm=MASTER_ACCESS)
        CreatedPDMArrays = CreatedPDMArrays+1
        arr%itype        = REPLICATED
        arr%atype        = 'REAR'
      case('TDPD')
        !INITIALIZE a Tiled Distributed Pseudo Dense array
        it               = TILED_DIST ! for array_init_tiled routine
        if(present(tdims))then
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=MASTER_ACCESS,tdims=tdims,ps_d=.true.)
        else
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=MASTER_ACCESS,ps_d=.true.)
        endif
        arr%itype        = DENSE ! back to dense after init
        CreatedPDMArrays = CreatedPDMArrays+1
      case('REPD')
        !INITIALIZE a REplicated Pseudo Dense array
        arr              = array_init_replicated(dims,nmodes,pdm=MASTER_ACCESS)
        CreatedPDMArrays = CreatedPDMArrays+1
        arr%itype        = DENSE
        arr%atype        = 'REPD'
      case default 
        call lsquit("ERROR(array_minit): atype not known",-1)
      end select
    endif
#else
    arr=array_init(dims,nmodes)
    arr%atype='LDAR'
#endif
    arr%initialized=.true.
  end function array_minit

  function array_ainit( dims, nmodes, local, atype, tdims )result(arr)
    !> the output array
    type(array) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer, intent(in)              :: nmodes, dims(nmodes)
    integer, intent(in),optional     :: tdims(nmodes)
    logical, intent(in),optional     :: local
    character(4),intent(in),optional :: atype
    character(4)  :: at
    integer       :: it
    logical :: loc
 
    ! Sanity check
    if(arr%initialized)call lsquit("ERROR(array_ainit):array already initialized",-1) 
    do i=1, nmodes
      if (dims(i) == 0) call lsquit("ERROR(array_minit): 0 dimendion not allowed",-1)
    end do
 
    !set defaults
    loc = .true.
    at  = 'LDAR'

    if(present(local)) loc = local
    if(present(atype)) at  = atype

#ifdef VAR_MPI
    if(loc) then
      select case(at)
      case('LDAR','REAR','REPD','TDAR','TDPD')
        !if local recast to a local dense array
        arr=array_init_standard(dims,nmodes,pdm=NO_PDM_ACCESS)
        arr%atype='LDAR'
      !case('TDAR','TDPD')
      !  arr=array_init_tiled(dims,nmodes,pdm=NO_PDM_ACCESS)
      !  arr%atype='LTAR'
      case default
        call lsquit("ERROR(array_minit): atype not known",-1)
      end select
    else
      select case(at)
      case('LDAR')
        !INITIALIZE a Local Dense ARray
        arr              = array_init_standard(dims,nmodes,pdm=ALL_ACCESS)
        arr%atype        = 'LDAR'
      case('TDAR')
        !INITIALIZE a Tiled Distributed ARray
        it               = TILED_DIST
        if(present(tdims))then
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=ALL_ACCESS,tdims=tdims)
        else
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=ALL_ACCESS)
        endif
        CreatedPDMArrays = CreatedPDMArrays+1
      case('REAR')
        !INITIALIZE a REplicated ARray
        arr              = array_init_replicated(dims,nmodes,pdm=ALL_ACCESS)
        CreatedPDMArrays = CreatedPDMArrays+1
        arr%itype        = REPLICATED
        arr%atype        = 'REAR'
      case('TDPD')
        !INITIALIZE a Tiled Distributed Pseudo Dense array
        it               = TILED_DIST ! for array_init_tiled routine
        if(present(tdims))then
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=ALL_ACCESS,tdims=tdims,ps_d=.true.)
        else
          arr            = array_init_tiled(dims,nmodes,at,it,pdm=ALL_ACCESS,ps_d=.true.)
        endif
        arr%itype        = DENSE ! back to dense after init
        CreatedPDMArrays = CreatedPDMArrays+1
      case('REPD')
        !INITIALIZE a REplicated Pseudo Dense array
        arr              = array_init_replicated(dims,nmodes,pdm=ALL_ACCESS)
        CreatedPDMArrays = CreatedPDMArrays+1
        arr%itype        = DENSE
        arr%atype        = 'REPD'
      case default 
        call lsquit("ERROR(array_minit): atype not known",-1)
      end select
    endif
#else
    arr=array_init_standard(dims,nmodes,NO_PDM_ACCESS)
    arr%atype='LDAR'
#endif
    arr%initialized=.true.
  end function array_ainit

  !> \author Patrick Ettenhuber
  !> \date September 2012
  !> \brief MAIN ARRAY INITIALIZATION ROUTINE
  function array_init(dims,nmodes,arr_type,pdm,tdims)result(arr)
    implicit none
    !> output array
    type(array) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer, intent(in) :: nmodes, dims(nmodes)
    !> integer specifying the type of array (a list of possible types is found
    !> at the beginning of array_memory.f90
    integer, optional :: arr_type
    !> if tiled then the size of the tile in each mode can be specified explicitly 
    integer, optional :: tdims(nmodes)
    !> specifies the type of access to the array (NO_PDM_ACCESS,MASTER_ACCESS,ALL_ACCESS)
    integer, optional :: pdm
    integer :: sel_type,pdmtype,it
    logical :: zeros_in_tiles,wcps
    !choose which kind of array

    if(arr%initialized)call lsquit("ERROR(array_init):array already initialized",-1) 

    !DEFAULTS
    it     = DENSE
    pdmtype = NO_PDM_ACCESS !NO PDM

    !OPTIONAL SPECIFICATIONS
    if(present(arr_type))    it      = arr_type
    if(present(pdm))         pdmtype = pdm

    !CHECK INPUT
    if(pdmtype>=3)call lsquit("ERROR(array_init):WRONG CHOICE IN PDMTYPE",DECinfo%output)

    !EXPERIMENTAL, THIS IS NOT RECOMMENDED!!!!!!!!!!!!!!:
    !instead of modulo dimensions in the rims of an array
    !use the same dimensions as in full tiles and fill with
    !zeros
    zeros_in_tiles = .false.
    ArraysCreated = ArraysCreated+1
    
    !select corresponding routine
    select case(it)
      case(DENSE)
        arr=array_init_standard(dims,nmodes,pdmtype)
        arr%atype = 'LDAR'
      case(REPLICATED)
        arr=array_init_replicated(dims,nmodes,pdmtype)
        arr%atype = 'REAR'
        CreatedPDMArrays = CreatedPDMArrays+1
      case(TILED)
        if(present(tdims))arr=array_init_tiled(dims,nmodes,'TIAR',it,pdmtype,tdims,zeros_in_tiles)
        if(.not.present(tdims))arr=array_init_tiled(dims,nmodes,'TIAR',it,pdmtype)
      case(TILED_DIST)
        if(present(tdims))arr=array_init_tiled(dims,nmodes,'TDAR',it,pdmtype,tdims,zeros_in_tiles)
        if(.not.present(tdims))arr=array_init_tiled(dims,nmodes,'TDAR',it,pdmtype)
        CreatedPDMArrays = CreatedPDMArrays+1
    end select
    arr%access_type   = pdmtype
    arr%itype       = it
    arr%initialized = .true.
  end function array_init



  !> \author Patrick Ettenhuber adpted from Marcin Ziolkowski
  !> \date September 2012
  !> \brief get mode index from composite index
  function array_init_standard(dims,nmodes,pdm) result(arr)
    implicit none
    integer, intent(in)   :: nmodes,dims(nmodes),pdm
    type(array)           :: arr
    logical               :: master
    integer               :: i,addr,tdimdummy(nmodes)
    integer,pointer       :: buf(:)
    integer(kind=long)    :: nelms
    integer(kind=ls_mpik) :: pc_nnodes,me

    pc_nnodes = 1
    master    = .true.
    me        = 0
#ifdef VAR_MPI
    if( lspdm_use_comm_proc ) then
      me        = infpar%pc_mynum
      pc_nnodes = infpar%pc_nodtot
      master    = (infpar%parent_comm == MPI_COMM_NULL)
    endif
#endif
    
    
    !find space in the persistent array
    p_arr%curr_addr_on_node = get_free_address(.true.)
    addr                    = p_arr%curr_addr_on_node
    p_arr%arrays_in_use     = p_arr%arrays_in_use + 1

    !SET MODE
    p_arr%a(addr)%mode      = nmodes

    !SET DIMS
    call arr_set_dims(p_arr%a(addr),dims,nmodes)

    !SET ARRAY TYPE
    p_arr%a(addr)%itype     = DENSE

    !SET INIT TYPE
    !default
    p_arr%a(addr)%access_type = NO_PDM_ACCESS
    !if one uses comm threads the following replace the access_type
    if( pdm == MASTER_ACCESS .and. lspdm_use_comm_proc )&
    & p_arr%a(addr)%access_type = MASTER_ACCESS
    if( pdm == ALL_ACCESS .and. lspdm_use_comm_proc )&
    & p_arr%a(addr)%access_type = ALL_ACCESS

    !SET IF ALLOCATED WITH COMM PROCS
    p_arr%a(addr)%allocd_w_c_p = lspdm_use_comm_proc

    !SET NELMS
    nelms=1
    do i=1,nmodes
      nelms=nelms*dims(i)
    enddo
    p_arr%a(addr)%nelms=nelms

    !put 0 in tdim, since for the replicated array it is not important
    tdimdummy=0
    call arr_set_tdims(p_arr%a(addr),tdimdummy,nmodes)

    call mem_alloc(buf,pc_nnodes)
    buf = 0

    !if master init only master has to init the addresses addresses before
    !pdm syncronization
    if(master .and. p_arr%a(addr)%access_type==MASTER_ACCESS .and. lspdm_use_comm_proc)then
      call arr_set_addr(p_arr%a(addr),buf,pc_nnodes,.true.)
#ifdef VAR_MPI
      call pdm_array_sync(infpar%pc_comm,JOB_INIT_ARR_PC,p_arr%a(addr),loc_addr=.true.)
#endif
    endif



    !if ALL_ACCESS all have to have the addresses allocated
    if(p_arr%a(addr)%access_type==ALL_ACCESS.and. lspdm_use_comm_proc)&
       &call arr_set_addr(p_arr%a(addr),buf,pc_nnodes,.true.)

    !SET THE ADDRESSES ON ALL NODES     
    buf(me+1)=addr 
#ifdef VAR_MPI
    if( lspdm_use_comm_proc )call lsmpi_allreduce(buf,pc_nnodes,infpar%pc_comm)
#endif

    call arr_set_addr(p_arr%a(addr),buf,pc_nnodes,.true.)

    call mem_dealloc(buf)

    !ALLOCATE STORAGE SPACE FOR THE ARRAY
    call memory_allocate_array_dense(p_arr%a(addr))

    !RETURN THE CURRENLY ALLOCATE ARRAY
    arr=p_arr%a(addr)

  end function array_init_standard

  

  subroutine array_free_standard(arr)
    implicit none
    type(array), intent(inout) :: arr
    integer :: me
    logical :: parent
 
    me     = 0
    parent = .true.
#ifdef VAR_MPI
    if( lspdm_use_comm_proc )then
      me = infpar%pc_mynum
      if(parent) call pdm_array_sync(infpar%pc_comm,JOB_FREE_ARR_STD,arr,loc_addr=.true.)
    endif
#endif
    p_arr%free_addr_on_node(arr%addr_loc(me+1))=.true.
    p_arr%arrays_in_use = p_arr%arrays_in_use - 1 
    call array_free_basic(p_arr%a(arr%addr_loc(me+1))) 
    call array_nullify_pointers(arr)

  end subroutine array_free_standard

  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief free a replicated matrix wich is pseudo dense for the routines,
  !> wrapper for use without mpi
  !> without mpi
  subroutine array_free_rpseudo_dense(arr,local)
    !> the array to free
    type(array) :: arr
    logical, intent(in) :: local
#ifdef VAR_MPI
    if(.not.local)then
      arr%itype=REPLICATED
    endif
#endif
    call array_free(arr)
  end subroutine array_free_rpseudo_dense

  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief free a td matrix wich is pseudo dense for the routines,
  !> wrapper for use without mpi
  !> without mpi
  subroutine array_free_tdpseudo_dense(arr,local)
    !> the array to free
    type(array) :: arr
    logical, intent(in) :: local
#ifdef VAR_MPI
    if( .not. local ) then
      arr%itype=TILED_DIST
    endif
#endif
    call array_free(arr)
  end subroutine array_free_tdpseudo_dense

  !> \brief array freeing routine, give an arbitrary array and all allocated
  !memory associated with the array will be freed
  !> \author Patrick Ettenhuber
  !> \Date probably september 2012
  subroutine array_free(arr)
    implicit none
    !> array to free
    type(array),intent(inout) :: arr

    if(.not.arr%initialized)call lsquit("ERROR(array_free):array not initialized",-1) 

    select case(arr%atype)
    case('LDAR')
      call array_free_standard(arr)
    case('TDAR','REAR','TDPD','REPD','RTAR')
      call array_free_pdm(arr)
      DestroyedPDMArrays = DestroyedPDMArrays + 1
    case('TIAR')
        call lsquit("ERROR(array_free): local tiled not maintained",-1)
    case default 
        call lsquit("ERROR(array_free): atype not known",-1)
    end select
    arr%initialized = .false. 
  end subroutine array_free


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY CONVERSION ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \brief deallocate the dense part of an array and change the type to tiled
  !distributed array. in the case of non mpi builds, this routine does nothing,
  !so that implementations still run
  !> \author Patrick Ettenhuber
  !> \date January 2013
  subroutine array_change_itype_to_td(arr,local)
    implicit none
    !> array to change the array type
    type(array),intent(inout) :: arr
    logical, intent(in) :: local
#ifdef VAR_MPI
    if( .not. local )then
      arr%itype=TILED_DIST
      if(associated(arr%elm1))then
        call arr_deallocate_dense(arr)
      endif
    endif
#else
    return
#endif
  end subroutine array_change_itype_to_td

  !> \brief change the array type to replicated, no action in case of
  !non-mpi-build
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine array_change_atype_to_rep(arr,local)
    implicit none
    !> array to change the array type
    type(array),intent(inout) :: arr
    logical, intent(in) :: local
#ifdef VAR_MPI
    if(.not.local)then
      arr%itype=REPLICATED
    endif
#else
    return
#endif
  end subroutine array_change_atype_to_rep

  !> \brief change the array type to replicated, no action in case of
  !non-mpi-build
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine array_change_atype_to_d(arr)
    implicit none
    !> array to change the array type
    type(array),intent(inout) :: arr
#ifdef VAR_MPI
    arr%itype=DENSE
#else
    return
#endif
  end subroutine array_change_atype_to_d

  !> \brief copy the tiled data of an array to its dense part, if change is
  !true, also change the %itype to DENSE, but NO deallocation of the tiled
  !distributed part
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine array_cp_tiled2dense(arr,change,order)
    implicit none
    !> array to copy data from the tiled to its dense part
    type(array),intent(inout) :: arr
    !> logical to specify whether to change the atype
    logical :: change
    !> if order is given the dense part will be reordered with respect to the
    !tiled distributed part
    integer,intent(in),optional:: order(arr%mode)
    logical :: pdm
    pdm=.false.
    if(arr%itype/=DENSE)then
      if(.not.associated(arr%elm1))then
        call memory_allocate_array_dense(arr)
      else
        call lsquit("ERROR(array_cp_tiled2dense):dense is already allocated,&
        & please make sure you are not doing someting stupid",DECinfo%output)
      endif
      if(arr%access_type>0)pdm=.true.
      if(.not.present(order))call cp_tileddata2fort(arr,arr%elm1,arr%nelms,pdm)
      if(present(order))call cp_tileddata2fort(arr,arr%elm1,arr%nelms,pdm,order)
      if(change)arr%itype=DENSE
    endif
  end subroutine array_cp_tiled2dense

  !> \brief copy the dense part of an array to its tiled distributed part and
  !delallocate the dense part afterwards, if desired change the atype. no
  !action if non-mpi build.
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine array_mv_dense2tiled(arr,change)
    implicit none
    !> array to copy dense to tiled and deallocate dense part
    type(array),intent(inout) :: arr
    logical, intent(in) :: change
    logical :: pdm
#ifdef VAR_MPI
    pdm=.false.
    if(.not.associated(arr%elm1))then
      call lsquit("ERROR(array_cp_dense2tiled):dense is NOT allocated,&
      & please make sure you are not doing someting stupid",DECinfo%output)
    endif
    if(arr%access_type>0)pdm=.true.
    if(change)arr%itype=TILED_DIST
    call array_convert_fort2arr(arr%elm1,arr,arr%nelms)
    call arr_deallocate_dense(arr)
#else
    return
#endif
  end subroutine array_mv_dense2tiled


  !\brief all the following wrappers are necessary to use the conversion routine
  !in an interface for different shapes of the fortran array  
  subroutine array_convert_arr2fort_wrapper1(arr,fort,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(inout) :: fort(:)
    integer, intent(in), optional :: order(arr%mode)
    call array_convert_arr2fort(arr,fort,arr%nelms,order=order)
  end subroutine array_convert_arr2fort_wrapper1
  subroutine array_convert_arr2fort_wrapper2(arr,fort,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(inout) :: fort(:,:)
    integer, intent(in), optional :: order(arr%mode)
    call array_convert_arr2fort(arr,fort,arr%nelms,order=order)
  end subroutine array_convert_arr2fort_wrapper2
  subroutine array_convert_arr2fort_wrapper3(arr,fort,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(inout) :: fort(:,:,:)
    integer, intent(in), optional :: order(arr%mode)
    call array_convert_arr2fort(arr,fort,arr%nelms,order=order)
  end subroutine array_convert_arr2fort_wrapper3
  subroutine array_convert_arr2fort_wrapper4(arr,fort,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(inout) :: fort(:,:,:,:)
    integer, intent(in), optional :: order(arr%mode)
    call array_convert_arr2fort(arr,fort,arr%nelms,order=order)
  end subroutine array_convert_arr2fort_wrapper4

  !\brief all the following wrappers are necessary to use the conversion routine
  !in an interface for different shapes of the fortran array  
  subroutine array_convert_fort2arr_wrapper1(fortarr,arr,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(in) :: fortarr(arr%nelms)
    integer, intent(in),optional :: order(arr%mode)
    call array_convert_fort2arr(fortarr,arr,arr%nelms,order=order)
  end subroutine array_convert_fort2arr_wrapper1
  subroutine array_convert_fort2arr_wrapper2(fortarr,arr,order)
    implicit none
    real(realk), intent(in) :: fortarr(:,:)
    type(array), intent(inout) :: arr
    integer, intent(in),optional :: order(arr%mode)
    call array_convert_fort2arr(fortarr,arr,arr%nelms,order=order)
  end subroutine array_convert_fort2arr_wrapper2
  subroutine array_convert_fort2arr_wrapper3(fortarr,arr,order)
    implicit none
    real(realk), intent(in) :: fortarr(:,:,:)
    type(array), intent(inout) :: arr
    integer, intent(in),optional :: order(arr%mode)
    call array_convert_fort2arr(fortarr,arr,arr%nelms,order=order)
  end subroutine array_convert_fort2arr_wrapper3
  subroutine array_convert_fort2arr_wrapper4(fortarr,arr,order)
    implicit none
    real(realk), intent(in) :: fortarr(:,:,:,:)
    type(array), intent(inout) :: arr
    integer, intent(in),optional :: order(arr%mode)
    call array_convert_fort2arr(fortarr,arr,arr%nelms,order=order)
  end subroutine array_convert_fort2arr_wrapper4


  !> \brief put data of a fortan array into an arbitrary array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_convert_fort2arr(fortarr,arr,nelms,order)
    implicit none
    !> the fortran array with the data
    real(realk), intent(in) :: fortarr(*)
    !> the array which should contain the data after the operation
    type(array), intent(inout) :: arr
    !> number of elements to copy from the fortan array to the array
    integer(kind=8), intent(in) :: nelms
    !> if the array should have a different ordering than the fortran array,
    ! this can be specified with order
    integer, intent(in),optional :: order(arr%mode)
    !> checkmem is outdated
    real(realk) :: tilemem,MemFree
    integer :: i,o(arr%mode),fullfortdims(arr%mode)
    real(realk) :: nrm
    do  i=1,arr%mode
      o(i)=i
    enddo
    if(present(order))o=order
    do  i=1,arr%mode
      fullfortdims(o(i))=arr%dims(i)
    enddo
    
    if(nelms/=arr%nelms)call lsquit("ERROR(array_convert_fort2arr):array&
    &dimensions are not the same",-1)
    select case(arr%itype)
      case(DENSE)
        call dcopy(int(nelms),fortarr,1,arr%elm1,1)
      case(TILED)
        call cp_data2tiled_lowmem(arr,fortarr,arr%dims,arr%mode)
      case(TILED_DIST)
        !if enough memory is available the lower one should be faster      
        if(arr%access_type==ALL_ACCESS)then
          do i=1,arr%nlti
            call tile_from_fort(1.0E0_realk,fortarr,fullfortdims,arr%mode,&
                               &0.0E0_realk,arr%ti(i)%t,arr%ti(i)%gt,arr%tdim,o)
          enddo
        else
          !call cp_data2tiled_lowmem(arr,fortarr,arr%dims,arr%mode)
          call cp_data2tiled_intiles(arr,fortarr,arr%dims,arr%mode,o)
        endif
      case default
        call lsquit("ERROR(array_convert_fort2arr) the array type is not implemented",-1)
    end select
  end subroutine array_convert_fort2arr
  !> \brief change the init type for a fortan array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine change_access_type(arr,totype)
    implicit none
    !> array to chage the init type
    type(array),intent(inout) :: arr
    !> type to change it to
    integer,intent(in) :: totype
    if(arr%itype==TILED_DIST.or.arr%itype==REPLICATED.or.&
         &totype==TILED_DIST.or.totype==REPLICATED)then
      call change_access_type_td(arr,totype)
    else
      call lsquit("ERROR(change_access_type): what you want to do is not implemented",-1)
    endif
  end subroutine change_access_type

  
  !> \brief put data of an arbitrary array into a basic fortan type array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_convert_arr2fort(arr,fort,nelms,order)
    implicit none
    !> array with the data at the beginning
    type(array), intent(inout) :: arr
    !> fortan array to contain the data in the end
    real(realk), intent(inout) :: fort(*)
    !> number of elements to convert, must be the same as elements in the array
    integer(kind=8), intent(in) :: nelms
    !> if the fortan array has a different order than the array this specifies
    !the reordering
    integer, intent(in), optional :: order(arr%mode)
    if(nelms/=arr%nelms)call lsquit("ERROR(array_convert_arr2fort):array&
    &dimensions are not the same",DECinfo%output)
    select case(arr%itype)
      case(DENSE)
        call dcopy(int(nelms),arr%elm1,1,fort,1)
      case(TILED)
        if(present(order))call cp_tileddata2fort(arr,fort,nelms,.false.,order)
        if(.not.present(order))call cp_tileddata2fort(arr,fort,nelms,.false.)
      case(TILED_DIST)
        if(present(order))call cp_tileddata2fort(arr,fort,nelms,.true.,order)
        if(.not.present(order))call cp_tileddata2fort(arr,fort,nelms,.true.)
    end select
  end subroutine array_convert_arr2fort


  
  !> \brief convert an old array2 structure to an array structure
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_convert_array22array(arraytwo,arr)
    implicit none
    !> array2 input
    type(array2),intent(in) :: arraytwo
    !> array output
    type(array), intent(inout) :: arr
    integer(kind=8) :: nel
    if(arr%mode/=2)then
      call lsquit("ERROR(array_convert_array22array):wrong mode in arr&
      &input)",DECinfo%output)
    endif
    nel = int(arraytwo%dims(1)*arraytwo%dims(2),kind=8)
    if(arr%nelms/=nel)then
      call lsquit("ERROR(array_convert_array22array):number of elements in arr&
      &input)",DECinfo%output)
    endif
    call array_convert_fort2arr(arraytwo%val,arr,nel)
  end subroutine array_convert_array22array


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY UTILITIES !!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine array_cp_data(from_arr,to_arr)
    implicit none
    type(array),intent(in) :: from_arr
    type(array),intent(inout) :: to_arr
    integer :: i
    real(realk) :: tilemem,MemFree
    if(from_arr%nelms/=to_arr%nelms)then
      call lsquit("ERROR(array_cp_data):arrays need the same number of& 
      & elements",DECinfo%output)
    endif
    
    select case(from_arr%itype)
      case(DENSE)
        select case(to_arr%itype)
          case(DENSE)
            call dcopy(int(from_arr%nelms),from_arr%elm1,1,to_arr%elm1,1)
          case(REPLICATED)
            call dcopy(int(from_arr%nelms),from_arr%elm1,1,to_arr%elm1,1)
            call array_sync_replicated(to_arr)
          case(TILED_DIST)
            call cp_data2tiled_intiles(to_arr,from_arr%elm1,from_arr%dims,from_arr%mode)
          case default
            call lsquit("ERROR(array_cp_data):operation not yet&
            & implemented",DECinfo%output)
        end select
      case(REPLICATED)
        select case(to_arr%itype)
          case(DENSE)
            call dcopy(int(from_arr%nelms),from_arr%elm1,1,to_arr%elm1,1)
          case(REPLICATED)
            call dcopy(int(from_arr%nelms),from_arr%elm1,1,to_arr%elm1,1)
            call array_sync_replicated(to_arr)
          case(TILED_DIST)
            call cp_data2tiled_intiles(to_arr,from_arr%elm1,from_arr%dims,from_arr%mode)
          case default
            call lsquit("ERROR(array_cp_data):operation not yet&
            & implemented",DECinfo%output)
        end select
      case(TILED_DIST)
        select case(to_arr%itype)
          case(DENSE)
            call cp_tileddata2fort(from_arr,to_arr%elm1,from_arr%nelms,.true.)
          case(TILED_DIST)
            call array_cp_tiled(from_arr,to_arr)
          case default
            call lsquit("ERROR(array_cp_data):operation not yet&
            & implemented",DECinfo%output)
        end select
      case default
        call lsquit("ERROR(array_cp_data):operation not yet&
        & implemented",DECinfo%output)
    end select

  end subroutine array_cp_data

  subroutine array_zero(zeroed)
    implicit none
    type(array) :: zeroed
    integer :: i
    
    select case(zeroed%itype)
      case(DENSE)
        zeroed%elm1=0.0E0_realk
      case(REPLICATED)
        zeroed%elm1=0.0E0_realk
        call array_sync_replicated(zeroed)
      case(TILED)
        if (zeroed%atype=='RTAR') then
          call array_zero_tiled_dist(zeroed)
        else
          do i=1,zeroed%ntiles
            zeroed%ti(i)%t=0.0E0_realk
          enddo
        end if
      case(TILED_DIST)
        call array_zero_tiled_dist(zeroed)
      case default
        call lsquit("ERROR(array_zero):not yet implemented",-1)
    end select

  end subroutine array_zero

  subroutine array_print_tile_norm(arr,globtinr,nrm,returnsquared)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtinr
    real(realk),intent(inout),optional::nrm
    logical,intent(in),optional :: returnsquared
    real(realk)::norm
    integer :: loctinr,i,j,on
    logical :: squareback
    squareback=.false.
    if(present(returnsquared))squareback=returnsquared
    norm=0.0d0
    if(globtinr>arr%ntiles)then
      call lsquit("ERROR(array_print_tile_norm):tile does not exist",DECinfo%output)
    endif 
    select case(arr%itype)
      case(DENSE)
        print *,"WARNING INVALID OPTION(array_print_tile_norm):no tiles in dense array"     
        do i=1,arr%nelms
          norm=norm+arr%elm1(i)*arr%elm1(i)
        enddo
        on = 1
      case(TILED)
        do j=1,arr%ti(globtinr)%e
          norm=norm + arr%ti(globtinr)%t(j) * arr%ti(globtinr)%t(j)
        enddo
        on = 1
      case(TILED_DIST)
        call array_tiled_pdm_print_ti_nrm(arr,globtinr,on,norm)
    end select
    if(.not.squareback)norm=sqrt(norm)
    if(present(nrm))nrm=norm
    if(.not.present(nrm))then
    if(.not.squareback)write(DECinfo%output,'("LOCAL TILE NORM ON",I3,f20.15)')on,norm
    if(squareback)write(DECinfo%output,'("LOCAL TILE NORM^2 ON",I3,f20.15)')on,norm
    endif
  end subroutine array_print_tile_norm

  subroutine array_print_norm_nrm(arr,nrm,returnsquared)
    implicit none
    real(realk),intent(inout),optional :: nrm
    type(array),intent(in) :: arr
    logical,intent(in),optional :: returnsquared
    real(realk)::norm
    integer(kind=8) :: i,j
    logical :: squareback
    squareback=.false.
    if(present(returnsquared))squareback=returnsquared
    norm=0.0d0
    select case(arr%itype)
      case(DENSE)
        do i=1,arr%nelms
          norm=norm+arr%elm1(i)*arr%elm1(i)
        enddo
      case(REPLICATED)
        norm = array_print_norm_repl(arr)
      case(TILED)
        do i=1,arr%nlti
          do j=1,arr%ti(i)%e
            !ISNAN is gfort intrinsic, so this options should be commented out in
            !general
            !if(ISNAN(arr%ti(i)%t(j)))then
            !  write(*,'("NaN detected in norm_t, tile:",I5," element: ",I5)')i,j
            !  stop 1
            !endif
            norm=norm + arr%ti(i)%t(j) * arr%ti(i)%t(j)
          enddo
        enddo
      case(TILED_DIST)
        norm=array_tiled_pdm_get_nrm2(arr)
    end select
    if(.not.squareback)norm = sqrt(norm)
    if(.not.squareback.and..not.present(nrm))print *,"NORM:",norm
    if(squareback.and..not.present(nrm))print *,"NORM^2:",norm
    if(present(nrm))nrm=norm
  end subroutine array_print_norm_nrm
  subroutine array_print_norm_customprint(arr,msg,returnsquared)
    implicit none
    character(ARR_MSG_LEN),intent(in) :: msg
    type(array),intent(in) :: arr
    logical,intent(in),optional :: returnsquared
    real(realk)::norm
    integer(kind=8) :: i,j
    logical :: squareback
    squareback=.false.
    if(present(returnsquared))squareback=returnsquared
    norm=0.0d0
    select case(arr%itype)
      case(DENSE)
        do i=1,arr%nelms
          norm=norm+arr%elm1(i)*arr%elm1(i)
        enddo
      case(TILED)
        do i=1,arr%nlti
          do j=1,arr%ti(i)%e
            norm=norm + arr%ti(i)%t(j) * arr%ti(i)%t(j)
          enddo
        enddo
      case(TILED_DIST)
        norm=array_tiled_pdm_get_nrm2(arr)
    end select
    if(.not.squareback)norm = sqrt(norm)
    print *,msg,norm
  end subroutine array_print_norm_customprint


  !> \brief Master routine for getting memory information in different shapes
  !> \autonr Patrick Ettenhuber
  subroutine array_print_mem_info(output,print_all_nodes,allaccess,reducetocheck)
    implicit none
    !> integer controling the output, if there should be any
    integer, intent(in) :: output
    !> optional logigal whether all nodes should print their information
    logical, intent(in), optional :: print_all_nodes
    !> optional logical if all nodes access this routine at the same time
    logical, intent(in), optional :: allaccess
    !> optional logigal stating that the informaion should be gathered on master
    integer, intent(inout), optional :: reducetocheck
    logical :: alln,red,master
    integer :: nnod
    real(realk),pointer :: red_info(:)
    alln = .false.
    red = .false.
    master=.true.
#ifdef VAR_MPI
    if(infpar%lg_mynum/=0)master=.false.
    nnod=infpar%lg_nodtot
#endif
    if(present(print_all_nodes))alln=print_all_nodes
    if(present(reducetocheck))then
      red=.true.
      if(red.and.master)then
#ifdef VAR_MPI
      call mem_alloc(red_info,nnod*8)
#else
      call mem_alloc(red_info,8)
#endif
      endif
    endif
    if(alln)then
      if(present(allaccess).and.red)call print_mem_per_node(output,allaccess,red_info)
      if(.not.present(allaccess).and.red)call print_mem_per_node(output,.false.,red_info)
      if(present(allaccess).and..not.red)call print_mem_per_node(output,allaccess)
      if(.not.present(allaccess).and..not.red)call print_mem_per_node(output,.false.)
    else
      call array_print_memory_currents(output)
    endif
    if(red.and.master)then
      if(abs(red_info(1))<1.0E-11_realk)then
        reducetocheck=0
      else
        reducetocheck=1
      endif
      call mem_dealloc(red_info)
    endif

  end subroutine array_print_mem_info

  subroutine print_norm_fort_wrapper1_nrm(fort,nelms,nrm,square)
    implicit none
    real(realk),intent(in) :: fort(:)
    integer(kind=8),intent(in) ::  nelms
    real(realk),intent(out),optional :: nrm
    logical,intent(in),optional :: square
    if(present(nrm))call print_norm_fort_nrm(fort,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(fort,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(fort,nelms)
  end subroutine print_norm_fort_wrapper1_nrm
  subroutine print_norm_fort_wrapper2_nrm(fort,nelms,nrm,square)
    implicit none
    real(realk),intent(in) :: fort(:,:)
    integer(kind=8),intent(in) ::  nelms
    real(realk),intent(out),optional :: nrm
    logical,intent(in),optional :: square
    if(present(nrm))call print_norm_fort_nrm(fort,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(fort,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(fort,nelms)
  end subroutine print_norm_fort_wrapper2_nrm
  subroutine print_norm_fort_wrapper3_nrm(fort,nelms,nrm,square)
    implicit none
    real(realk),intent(in) :: fort(:,:,:)
    integer(kind=8),intent(in) ::  nelms
    real(realk),intent(out),optional :: nrm
    logical,intent(in),optional :: square
    if(present(nrm))call print_norm_fort_nrm(fort,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(fort,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(fort,nelms)
  end subroutine print_norm_fort_wrapper3_nrm
  subroutine print_norm_fort_wrapper4_nrm(fort,nelms,nrm,square)
    implicit none
    real(realk),intent(in) :: fort(:,:,:,:)
    integer(kind=8),intent(in) ::  nelms
    real(realk),intent(out),optional :: nrm
    logical,intent(in),optional :: square
    if(present(nrm))call print_norm_fort_nrm(fort,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(fort,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(fort,nelms)
  end subroutine print_norm_fort_wrapper4_nrm
  subroutine print_norm_fort_wrapper1_customprint(fort,nelms,msg,square)
    implicit none
    real(realk),intent(in) :: fort(:)
    integer(kind=8),intent(in) ::  nelms
    character(ARR_MSG_LEN),intent(in) :: msg
    logical,intent(in),optional :: square
    if(.not.present(square))call print_norm_fort_customprint(fort,nelms,msg)
    if(present(square))call print_norm_fort_customprint(fort,nelms,msg,square)
  end subroutine print_norm_fort_wrapper1_customprint
  subroutine print_norm_fort_wrapper2_customprint(fort,nelms,msg,square)
    implicit none
    real(realk),intent(in) :: fort(:,:)
    integer(kind=8),intent(in) ::  nelms
    character(ARR_MSG_LEN),intent(in) :: msg
    logical,intent(in),optional :: square
    if(.not.present(square))call print_norm_fort_customprint(fort,nelms,msg)
    if(present(square))call print_norm_fort_customprint(fort,nelms,msg,square)
  end subroutine print_norm_fort_wrapper2_customprint
  subroutine print_norm_fort_wrapper3_customprint(fort,nelms,msg,square)
    implicit none
    real(realk),intent(in) :: fort(:,:,:)
    integer(kind=8),intent(in) ::  nelms
    character(ARR_MSG_LEN),intent(in) :: msg
    logical,intent(in),optional :: square
    if(.not.present(square))call print_norm_fort_customprint(fort,nelms,msg)
    if(present(square))call print_norm_fort_customprint(fort,nelms,msg,square)
  end subroutine print_norm_fort_wrapper3_customprint
  subroutine print_norm_fort_wrapper4_customprint(fort,nelms,msg,square)
    implicit none
    real(realk),intent(in) :: fort(:,:,:,:)
    integer(kind=8),intent(in) ::  nelms
    character(ARR_MSG_LEN),intent(in) :: msg
    logical,intent(in),optional :: square
    if(.not.present(square))call print_norm_fort_customprint(fort,nelms,msg)
    if(present(square))call print_norm_fort_customprint(fort,nelms,msg,square)
  end subroutine print_norm_fort_wrapper4_customprint

  subroutine array2_print_norm_nrm(arrtwo,nrm,square)
    implicit none
    type(array2),intent(in) :: arrtwo
    real(realk),intent(inout), optional :: nrm
    integer(kind=8) :: nelms
    logical,intent(in),optional :: square
    nelms = int(arrtwo%dims(1)*arrtwo%dims(2),kind=8)
    if(present(nrm).and..not.present(square))call print_norm_fort_nrm(arrtwo%val,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(arrtwo%val,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(arrtwo%val,nelms)
  end subroutine array2_print_norm_nrm
  subroutine array2_print_norm_customprint(arrtwo,msg,square)
    implicit none
    type(array2),intent(in) :: arrtwo
    character(ARR_MSG_LEN), intent(in) :: msg
    integer(kind=8) :: nelms
    logical,intent(in),optional :: square
    nelms = int(arrtwo%dims(1)*arrtwo%dims(2),kind=8)
    if(.not.present(square))call print_norm_fort_customprint(arrtwo%val,nelms,msg)
    if(present(square))call print_norm_fort_customprint(arrtwo%val,nelms,msg,square)
  end subroutine array2_print_norm_customprint
  subroutine array4_print_norm_nrm(arrf,nrm,square)
    implicit none
    type(array4),intent(in) :: arrf
    real(realk),intent(inout), optional :: nrm
    logical,intent(in),optional :: square
    integer(kind=8) :: nelms
    nelms = int(arrf%dims(1)*arrf%dims(2)*arrf%dims(3)*arrf%dims(4),kind=8)
    if(present(nrm).and..not.present(square))call print_norm_fort_nrm(arrf%val,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(arrf%val,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(arrf%val,nelms)
  end subroutine array4_print_norm_nrm
  subroutine array4_print_norm_customprint(arrf,msg,square)
    implicit none
    type(array4),intent(in) :: arrf
    character(ARR_MSG_LEN),intent(in):: msg
    logical,intent(in),optional :: square
    integer(kind=8) :: nelms
    nelms = int(arrf%dims(1)*arrf%dims(2)*arrf%dims(3)*arrf%dims(4),kind=8)
    if(.not.present(square))call print_norm_fort_customprint(arrf%val,nelms,msg)
    if(present(square))call print_norm_fort_customprint(arrf%val,nelms,msg,square)
  end subroutine array4_print_norm_customprint
  subroutine matrix_print_norm_nrm(mat,nrm,square)
    implicit none
    type(matrix),intent(in) :: mat
    real(realk),intent(inout), optional :: nrm
    logical,intent(in),optional :: square
    integer(kind=8) :: nelms
    nelms = int(mat%nrow*mat%ncol,kind=8)
    if(present(nrm).and..not.present(square))call print_norm_fort_nrm(mat%elms,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(mat%elms,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(mat%elms,nelms)
  end subroutine matrix_print_norm_nrm

  subroutine print_norm_fort_nrm(fort,nelms,nrm,returnsquared)
    implicit none
    real(realk),intent(in) :: fort(*)
    integer(kind=8),intent(in) ::  nelms
    real(realk),intent(out),optional :: nrm
    logical,intent(in),optional :: returnsquared
    integer(kind=8) :: i
    real(realk) :: norm
    logical :: squareback
    squareback=.false.
    if(present(returnsquared))squareback=returnsquared
    norm=0.0E0_realk
    do i=1,nelms
      norm = norm + fort(i) * fort(i)
    enddo
    if(.not.squareback)norm = sqrt(norm)
    if(.not.present(nrm).and..not.squareback)print *,"NORM:",norm
    if(.not.present(nrm).and.squareback)print *,"NORM^2:",norm
    if(present(nrm))nrm=norm
  end subroutine print_norm_fort_nrm
  subroutine print_norm_fort_customprint(fort,nelms,string,returnsquared)
    implicit none
    real(realk),intent(in) :: fort(*)
    integer(kind=8),intent(in) ::  nelms
    character(ARR_MSG_LEN),intent(in) :: string
    logical,intent(in),optional :: returnsquared
    integer(kind=8) :: i
    real(realk) :: norm
    logical :: squareback
    squareback=.false.
    if(present(returnsquared))squareback=returnsquared
    norm=0.0E0_realk
    do i=1,nelms
      norm = norm + fort(i) * fort(i)
    enddo
    if(.not.squareback)norm = sqrt(norm)
    print *,string,norm
  end subroutine print_norm_fort_customprint



  subroutine array_scale(arr,sc)
    implicit none
    type(array) :: arr
    real(realk) :: sc
    
    select case(arr%itype)
    case(DENSE)
      call dscal(int(arr%nelms),sc,arr%elm1,1)
    case(REPLICATED)
      call dscal(int(arr%nelms),sc,arr%elm1,1)
      call array_sync_replicated(arr)
    case(TILED_DIST)
      call array_scale_td(arr,sc)
    case default
      call lsquit("ERROR(array_scale):not yet implemented",DECinfo%output)
    end select
  end subroutine array_scale
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY TESTCASES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> \brief Test the array structure 
  !> \author Patrick Ettenhuber
  !> \date September 2012
  subroutine test_array_struct()
    implicit none

    type(array) :: test,test2
    real(realk),pointer :: dummy1(:),tileget(:),dummy2(:)
    real(realk),pointer :: tileget2(:)
    real(realk) :: normher,ref,ref2,ref3
    integer(kind=long) :: testint
    logical :: master
    integer :: no,nv,nb,na,i,j,succ,to_get_from,ti,midx(4)
    integer(kind=ls_mpik) :: sender, recver, nnod, rnk, me
    character(len=7) :: teststatus
    character(ARR_MSG_LEN) :: msg
    master = .true.
    nnod   = 1_ls_mpik
    me     = 0
#ifdef VAR_MPI
    me = infpar%lg_mynum
    if(me /= 0) then
      master =.false.
    endif
    nnod = infpar%lg_nodtot
    if(nnod < 3) print*,"WARNING(test_array_struct): not enough MPI processes to test all features"
#endif
    nb =  21
    nv =  18
    no =  12
    na =  7
 

#ifdef VAR_MPI
    if(master)then
      write(DECinfo%output,*)"TESTING PDM TILED ARRAY ALLOCATIONS"
      write(DECinfo%output,'(" Using",f8.3," GB of mem for the testarray")')&
      &(nv*no*(nv+nb)*8.0E0_realk)/(1024.E0_realk*1024.E0_realk*1024.E0_realk)
      testint=2

      call mem_alloc(dummy1,nb*na*nv*no)
      call mem_alloc(dummy2,nb*na*nv*no)
      call random_number(dummy1)
      call print_norm(dummy1,int(nb*na*nv*no,kind=8),ref)
      write(DECinfo%output,'("REFERENCE NORM:",f19.12)')ref
     

      !DIFFERENT ALLOCATION AND DEALLOCATION STEPS
      write (DECinfo%output,*) ""
      write (DECinfo%output,*) ""
      write (DECinfo%output,*) "TESTING SIMPLE ARRAY FUNCTIONS - MASTER DIRECTED"
      write (DECinfo%output,*) ""
      write (DECinfo%output,*)"ALLOC-DEALLOC TESTS"
      print *,"alloc dealloc tests"
      teststatus="SUCCESS"
      test=array_init([nv,na,nv,nb],4,TILED_DIST,MASTER_ACCESS,[nv,no-1,1,2])
      test2=array_init([na,nb,nv,no],4,TILED_DIST,MASTER_ACCESS,[nv,no-1,1,2])
      call array_free(test2)
      test2=array_init([no,no+1,no-1,no+1],4,TILED_DIST,MASTER_ACCESS,[no,no-1,nv,nb])
      call array_free(test)
      call array_free(test2)
      call array_print_mem_info(DECinfo%output,.true.,.false.,succ)
      if(succ/=0)teststatus=" FAILED"
      test2=array_init([nb,no,nv,no+1],4,TILED_DIST,MASTER_ACCESS,[nb,2,3,4])
      write (DECinfo%output,'(" ALLOC-DEALLOC TESTS: ",A7)')teststatus  
      print *,"DIFFERENT ALLOCATION AND DEALLOCATION STEPS: ",teststatus

      !ALLOCATING A FULL MATRIX AND PUT IT TO DISTRIBUTED MEMORY
      !check for errors via norm
      write(DECinfo%output,*)""
      write(DECinfo%output,*)""
      teststatus="SUCCESS"
      test=array_init([nb,na,nv,no],4,TILED_DIST,MASTER_ACCESS,[nb,na-1,3,no/2])
      write (DECinfo%output,*) "CONVERT PREVIOUS ARRAY TO PDM TILED" 
      call array_convert(dummy1,test,[1,2,3,4])
      call print_norm(test,normher)
      write(DECinfo%output,'("NORM OF PDM ARRAY  : ",f20.15)')normher
      if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      write (DECinfo%output,'("CNVRT: NORM, TEST STATUS:",f19.10," : ",A7)')normher,teststatus
      print *,"ALLOCATING A FULL MATRIX AND PUT IT TO DISTRIBUTED MEMORY: ",normher,teststatus


      !GET A TILE OF A PDM ARRAY
      !calculate how many elements are in the desired tile, and allocate the
      !respective amount of memory in a fortran array
      write(DECinfo%output,*)""
      write(DECinfo%output,*)""
      write(DECinfo%output,*)"TESTING MPI_GET"
      do ti = 1, test%ntiles
        to_get_from = get_residence_of_tile(ti,test)
        if(to_get_from /= me .or. nnod==1)then
         testint = ti
         exit
        endif 
      enddo
      call get_tile_dim(j,test,testint)
      print *,"trying to get",testint," with size", j
      call mem_alloc(tileget,j)
      call random_number(tileget)
      !initiatilize with some weird number here 10 and after get compare the
      !norms of the tile and the local fortran array
      teststatus="SUCCESS"
      if(nnod>1) call array_print_tile_norm(test,ti,ref)
      write(DECinfo%output,'("NORM OF TILE IN ARRAY   : ",f20.15)')ref
      if(nnod>1) call print_norm(tileget,int(j,kind=8),normher)
      write(DECinfo%output,'("NORM OF FORT BEFORE GET : ",f20.15)')normher
      if(nnod>1) call array_get_tile(test,ti,tileget,j)
      if(nnod>1) call print_norm(tileget,int(j,kind=8),normher)
      write(DECinfo%output,'("NORM OF FORT AFTER GET  : ",f20.15)')normher
      if(nnod>1)then
        if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      else
        print *,"GET A TILE OF A PDM ARRAY has been skipped"
      endif
      write (DECinfo%output,'("GET: NORM, TEST STATUS:  ",f20.15," : ",A7)')normher,teststatus
      call mem_dealloc(tileget)
      print *,"GET A TILE OF A PDM ARRAY: ",normher,teststatus

      write(DECinfo%output,*)""
      write(DECinfo%output,*)""
      write(DECinfo%output,*)"TESTING MPI_PUT"
      teststatus="SUCCESS"
      do ti = test%ntiles,1, -1
        to_get_from = get_residence_of_tile(ti,test)
        if(to_get_from /= me .or. nnod==1)then
         testint = ti
         exit
        endif 
      enddo
      call get_tile_dim(j,test,testint)
      call mem_alloc(tileget,j)
      call random_number(tileget)
      !initiatilize with some weird number here 10 and after get compare the
      !norms of the tile and the local fortran array
      if(nnod>1) call array_print_tile_norm(test,ti,normher)
      write(DECinfo%output,'("NORM OF TILE BEFORE PUT : ",f20.15)')normher
      if(nnod>1) call print_norm(tileget,int(j,kind=8),ref)
      write(DECinfo%output,'("NORM OF FORT TO PUT     : ",f20.15)')ref
      if(nnod>1) call array_put_tile(test,ti,tileget,j)
      if(nnod>1) call array_print_tile_norm(test,ti,normher)
      write(DECinfo%output,'("NORM OF TILE AFTER PUT  : ",f20.15)')normher
      if(nnod>1)then 
        if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      else
        print *,"PUT A TILE OF A PDM ARRAY has been skipped"
      endif
      write (DECinfo%output,'("PUT: NORM, TEST STATUS:  ",f20.15," : ",A7)')normher,teststatus
      print *,"TESTING MPI_PUT",normher,teststatus
     

      !GET A TILE FROM YOURSELF AS CHECK
      !testint=1
      !call mem_dealloc(tileget)
      !j=get_tileinfo_nels(test,testint)
      !call mem_alloc(tileget,j)
      !call array_get_tile(test,1,tileget,j)
      !call print_norm(tileget,j)
      !call array_print_tile_norm(test,1)
      !call sleep(4)


      !CHECK MPI PUT AND ACCUMULATE IN THE SAME WAY
      !add three to the current tile and print the norm
      write(DECinfo%output,*)""
      write(DECinfo%output,*)""
      write(DECinfo%output,*)"TESTING MPI_ACCUMULATE"
      teststatus="SUCCESS"
      if(nnod>1)then
        do i=1,j
          tileget(i)=tileget(i)+3.0E0_realk
        enddo
        call print_norm(tileget,int(j,kind=8),ref)
      endif
      write(DECinfo%output,'("NORM LOCAL ACCUMULATION : ",f20.15)')ref
      !initialize the local tile with 3 and accumulate it --> compare norm
      tileget=3.0E0_realk
      if(nnod>1) call print_norm(tileget,int(j,kind=8),normher)
      write(DECinfo%output,'("NORM OF FORT TO ADD:      ",f20.15)')normher
      if(nnod>1) call array_accumulate_tile(test,ti,tileget,j)
      if(nnod>1) call array_print_tile_norm(test,ti,normher)
      write(DECinfo%output,'("NORM REMOTE ACCUMULATION: ",f20.15)')normher
      !use the tile with three in it, print its norm put and compare norms
      if(nnod>1)then 
        if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      else
        print *,"ACCUMULATE A TILE OF A PDM ARRAY has been skipped"
      endif
      write (DECinfo%output,'("ACC: NORM, TEST STATUS:  ",f20.15," : ",A7)')normher,teststatus
      print *,"TESTING MPI_ACCUMULATE: ",normher,teststatus


      call array_free(test)
      test=array_init([nb,na,nv,no],4,TILED_DIST,MASTER_ACCESS,[0,0,0,0])
      write(DECinfo%output,*)""
      write(DECinfo%output,*)""
      write(DECinfo%output,*)"TESTING CONVERSION TO FORT"
      call array_convert(dummy1,test)
      teststatus="SUCCESS"
      call print_norm(dummy1,int(nb*na*nv*no,kind=8),ref)
      write(DECinfo%output,'("NORM OF DENSE ARRAY:      ",f20.15)')ref
      call print_norm(dummy1,int(no*nv*na*nb,kind=8),normher)
      write(DECinfo%output,'("NORM OF PDM ARRAY :       ",f20.15)')normher
      dummy2=1.0E13_realk
      call array_convert(test,dummy2)
      call print_norm(dummy2,int(no*nv*na*nb,kind=8),normher)
      write(DECinfo%output,'("NORM OF CONTRACTED ARRAY: ",f20.15)')normher
      if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      write (DECinfo%output,'("CTR: NORM, TEST STATUS:  ",f20.15," : ",A7)')normher,teststatus
      teststatus="SUCCESS"
      do i=1,no*nv*na*nb
        if(abs(dummy1(i)-dummy2(i))>1.0E-12)then
          print *,"element",i,dummy1(i),dummy2(i)
          teststatus=" FAILED"
        endif
      enddo
      write (DECinfo%output,'("ORDER: TEST STATUS:                              ",A7)')teststatus
      print *,"TESTING CONVERSION TO FORT: ",teststatus



      call mem_dealloc(dummy1)
      call mem_dealloc(dummy2)
      call mem_dealloc(tileget)
      !call array_print_mem_info(DECinfo%output,.true.,.false.)
      call array_free(test)
      call array_free(test2)

      teststatus="SUCCESS"
      call array_print_mem_info(DECinfo%output,.true.,.false.,succ)
      !call array_print_mem_info(DECinfo%output,.true.,.false.)
      if(succ/=0)teststatus=" FAILED"
       write (DECinfo%output,'("FIRST HALF ALLOCATION: ",A7)')teststatus  
    endif
    !get the slaves into this routine
    if(master)then
      print *,"MASTER GETTING SLAVES"
      call ls_mpibcast(ARRAYTEST,infpar%master,infpar%lg_comm)
      write (DECinfo%output,*)""
      write (DECinfo%output,*)""
      write (DECinfo%output,*)"TESTING PARALLEL ACCESS TO THE SAME ROUTINES"
      write (DECinfo%output,*)""
    else
      print *,"SLAVE ARRIVED",infpar%lg_mynum
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!ALL OF THE SLAVES WILL BE HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    testint=2

    !initialize a matrix
    teststatus="SUCCESS"
    if(master) write (DECinfo%output,*)"ALLOC-DEALLOC TESTS"
    test=array_init([nb,nb+2,nb+3,nb+4],4,TILED_DIST,ALL_ACCESS,[nb,nb+2,40,2])
    test2=array_init([no+3,no+2,no+1,no],4,TILED_DIST,ALL_ACCESS,[no,40,40,10])
    call array_free(test)
    call array_free(test2)
    test2=array_init([nb,na,nv,no],4,TILED_DIST,ALL_ACCESS,[nb,nv-1,1,2])
    call array_free(test2)
    call array_print_mem_info(DECinfo%output,.true.,.true.,succ)
    if(succ/=0)teststatus=" FAILED"
    test2=array_init([nb,no,nv,no+1],4,TILED_DIST,ALL_ACCESS,[nb,2,3,4])
    if(master) write (DECinfo%output,'(" ALLOC-DEALLOC TESTS: ",A7)')teststatus  
    if(master) write (DECinfo%output,*)"DONE -- NOW COMMUNICATION"
    if(master) write(DECinfo%output,*)""
    if(master) write(DECinfo%output,*)""
    print *,"ALL-INIT ALLOC-DEALLOC TESTS",teststatus
    !call lsmpi_barrier(infpar%lg_comm)

    !IF MY RANK IS NNOD-1, PUT A MATRIX CONTAINING 10 the first tile not on the
    !current rank
    teststatus="SUCCESS"
    rnk = nnod - 1
    do ti = test%ntiles,1, -1
      to_get_from = get_residence_of_tile(ti,test)
      if(to_get_from /= nnod-1 .and. to_get_from/=nnod-2)then
       testint = ti
       exit
      endif 
    enddo


    if((infpar%lg_mynum==rnk.or.master).and. nnod > 2)then

      recver=rnk

      if(.not.master)then
        call get_tile_dim(j,test2,ti)
        call mem_alloc(tileget,j)
        tileget = 1.0E1_realk
        call array_put_tile(test2,ti,tileget,j)
        call print_norm(tileget,int(j,kind=8),normher)
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        call mem_dealloc(tileget)

      else
        call ls_mpisendrecv(ref,infpar%lg_comm,recver,infpar%master)
        write(DECinfo%output,'("NORM PARALLEL 3LPN: ",f20.15)')ref
      endif

    else
      if(master) print*,"WARNING: skipping test NORM PARALLEL 3LPN, not enough nodes"
    endif


    !BEFORE rank NNOD - 2  CAN GET THE TILE
    rnk = nnod - 2

    
    call lsmpi_barrier(infpar%lg_comm)
    if((infpar%lg_mynum==rnk.or.master).and.nnod>2)then
      recver=rnk
      if(.not.master)then
        call get_tile_dim(j,test2,ti)
        call mem_alloc(tileget,j)
        call array_get_tile(test2,ti,tileget,j)
        call print_norm(tileget,int(j,kind=8),normher)
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        do i=1,j
          tileget(i) = tileget(i) + 2.4E0_realk
        enddo
        call print_norm(tileget,int(j,kind=8),normher)
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        tileget = 2.4E0_realk
        call get_midx(ti,midx,test2%ntpm,test2%mode)
        call array_accumulate_tile(test2,midx,tileget,j)
        call mem_dealloc(tileget)
      else
        teststatus="SUCCESS"
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        write(DECinfo%output,'("NORM PARALLEL 2LGN: ",f20.15)')normher
        if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
        write (DECinfo%output,'("PUT-GET: NORM, TEST STATUS: ",f19.10," : ",A7)')normher,teststatus

        call ls_mpisendrecv(ref,infpar%lg_comm,recver,infpar%master)
        write(DECinfo%output,'("NORM PARALLEL 2LAC: ",f20.15)')ref
      endif
    endif

    !BE CAREFUL ABOUT WHETER THE INFORMATION IS ALREADY TRANSMITTED --> AT
    !CRITICAL POINTS INSERT BARRIER STATEMENTS TO SYNCHONIZE THE NODES 
    call lsmpi_barrier(infpar%lg_comm)

    if(nnod>2)call array_print_tile_norm(test2,ti,normher)
    call array_free(test2)
    if(master)then
       teststatus="SUCCESS"
       write(DECinfo%output,'("NORM PARALLEL WORK: ",f20.15)')normher
       if(nnod>2)then
         if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
       else
         print *,"AS TEST WAS SKIPPED WE DO NOT CHECK FOR THE RESULT"
       endif
       write (DECinfo%output,'("ACC2    : NORM, TEST STATUS: ",f19.10," : ",A7)')ref,teststatus
       print *,"ACC2    : NORM, TEST STATUS:",teststatus
    endif
    if(master) write (DECinfo%output,*)""

  
    !test extracting a tile with a different ordering than the dense matrix, put
    !that into pdm, get these tiles on each node and put them in reversed
    !reordering back into the full array, check norms and order
    teststatus="SUCCESS"
    call lsmpi_barrier(infpar%lg_comm)
    test2=array_init([no-4,nv+3,nv/7,no],4,TILED_DIST,ALL_ACCESS,[no-4,nv+3,5,2])
    test=array_init([nv/7,nv+3,no,no-4],4,TILED_DIST,ALL_ACCESS)
    call memory_allocate_array_dense(test)
    call random_number(test%elm1)
    call lsmpi_allreduce(test%elm1,test%nelms,infpar%lg_comm)
    if(infpar%lg_mynum==0)then
      write (msg,*)"local test norm master"
      call print_norm(test%elm1,test%nelms,msg)
    endif
    if(infpar%lg_nodtot>1)then
      rnk = 1
    else
      rnk = 0
    endif
    if(me==rnk)then
      write (msg,*)"local test norm slave"
      call print_norm(test%elm1,test%nelms,msg)
    endif
    call print_norm(test%elm1,test%nelms,ref)
    call array_convert(test%elm1,test2,[4,2,1,3])
    call lsmpi_barrier(infpar%lg_comm)
    call print_norm(test2,normher)
    print *,"convert",ref,normher
    call array_mv_dense2tiled(test,.false.)
    call memory_allocate_array_dense(test2)
    call lsmpi_barrier(infpar%lg_comm)
    test2%elm1=0.0E0_realk
    do i=1,test2%ntiles
      call get_tile_dim(j,test2,i)
      call mem_alloc(tileget,j)
      call array_get_tile(test2,i,tileget,j)
      call tile_in_fort(1.0E0_realk,tileget,i,test2%tdim,0.0E0_realk,&
                        &test2%elm1,test%dims,4,[3,2,4,1])
      call mem_dealloc(tileget)
    enddo
    call lsmpi_barrier(infpar%lg_comm)
    call array_cp_tiled2dense(test,.true.)
    call lsmpi_barrier(infpar%lg_comm)
    if(infpar%lg_mynum==0)then
      write (msg,*)"local test 2 norm master"
      call print_norm(test2%elm1,test2%nelms,msg)
    endif
    if(infpar%lg_nodtot>1)then
      rnk = 1
    else
      rnk = 0
    endif
    if(infpar%lg_mynum==rnk)then
      write (msg,*)"local test 2 norm slave"
      call print_norm(test2%elm1,test2%nelms,msg)
    endif
    call print_norm(test2%elm1,test2%nelms,normher)
    do i=1,test%nelms
      if(abs(test%elm1(i)-test2%elm1(i))>1.0E-12)then
        teststatus=" FAILED"
      endif
    enddo
    call arr_deallocate_dense(test)
    call arr_deallocate_dense(test2)
    test%itype=TILED_DIST
    test2%itype=TILED_DIST
    call array_free(test)
    call array_free(test2)
    if(master)then
       write(DECinfo%output,'("PDM REORDERINGS: ",f20.15)')normher
       if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
       write (DECinfo%output,'("PDMR    : NORM, TEST STATUS: ",f19.10," : ",A7)')ref,teststatus
    endif
    if(master) write (DECinfo%output,*)""
#endif

  end subroutine test_array_struct 

end module tensor_interface_module


#ifdef VAR_MPI
subroutine get_slaves_to_array_test()
  use precision
  use tensor_interface_module,only:test_array_struct
  implicit none
  call test_array_struct()
end subroutine get_slaves_to_array_test
#endif
