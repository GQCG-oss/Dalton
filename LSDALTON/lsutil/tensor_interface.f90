!> @file
!> Operations for general arrays
!> \author Patrick Ettenhuber

module tensor_interface_module


  ! Outside DEC directory
  use memory_handling
  use precision
  use files!,only: lsopen,lsclose
  use LSTIMING!,only:lstimer
  use manual_reorderings_module 
  use lspdm_tensor_operations_module


  ! DEC DEPENDENCIES (within deccc directory)
  ! *****************************************
  !use dec_fragment_utils
  !use array4_simple_operations!,only:array_reorder_4d
  !use array3_simple_operations!,only:array_reorder_3d
  !use array_memory_manager
  !use dec_pdm_module


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
    arr_out%atype = arr_in%atype
    arr_out%nelms = arr_in%nelms
    arr_out%ntiles = arr_in%ntiles
    arr_out%init_type = arr_in%init_type
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


!  subroutine array_contract_pref(pre1,A,B,cmA,cmB,NM2C,pre2,C,order)
!    implicit none
!    real(realk), intent(in) :: pre1,pre2
!    type(array), intent(in) :: A, B
!    type(array), intent(inout) :: C
!    integer, intent(in) :: NM2C
!    integer, intent(in) :: cmA(NM2C),cmB(NM2C)
!    integer, optional, intent(in) :: order(C%mode)
!    integer :: i,j,car,cad,cbr,cbd,k
!    integer :: cmA_o(NM2C)!,cmB_o(NM2C)
!    integer, pointer :: ar(:),ad(:),br(:),bd(:)
!
!    !Categorize contractions
!    car=0;cad=0;cbr=0;cbd=0
!    do i=1,NM2C
!      do j=i+1,NM2C
!        if(cmA(i)<cmA(j))then
!          cmA_o(i)=cmA(j)
!        else
!          cmA_o(i)=cmA(i)
!        endif
!      enddo
!    enddo
!    print *,cmA_o,cmA
!!    print *,cmB_o,cmB
!    print *,"not yet further implemented"
!    stop 0
!
!
!  end subroutine array_contract_pref

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
    integer :: ti,nel
    select case(x%atype)
    case(DENSE)
      select case(y%atype)
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
          call add_tile_to_fort(buffer,ti,y%tdim,x%elm1,x%dims,x%mode,b)
        enddo
        call mem_dealloc(buffer)
      case default
        print *,x%atype,y%atype
        call lsquit("ERROR(array_add):not yet implemented",DECinfo%output)
      end select
    case(REPLICATED)
      select case(y%atype)
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
          call add_tile_to_fort(buffer,ti,y%tdim,x%elm1,x%dims,x%mode,b)
        enddo
        call mem_dealloc(buffer)
        call array_sync_replicated(x)
      case default
        print *,x%atype,y%atype
        call lsquit("ERROR(array_add):not yet implemented",DECinfo%output)
      end select
    case(TILED_DIST)
      select case(y%atype)
      case(TILED_DIST)
        call array_add_par(x,b,y)
      case default
        print *,x%atype,y%atype
        call lsquit("ERROR(array_add):not yet implemented",DECinfo%output)
      end select
    case default
      print *,x%atype,y%atype
      call lsquit("ERROR(array_add):not yet implemented",DECinfo%output)
    end select
  end subroutine array_add_normal

  ! x = x + b * y
  !> \brief add a scaled fortran-array to an array. The data may have arbitrary
  !distribution for the array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_add_fullfort2arr(arrx,b,fortarry,order)
    implicit none
    !> full fortan arra´y, this corresponds to y
    real(realk), intent(in) :: fortarry(*)
    !> scaling factor for fortran array
    real(realk), intent(in) :: b
    !> array which is overwritten
    type(array), intent(inout) :: arrx
    !> order of the fortran array with respect to the array
    integer, intent(in),optional :: order(arrx%mode)
    !> check if there is enough memory to send a full tile, this will die out
    integer :: i
    real(realk) :: MemFree,tilemem
    select case(arrx%atype)
      case(DENSE)
        call daxpy(int(arrx%nelms),b,fortarry,1,arrx%elm1,1)
      case(TILED)
        call lsquit("ERROR(array_add_fullfort2arr):not implemented",-1)
      case(TILED_DIST)
        !if(present(order))call lsquit("ERROR(array_add_fullfort2arr:not implemented",-1)
        !call add_data2tiled_lowmem(arrx,b,fortarry,arrx%dims,arrx%mode)
        if(present(order))call add_data2tiled_intiles(arrx,b,fortarry,arrx%dims,arrx%mode,order)
        if(.not.present(order))call add_data2tiled_intiles(arrx,b,fortarry,arrx%dims,arrx%mode)
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
    select case(arry%atype)
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
    if(left%atype==TILED.or.left%atype==TILED_DIST.or.right%atype==TILED&
    &.or.right%atype==TILED_DIST)call lsquit("ERROR(array_contract_outer_&
    &indices_rl):not yet implemented for tiled/PDM",DECinfo%output)

    !do that only if type dense --> the other routines are to come
    select case(res%atype)
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
        if(res%init_type==MASTER_INIT)call array_sync_replicated(res)
      case(TILED)
              call lsquit("ERROR(array_contract_outer_indices_rl):not yet implemented for tiled",DECinfo%output)
      case(TILED_DIST)
              call lsquit("ERROR(array_contract_outer_indices_rl):not yet implemented for PDM",DECinfo%output)
      case(SCALAPACK)
              call lsquit("scalapack for arrays not yet implemented",DECinfo%output)
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
    if(left%atype==TILED.or.left%atype==TILED_DIST.or.right%atype==TILED&
    &.or.right%atype==TILED_DIST)call lsquit("ERROR(array_contract_outer_&
    &indices_ll):not yet implemented for tiled/PDM",DECinfo%output)

    !do that only if type dense --> the other routines are to come
    select case(res%atype)
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
        if(res%init_type==MASTER_INIT)call array_sync_replicated(res)
      case(TILED)
              call lsquit("ERROR(array_contract_outer_indices_ll):not yet implemented for tiled",DECinfo%output)
      case(TILED_DIST)
              call lsquit("ERROR(array_contract_outer_indices_ll):not yet implemented for PDM",DECinfo%output)
      case(SCALAPACK)
              call lsquit("scalapack for arrays not yet implemented",DECinfo%output)
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
    if(left%atype==TILED.or.left%atype==TILED_DIST.or.right%atype==TILED&
    &.or.right%atype==TILED_DIST)call lsquit("ERROR(array_contract_outer_&
    &indices_lr):not yet implemented for tiled/PDM",DECinfo%output)

    !do that only if type dense --> the other routines are to come
    select case(res%atype)
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
        if(res%init_type==MASTER_INIT)call array_sync_replicated(res)
      case(TILED)
              call lsquit("ERROR(array_contract_outer_indices_lr):not yet implemented for tiled",DECinfo%output)
      case(TILED_DIST)
              call lsquit("ERROR(array_contract_outer_indices_lr):not yet implemented for PDM",DECinfo%output)
      case(SCALAPACK)
              call lsquit("scalapack for arrays not yet implemented",DECinfo%output)
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
    if(left%atype==TILED.or.left%atype==TILED_DIST.or.right%atype==TILED&
    &.or.right%atype==TILED_DIST)call lsquit("ERROR(array_contract_outer_&
    &indices_rr):not yet implemented for tiled/PDM",DECinfo%output)

    select case(res%atype)
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
        if(res%init_type==MASTER_INIT)call array_sync_replicated(res)
      case(TILED)
              call lsquit("ERROR(array_contract_outer_indices_rr):not yet implemented for tiled",DECinfo%output)
      case(TILED_DIST)
              call lsquit("ERROR(array_contract_outer_indices_rr):not yet implemented for PDM",DECinfo%output)
      case(SCALAPACK)
              call lsquit("scalapack for arrays not yet implemented",DECinfo%output)
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
    if(arr1%init_type==MASTER_INIT)dest=0
    if(present(opt_par))dest=opt_par

    select case(arr1%atype)
    case(DENSE)
      select case(arr2%atype)
      case(DENSE)
        res=ddot(int(arr1%nelms),arr1%elm1,1,arr2%elm1,1)
      case default
        call lsquit("ERROR(array_ddot):operation not yet&
        & implemented",DECinfo%output)
      end select
            
    case(TILED_DIST)
      select case(arr2%atype)
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
  !> \brief initializes a tiled dist matrix from master, wrapper for use without mpi
  function array_minit_td(dims,nmodes)result(arr)
    !> the output array
    type(array) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer, intent(in) :: nmodes, dims(nmodes)
#ifdef VAR_LSMPI
    arr=array_init_tiled(dims,nmodes,MASTER_INIT)
    CreatedPDMArrays = CreatedPDMArrays+1
    arr%atype=TILED_DIST
#else
    arr=array_init(dims,nmodes)
#endif
  end function array_minit_td
  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief initializes a tiled dist matrix wich has local dense memory 
  !>  allocated from master, wrapper for use without mpi
  function array_minit_tdpseudo_dense(dims,nmodes)result(arr)
    !> the output array
    type(array) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer, intent(in) :: nmodes, dims(nmodes)
#ifdef VAR_LSMPI
    arr=array_init_tiled(dims,nmodes,MASTER_INIT)
    CreatedPDMArrays = CreatedPDMArrays+1
    call memory_allocate_array_dense(arr)
    arr%atype=DENSE
#else
    arr=array_init(dims,nmodes)
#endif
  end function array_minit_tdpseudo_dense
  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief initializes a rep matrix wich has type dense, wrapper for use
  !> without mpi
  function array_minit_rpseudo_dense(dims,nmodes)result(arr)
    !> the output array
    type(array) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer, intent(in) :: nmodes, dims(nmodes)
#ifdef VAR_LSMPI
    arr=array_init_replicated(dims,nmodes,MASTER_INIT)
    CreatedPDMArrays = CreatedPDMArrays+1
    arr%atype=DENSE
#else
    arr=array_init(dims,nmodes)
#endif
  end function array_minit_rpseudo_dense

  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief free a replicated matrix wich is pseudo dense for the routines,
  !> wrapper for use without mpi
  !> without mpi
  subroutine array_free_rpseudo_dense(arr)
    !> the array to free
    type(array) :: arr
#ifdef VAR_LSMPI
    arr%atype=REPLICATED
#endif
    call array_free(arr)
  end subroutine array_free_rpseudo_dense

  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief free a td matrix wich is pseudo dense for the routines,
  !> wrapper for use without mpi
  !> without mpi
  subroutine array_free_tdpseudo_dense(arr)
    !> the array to free
    type(array) :: arr
#ifdef VAR_LSMPI
    arr%atype=TILED_DIST
#endif
    call array_free(arr)
  end subroutine array_free_tdpseudo_dense

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
    !> specifies the type of access to the array (NO_PDM,MASTER_INIT,ALL_INIT)
    integer, optional :: pdm
    integer :: sel_type,pdmtype,atype
    logical :: zeros_in_tiles
    !choose which kind of array
    if(present(arr_type))then
      atype=arr_type
    else
      atype=DENSE
    endif

    !EXPERIMENTAL, THIS IS NOT RECOMMENDED!!!!!!!!!!!!!!:
    !instead of modulo dimensions in the rims of an array
    !use the same dimensions as in full tiles and fill with
    !zeros
    zeros_in_tiles = .false.
    ArraysCreated = ArraysCreated+1
    
    pdmtype=NO_PDM !NO PDM
    if(present(pdm))then
      pdmtype=pdm
      if(pdmtype>=3)call lsquit("ERROR(array_init):WRONG CHOICE IN PDMTYPE",DECinfo%output)
    endif
    !select corresponding routine
    select case(atype)
      case(DENSE)
        arr=array_init_standard(dims,nmodes)
      case(REPLICATED)
        arr=array_init_replicated(dims,nmodes,pdmtype)
        CreatedPDMArrays = CreatedPDMArrays+1
      case(TILED)
        if(present(tdims))arr=array_init_tiled(dims,nmodes,pdmtype,tdims,zeros_in_tiles)
        if(.not.present(tdims))arr=array_init_tiled(dims,nmodes,pdmtype)
      case(TILED_DIST)
        if(present(tdims))arr=array_init_tiled(dims,nmodes,pdmtype,tdims,zeros_in_tiles)
        if(.not.present(tdims))arr=array_init_tiled(dims,nmodes,pdmtype)
        CreatedPDMArrays = CreatedPDMArrays+1
      !case(SCALAPACK)
      !  call lsquit("scalapack for arrays not yet implemented",DECinfo%output)
      !  arr = array_init_scalapack(dims)
      !  CreatedPDMArrays = CreatedPDMArrays+1
    end select
    arr%init_type=pdmtype
    arr%atype=atype
  end function array_init



  !> \author Patrick Ettenhuber adpted from Marcin Ziolkowski
  !> \date September 2012
  !> \brief get mode index from composite index
  function array_init_standard(dims,nmodes) result(arr)
    implicit none
    integer, intent(in) :: nmodes,dims(nmodes)
    type(array) :: arr
    logical :: file_exist
    integer :: ierr,i
    integer(kind=long) :: nelms
    real(realk) :: t0,t1

    arr%mode=nmodes
    call arr_set_dims(arr,dims,nmodes)
    nelms=1
    do i=1,nmodes
      nelms=nelms*dims(i)
    enddo

    arr%nelms=nelms
    ! allocate the dense array 
    call memory_allocate_array_dense(arr)
    call ls_dzero(arr%elm1,size(arr%elm1))


    ! By default: Do not create files when arrays are stored in memory.
    arr%funit=0
    arr%filename = 'NoFilename'
    arr%address_counter=0
    arr%storing_type=0
    arr%nelements=0
    nullify(arr%address)


  end function array_init_standard

  !> \brief array freeing routine, give an arbitrary array and all allocated
  !memory associated with the array will be freed
  !> \author Patrick Ettenhuber
  !> \Date probably september 2012
  subroutine array_free(arr)
    implicit none
    !> array to free
    type(array),intent(inout) :: arr
    real(realk) :: t0,t1
    ArraysDestroyed = ArraysDestroyed + 1
    select case(arr%atype)
      case(DENSE)
        call array_free_basic(arr)
      case(REPLICATED)
        call array_free_pdm(arr)
        DestroyedPDMArrays = DestroyedPDMArrays + 1
      case(TILED)
        call array_free_basic(arr)
      case(TILED_DIST)
        call array_free_pdm(arr)
        DestroyedPDMArrays = DestroyedPDMArrays + 1
      case(SCALAPACK)
        call lsquit("SCALAPACK FOR ARRAY NOT YET IMPLEMENTED",DECinfo%output)
    end select
    !call print_memory_currents(DECinfo%output)
  end subroutine array_free


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY CONVERSION ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \brief deallocate the dense part of an array and change the type to tiled
  !distributed array. in the case of non mpi builds, this routine does nothing,
  !so that implementations still run
  !> \author Patrick Ettenhuber
  !> \date January 2013
  subroutine array_change_atype_to_td(arr)
    implicit none
    !> array to change the array type
    type(array),intent(inout) :: arr
#ifdef VAR_LSMPI
    arr%atype=TILED_DIST
    if(associated(arr%elm1))then
      call memory_deallocate_array_dense(arr)
    endif 
#else
    return
#endif
  end subroutine array_change_atype_to_td

  !> \brief change the array type to replicated, no action in case of
  !non-mpi-build
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine array_change_atype_to_rep(arr)
    implicit none
    !> array to change the array type
    type(array),intent(inout) :: arr
#ifdef VAR_LSMPI
    arr%atype=REPLICATED
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
#ifdef VAR_LSMPI
    arr%atype=DENSE
#else
    return
#endif
  end subroutine array_change_atype_to_d

  !> \brief copy the tiled data of an array to its dense part, if change is
  !true, also change the %atype to DENSE, but NO deallocation of the tiled
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
    if(arr%atype/=DENSE)then
      if(.not.associated(arr%elm1))then
        call memory_allocate_array_dense(arr)
      else
        call lsquit("ERROR(array_cp_tiled2dense):dense is already allocated,&
        & please make sure you are not doing someting stupid",DECinfo%output)
      endif
      if(arr%init_type>0)pdm=.true.
      if(.not.present(order))call cp_tileddata2fort(arr,arr%elm1,arr%nelms,pdm)
      if(present(order))call cp_tileddata2fort(arr,arr%elm1,arr%nelms,pdm,order)
      if(change)arr%atype=DENSE
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
#ifdef VAR_LSMPI
    pdm=.false.
    if(.not.associated(arr%elm1))then
      call lsquit("ERROR(array_cp_dense2tiled):dense is NOT allocated,&
      & please make sure you are not doing someting stupid",DECinfo%output)
    endif
    if(arr%init_type>0)pdm=.true.
    if(change)arr%atype=TILED_DIST
    call array_convert_fort2arr(arr%elm1,arr,arr%nelms)
    call memory_deallocate_array_dense(arr)
#else
    return
#endif
  end subroutine array_mv_dense2tiled

!  subroutine array_convert_atype(arr,atype,usertdim,pdm)
!    implicit none
!    type(array), intent(inout) :: arr
!    integer, intent(in) :: atype
!    integer, intent(in), optional::usertdim(arr%mode)
!    integer, intent(in),optional :: pdm
!    integer :: tdim(arr%mode)
!    real(realk),pointer :: tmp(:)
!    tdim=DEFAULT_TDIM
!    if(present(usertdim))tdim=usertdim
!    if(atype==TILED_DIST.and.arr%atype==DENSE.and..not.present(pdm))then
!      call lsquit("ERROR(array_convert_atype): in a conversion from DENSE to&
!      & TILED_DIST the arguments usertdim and pdm need to be given",DECinfo%output)
!    endif
!
!    select case(arr%atype)
!    case(DENSE)
!      select case(atype)
!      case(DENSE)
!        call lsquit("ERROR(array_convert_atype):cowardly refusing to convert&
!        & from DENSE to DENSE",DECinfo%output)
!      case(TILED)
!        !arr=array_init_tiled(arr%dims,arr%mode,NO_PDM,.false.,arr%tdim)
!        arr%atype=TILED
!        call array_convert_fort2arr(arr%elm1,arr,arr%nelms)
!        call memory_deallocate_array_dense(arr)
!        stop 0
!      case(TILED_DIST)
!        call array_convert2pdm(arr,tdim,pdm)
!      end select
!    case(TILED)
!      select case(arr%atype)
!      case(DENSE)
!        call memory_allocate_array_dense(arr)
!        call array_convert_arr2fort(arr,arr%elm1,arr%nelms)
!        call memory_deallocate_tile(arr)
!          stop 0
!      case(TILED)
!        call lsquit("ERROR(array_convert_atype):cowardly refusing to convert&
!        & from TILED to TILED",DECinfo%output)
!      case(TILED_DIST)
!        print *,"should be simple, just try"
!        if(.not.present(pdm))call lsquit("ERROR(array_convert_atype):&
!        & pdm access type should be given",-1)
!              stop 0
!      end select
!    case(TILED_DIST)
!      select case(arr%atype)
!      case(DENSE)
!              print *,"easy to implement, but did not yet have time"
!              stop 0
!      case(TILED)
!              print *,"not yet implemented, somehow the order of the tiles on&
!              & the local node has to be handled"
!              stop 0
!      case(TILED_DIST)
!        print *,"if you want to change the distribution think of a clever&
!        & routine and implement it here, else contract and redistribute"
!        call lsquit("ERROR(array_convert_atype):cowardly refusing to convert&
!        & from TILED_DIST to TILED_DIST",DECinfo%output)
!      end select
!    end select
!  end subroutine array_convert_atype


  !\brief all the following wrappers are necessary to use the conversion routine
  !in an interface for different shapes of the fortran array  
  subroutine array_convert_arr2fort_wrapper1(arr,fort,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(inout) :: fort(:)
    integer, intent(in), optional :: order(arr%mode)
    if(present(order))call array_convert_arr2fort(arr,fort,arr%nelms,order)
    if(.not.present(order))call array_convert_arr2fort(arr,fort,arr%nelms)
  end subroutine array_convert_arr2fort_wrapper1
  subroutine array_convert_arr2fort_wrapper2(arr,fort,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(inout) :: fort(:,:)
    integer, intent(in), optional :: order(arr%mode)
    if(present(order))call array_convert_arr2fort(arr,fort,arr%nelms,order)
    if(.not.present(order))call array_convert_arr2fort(arr,fort,arr%nelms)
  end subroutine array_convert_arr2fort_wrapper2
  subroutine array_convert_arr2fort_wrapper3(arr,fort,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(inout) :: fort(:,:,:)
    integer, intent(in), optional :: order(arr%mode)
    if(present(order))call array_convert_arr2fort(arr,fort,arr%nelms,order)
    if(.not.present(order))call array_convert_arr2fort(arr,fort,arr%nelms)
  end subroutine array_convert_arr2fort_wrapper3
  subroutine array_convert_arr2fort_wrapper4(arr,fort,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(inout) :: fort(:,:,:,:)
    integer, intent(in), optional :: order(arr%mode)
    if(present(order))call array_convert_arr2fort(arr,fort,arr%nelms,order)
    if(.not.present(order))call array_convert_arr2fort(arr,fort,arr%nelms)
  end subroutine array_convert_arr2fort_wrapper4

  !\brief all the following wrappers are necessary to use the conversion routine
  !in an interface for different shapes of the fortran array  
  subroutine array_convert_fort2arr_wrapper1(fortarr,arr,order)
    implicit none
    type(array), intent(inout) :: arr
    real(realk), intent(in) :: fortarr(arr%nelms)
    integer, intent(in),optional :: order(arr%mode)
    if(present(order))call array_convert_fort2arr(fortarr,arr,arr%nelms,order)
    if(.not.present(order))call array_convert_fort2arr(fortarr,arr,arr%nelms)
  end subroutine array_convert_fort2arr_wrapper1
  subroutine array_convert_fort2arr_wrapper2(fortarr,arr,order)
    implicit none
    real(realk), intent(in) :: fortarr(:,:)
    type(array), intent(inout) :: arr
    integer, intent(in),optional :: order(arr%mode)
    if(present(order))call array_convert_fort2arr(fortarr,arr,arr%nelms,order)
    if(.not.present(order))call array_convert_fort2arr(fortarr,arr,arr%nelms)
  end subroutine array_convert_fort2arr_wrapper2
  subroutine array_convert_fort2arr_wrapper3(fortarr,arr,order)
    implicit none
    real(realk), intent(in) :: fortarr(:,:,:)
    type(array), intent(inout) :: arr
    integer, intent(in),optional :: order(arr%mode)
    if(present(order))call array_convert_fort2arr(fortarr,arr,arr%nelms,order)
    if(.not.present(order))call array_convert_fort2arr(fortarr,arr,arr%nelms)
  end subroutine array_convert_fort2arr_wrapper3
  subroutine array_convert_fort2arr_wrapper4(fortarr,arr,order)
    implicit none
    real(realk), intent(in) :: fortarr(:,:,:,:)
    type(array), intent(inout) :: arr
    integer, intent(in),optional :: order(arr%mode)
    if(present(order))call array_convert_fort2arr(fortarr,arr,arr%nelms,order)
    if(.not.present(order))call array_convert_fort2arr(fortarr,arr,arr%nelms)
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
    integer :: i,o(arr%mode)
    real(realk) :: nrm
    do  i=1,arr%mode
      o(i)=i
    enddo
    if(nelms/=arr%nelms)call lsquit("ERROR(array_convert_fort2arr):array&
    &dimensions are not the same",-1)
    select case(arr%atype)
      case(DENSE)
        call dcopy(int(nelms),fortarr,1,arr%elm1,1)
      case(TILED)
        call cp_data2tiled_lowmem(arr,fortarr,arr%dims,arr%mode)
      case(TILED_DIST)
        !if enough memory is available the lower one should be faster      
        if(arr%init_type==ALL_INIT)then
          do i=1,arr%nlti
            if(present(order))call extract_tile_from_fort(fortarr,arr%mode,arr%ti(i)%gt,arr%dims,&
                                              &arr%tdim,arr%ti(i)%t,order)
            if(.not.present(order))call extract_tile_from_fort(fortarr,arr%mode,arr%ti(i)%gt,arr%dims,&
                                              &arr%tdim,arr%ti(i)%t,o)
          enddo
        else
          !call cp_data2tiled_lowmem(arr,fortarr,arr%dims,arr%mode)
          if(present(order))call cp_data2tiled_intiles(arr,fortarr,arr%dims,arr%mode,order)
          if(.not.present(order))call cp_data2tiled_intiles(arr,fortarr,arr%dims,arr%mode)
        endif
      case default
        call lsquit("ERROR(array_convert_fort2arr) the array type is not implemented",-1)
    end select
  end subroutine array_convert_fort2arr
  !> \brief change the init type for a fortan array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine change_init_type(arr,totype)
    implicit none
    !> array to chage the init type
    type(array),intent(inout) :: arr
    !> type to change it to
    integer,intent(in) :: totype
    if(arr%atype==TILED_DIST.or.arr%atype==REPLICATED.or.&
         &totype==TILED_DIST.or.totype==REPLICATED)then
      call change_init_type_td(arr,totype)
    else
      print *,"what you want to do is not implemented"
      stop 0
    endif
  end subroutine change_init_type

  
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
    select case(arr%atype)
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
    
    select case(from_arr%atype)
      case(DENSE)
        select case(to_arr%atype)
          case(DENSE)
            call dcopy(int(from_arr%nelms),from_arr%elm1,1,to_arr%elm1,1)
          case(REPLICATED)
            call dcopy(int(from_arr%nelms),from_arr%elm1,1,to_arr%elm1,1)
            call array_sync_replicated(to_arr)
          case(TILED_DIST)
            !tilemem=to_arr%tsize*8.0E0_realk/(1024.0E0_realk**3)
            !call get_currently_available_memory(MemFree)
            !call cp_data2tiled_lowmem(to_arr,from_arr%elm1,from_arr%dims,from_arr%mode)
            call cp_data2tiled_intiles(to_arr,from_arr%elm1,from_arr%dims,from_arr%mode)
          case default
            call lsquit("ERROR(array_cp_data):operation not yet&
            & implemented",DECinfo%output)
        end select
      case(REPLICATED)
        select case(to_arr%atype)
          case(DENSE)
            call dcopy(int(from_arr%nelms),from_arr%elm1,1,to_arr%elm1,1)
          case(REPLICATED)
            call dcopy(int(from_arr%nelms),from_arr%elm1,1,to_arr%elm1,1)
            call array_sync_replicated(to_arr)
          case(TILED_DIST)
            !tilemem=to_arr%tsize*8.0E0_realk/(1024.0E0_realk**3)
            !call get_currently_available_memory(MemFree)
            !call cp_data2tiled_lowmem(to_arr,from_arr%elm1,from_arr%dims,from_arr%mode)
            call cp_data2tiled_intiles(to_arr,from_arr%elm1,from_arr%dims,from_arr%mode)
          case default
            call lsquit("ERROR(array_cp_data):operation not yet&
            & implemented",DECinfo%output)
        end select
      case(TILED_DIST)
        select case(to_arr%atype)
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
    
    select case(zeroed%atype)
      case(DENSE)
        zeroed%elm1=0.0E0_realk
      case(REPLICATED)
        zeroed%elm1=0.0E0_realk
        call array_sync_replicated(zeroed)
      case(TILED)
        do i=1,zeroed%ntiles
          zeroed%ti(i)%t=0.0E0_realk
        enddo
      case(TILED_DIST)
        call array_zero_tiled_dist(zeroed)
      case default
        print *,"ERROR(array_zero):not yet implemented"
        stop 0
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
    select case(arr%atype)
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
    if(.not.squareback)write(DECinfo%output,'("LOCAL TILE NORM ON",I3,f20.15)'),on,norm
    if(squareback)write(DECinfo%output,'("LOCAL TILE NORM^2 ON",I3,f20.15)'),on,norm
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
    select case(arr%atype)
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
            !  write(*,'("NaN detected in norm_t, tile:",I5," element: ",I5)'),i,j
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
    select case(arr%atype)
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
#ifdef VAR_LSMPI
    if(infpar%lg_mynum/=0)master=.false.
    nnod=infpar%lg_nodtot
#endif
    if(present(print_all_nodes))alln=print_all_nodes
    if(present(reducetocheck))then
      red=.true.
      if(red.and.master)then
#ifdef VAR_LSMPI
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY TESTCASES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_array_reorderings()
    implicit none

    type(array) :: test,test2
    real(realk),pointer :: datatata(:),tileget(:)
    real(realk),pointer :: tileget2(:),datata(:)
    real(realk) :: normher,ref1,ref2
    integer(kind=long) :: testint
    logical :: master,rigorous
    integer :: no,nv,nb,na,a,b,c,d
    character(len=7) :: teststatus
    master = .true.
    rigorous=.true.
    nb =  82
    nv =  81
    no =  80
    na =  0
    !for the upper parameters not all reorderings work!!!!!
    !nb =  9
    !nv =  8
    !no =  7
    !na =  6

    write(DECinfo%output,'(" Using",f8.3," GB of mem for the testarray")'),&
    &(nv*no*(nv+nb)*8.0E0_realk)/(1024.E0_realk*1024.E0_realk*1024.E0_realk)
    call mem_alloc(datatata,nb*nv*no)
    call mem_alloc(datata,nb*nv*no)
    testint=2
    call random_number(datatata)
    call print_norm(datatata,int(nb*nv*no,kind=8),ref1)
    ref2=0.5E0_realk*ref1
    write(DECinfo%output,'("REFERENCE NORM1:",f19.10)'),ref1
    write(DECinfo%output,'("REFERENCE NORM2:",f19.10)'),ref2
    write (DECinfo%output,*),""
    write (DECinfo%output,*),"TESTING 3D REORDERINGS"
    write (DECinfo%output,*),"**********************"
    write (DECinfo%output,*),""
    write (DECinfo%output,*),"reorder"

    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[1,2,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    if(rigorous)then
      do a=1,nb
        do b=1,nv
          do c=1,no
            if(abs(datatata(a+(b-1)*nb+(c-1)*nv*nb)-datata(a+(b-1)*nb+(c-1)*nv*nb))&
              & >1.0E-11_realk)teststatus=" FAILED"
          enddo
        enddo
      enddo
    endif
    write (DECinfo%output,'(" r 123: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[1,3,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    if(rigorous)then
      do a=1,nb
        do b=1,nv
          do c=1,no
            if(abs(datatata(a+(b-1)*nb+(c-1)*nv*nb)-datata(a+(c-1)*nb+(b-1)*no*nb))&
              & >1.0E-11_realk)teststatus=" FAILED"
          enddo
        enddo
      enddo
    endif
    write (DECinfo%output,'(" r 132: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[2,1,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    if(rigorous)then
      do a=1,nb
        do b=1,nv
          do c=1,no
            if(abs(datatata(a+(b-1)*nb+(c-1)*nv*nb)-datata(b+(a-1)*nv+(c-1)*nv*nb))&
              & >1.0E-11_realk)teststatus=" FAILED"
          enddo
        enddo
      enddo
    endif
    write (DECinfo%output,'(" r 213: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[2,3,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    if(rigorous)then
      do a=1,nb
        do b=1,nv
          do c=1,no
            if(abs(datatata(a+(b-1)*nb+(c-1)*nv*nb)-datata(b+(c-1)*nv+(a-1)*nv*no))&
              & >1.0E-11_realk)teststatus=" FAILED"
          enddo
        enddo
      enddo
    endif
    write (DECinfo%output,'(" r 231: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[3,1,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    if(rigorous)then
      do a=1,nb
        do b=1,nv
          do c=1,no
            if(abs(datatata(a+(b-1)*nb+(c-1)*nv*nb)-datata(c+(a-1)*no+(b-1)*nb*no))&
              & >1.0E-11_realk)teststatus=" FAILED"
          enddo
        enddo
      enddo
    endif
    write (DECinfo%output,'(" r 312: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[3,2,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    if(rigorous)then
      do a=1,nb
        do b=1,nv
          do c=1,no
            if(abs(datatata(a+(b-1)*nb+(c-1)*nv*nb)-datata(c+(b-1)*no+(a-1)*nv*no))&
              & >1.0E-11_realk)teststatus=" FAILED"
          enddo
        enddo
      enddo
    endif
    write (DECinfo%output,'(" r 321: ",f19.10," : ",A7)'),normher,teststatus


    write (DECinfo%output,*),""
    write (DECinfo%output,*),"reorder and add"
    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[1,2,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra123: ",f19.10," : ",A7)'),normher,teststatus
    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[1,3,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra132: ",f19.10," : ",A7)'),normher,teststatus

    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[2,1,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra213: ",f19.10," : ",A7)'),normher,teststatus
    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[2,3,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra231: ",f19.10," : ",A7)'),normher,teststatus

    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[3,1,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra312: ",f19.10," : ",A7)'),normher,teststatus
    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(1.0E0_realk,datatata,nb,nv,no,[3,2,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra321: ",f19.10," : ",A7)'),normher,teststatus


    write (DECinfo%output,*),""
    write (DECinfo%output,*),"scale and reorder"
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[1,2,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 123: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[1,3,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 132: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[2,1,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 213: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[2,3,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 231: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[3,1,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 312: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[3,2,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 321: ",f19.10," : ",A7)'),normher,teststatus


    write (DECinfo%output,*),""
    write (DECinfo%output,*),"scale, reorder and add"
    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[1,2,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra123: ",f19.10," : ",A7)'),normher,teststatus
    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[1,3,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra132: ",f19.10," : ",A7)'),normher,teststatus

    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[2,1,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra213: ",f19.10," : ",A7)'),normher,teststatus
    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[2,3,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra231: ",f19.10," : ",A7)'),normher,teststatus

    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[3,1,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra312: ",f19.10," : ",A7)'),normher,teststatus
    datata=0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_3d(0.5E0_realk,datatata,nb,nv,no,[3,2,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-11_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra321: ",f19.10," : ",A7)'),normher,teststatus

    call mem_dealloc(datatata)
    call mem_dealloc(datata)
    nb =  49
    nv =  37
    no =  31
    na =  29

    !nb =  30
    !nv =  30
    !no =  30
    !na =  30
    call mem_alloc(datatata,nb*na*nv*no)
    call mem_alloc(datata,nb*na*nv*no)
    call random_number(datatata)
    call print_norm(datatata,int(nb*na*nv*no,kind=8),ref1)
    ref2=0.5E0_realk*ref1
    write (DECinfo%output,*),""
    write(DECinfo%output,'("REFERENCE NORM3:",f19.10)'),ref1
    write(DECinfo%output,'("REFERENCE NORM4:",f19.10)'),ref2
    write (DECinfo%output,*),""
    write (DECinfo%output,*),"TESTING 4D REORDERINGS"
    write (DECinfo%output,*),"**********************"
    write (DECinfo%output,*),""
    write (DECinfo%output,*),"reorder"
    
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,2,3,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 1234: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,2,4,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 1243: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,3,2,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 1324: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,3,4,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 1342: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,4,2,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 1423: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,4,3,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 1432: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,1,3,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 2134: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,1,4,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 2143: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,3,1,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 2314: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,3,4,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 2341: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,4,1,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 2413: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,4,3,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 2431: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,1,2,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 3124: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,1,4,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 3142: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,2,1,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 3214: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,2,4,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 3241: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,4,1,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 3412: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,4,2,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 3421: ",f19.10," : ",A7)'),normher,teststatus
    
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,1,2,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 4123: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,1,3,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 4132: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,2,1,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 4213: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,2,3,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 4231: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,3,1,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 4312: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,3,2,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" r 4321: ",f19.10," : ",A7)'),normher,teststatus

    write (DECinfo%output,*),""
    write (DECinfo%output,*),"reorder and add"
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,2,3,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra1234: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,2,4,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra1243: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,3,2,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra1324: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,3,4,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra1342: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,4,2,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra1423: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[1,4,3,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra1432: ",f19.10," : ",A7)'),normher,teststatus

    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,1,3,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra2134: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,1,4,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra2143: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,3,1,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra2314: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,3,4,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra2341: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,4,1,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra2413: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[2,4,3,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra2431: ",f19.10," : ",A7)'),normher,teststatus

    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,1,2,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra3124: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,1,4,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra3142: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,2,1,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra3214: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,2,4,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra3241: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,4,1,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra3412: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[3,4,2,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra3421: ",f19.10," : ",A7)'),normher,teststatus
    
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,1,2,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra4123: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,1,3,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra4132: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,2,1,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra4213: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,2,3,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra4231: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,3,1,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra4312: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(1.0E0_realk,datatata,nb,na,nv,no,[4,3,2,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref1)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'(" ra4321: ",f19.10," : ",A7)'),normher,teststatus


    write (DECinfo%output,*),""
    write (DECinfo%output,*),"scale and reorder"
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,2,3,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 1234: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,2,4,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 1243: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,3,2,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 1324: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,3,4,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 1342: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,4,2,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 1423: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,4,3,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 1432: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,1,3,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 2134: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,1,4,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 2143: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,3,1,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 2314: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,3,4,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 2341: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,4,1,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 2413: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,4,3,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 2431: ",f19.10," : ",A7)'),normher,teststatus

    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,1,2,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 3124: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,1,4,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 3142: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,2,1,4],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 3214: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,2,4,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 3241: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,4,1,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 3412: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,4,2,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 3421: ",f19.10," : ",A7)'),normher,teststatus
    
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,1,2,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 4123: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,1,3,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 4132: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,2,1,3],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 4213: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,2,3,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 4231: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,3,1,2],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 4312: ",f19.10," : ",A7)'),normher,teststatus
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,3,2,1],0.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sr 4321: ",f19.10," : ",A7)'),normher,teststatus


    write (DECinfo%output,*),""
    write (DECinfo%output,*),"scale, reorder and add"
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,2,3,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra1234: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,2,4,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra1243: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,3,2,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra1324: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,3,4,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra1342: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,4,2,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra1423: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[1,4,3,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra1432: ",f19.10," : ",A7)'),normher,teststatus

    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,1,3,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra2134: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,1,4,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra2143: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,3,1,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra2314: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,3,4,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra2341: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,4,1,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra2413: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[2,4,3,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra2431: ",f19.10," : ",A7)'),normher,teststatus

    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,1,2,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra3124: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,1,4,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra3142: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,2,1,4],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra3214: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,2,4,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra3241: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,4,1,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra3412: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[3,4,2,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra3421: ",f19.10," : ",A7)'),normher,teststatus
    
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,1,2,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra4123: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,1,3,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra4132: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,2,1,3],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra4213: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,2,3,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra4231: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,3,1,2],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra4312: ",f19.10," : ",A7)'),normher,teststatus
    datata = 0.0E0_realk
    teststatus="SUCCESS"
    call array_reorder_4d(0.5E0_realk,datatata,nb,na,nv,no,[4,3,2,1],1.0E0_realk,datata)
    call print_norm(datata,int(nb*na*nv*no,kind=8),normher)
    if(abs(normher-ref2)>1.0E-10_realk)teststatus=" FAILED"
    write (DECinfo%output,'("sra4321: ",f19.10," : ",A7)'),normher,teststatus

    call mem_dealloc(datata)
    call mem_dealloc(datatata)
  end subroutine test_array_reorderings 

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
    integer :: no,nv,nb,na,i,j,succ
    integer(kind=ls_mpik) :: sender, recver
    character(len=7) :: teststatus
    character(ARR_MSG_LEN) :: msg
    master = .true.
#ifdef VAR_LSMPI
    if(infpar%lg_mynum /= 0) then
      master =.false.
    endif
#endif
    nb =  21
    nv =  18
    no =  12
    na =  7

#ifdef VAR_LSMPI
    if(master)then
      write(DECinfo%output,*),"TESTING PDM TILED ARRAY ALLOCATIONS"
      write(DECinfo%output,'(" Using",f8.3," GB of mem for the testarray")'),&
      &(nv*no*(nv+nb)*8.0E0_realk)/(1024.E0_realk*1024.E0_realk*1024.E0_realk)
      testint=2

      call mem_alloc(dummy1,nb*na*nv*no)
      call mem_alloc(dummy2,nb*na*nv*no)
      call random_number(dummy1)
      call print_norm(dummy1,int(nb*na*nv*no,kind=8),ref)
      write(DECinfo%output,'("REFERENCE NORM:",f19.12)'),ref
     

      !DIFFERENT ALLOCATION AND DEALLOCATION STEPS
      write (DECinfo%output,*), ""
      write (DECinfo%output,*), ""
      write (DECinfo%output,*), "TESTING SIMPLE ARRAY FUNCTIONS - MASTER DIRECTED"
      write (DECinfo%output,*), ""
      write (DECinfo%output,*),"ALLOC-DEALLOC TESTS"
      teststatus="SUCCESS"
      test=array_init([nv,na,nv,nb],4,TILED_DIST,MASTER_INIT,[nv,no-1,1,2])
      test2=array_init([na,nb,nv,no],4,TILED_DIST,MASTER_INIT,[nv,no-1,1,2])
      call array_free(test2)
      test2=array_init([no,no+1,no-1,no+1],4,TILED_DIST,MASTER_INIT,[no,no-1,nv,nb])
      call array_free(test)
      call array_free(test2)
      call array_print_mem_info(DECinfo%output,.true.,.false.,succ)
      if(succ/=0)teststatus=" FAILED"
      test2=array_init([nb,no,nv,no+1],4,TILED_DIST,MASTER_INIT,[nb,2,3,4])
      write (DECinfo%output,'(" ALLOC-DEALLOC TESTS: ",A7)'),teststatus  

      !ALLOCATING A FULL MATRIX AND PUT IT TO DISTRIBUTED MEMORY
      !check for errors via norm
      write(DECinfo%output,*),""
      write(DECinfo%output,*),""
      teststatus="SUCCESS"
      test=array_init([nb,na,nv,no],4,TILED_DIST,MASTER_INIT,[nb,na-1,3,no/2])
      write (DECinfo%output,*), "CONVERT PREVIOUS ARRAY TO PDM TILED" 
      call array_convert(dummy1,test,[1,2,3,4])
      call print_norm(test,normher)
      write(DECinfo%output,'("NORM OF PDM ARRAY  : ",f20.15)'),normher
      if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      write (DECinfo%output,'("CNVRT: NORM, TEST STATUS:",f19.10," : ",A7)'),normher,teststatus
      !GET A TILE OF A PDM ARRAY
      !calculate how many elements are in the desired tile, and allocate the
      !respective amount of memory in a fortran array
      write(DECinfo%output,*),""
      write(DECinfo%output,*),""
      write(DECinfo%output,*),"TESTING MPI_GET"
      testint=2
      call get_tile_dim(j,test,testint)
      call mem_alloc(tileget,j)
      call random_number(tileget)
      !initiatilize with some weird number here 10 and after get compare the
      !norms of the tile and the local fortran array
      teststatus="SUCCESS"
      call array_print_tile_norm(test,2,ref)
      write(DECinfo%output,'("NORM OF TILE IN ARRAY   : ",f20.15)'),ref
      call print_norm(tileget,int(j,kind=8),normher)
      write(DECinfo%output,'("NORM OF FORT BEFORE GET : ",f20.15)'),normher
      call array_get_tile(test,2,tileget,j)
      call print_norm(tileget,int(j,kind=8),normher)
      write(DECinfo%output,'("NORM OF FORT AFTER GET  : ",f20.15)'),normher
      if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      write (DECinfo%output,'("GET: NORM, TEST STATUS:  ",f20.15," : ",A7)'),normher,teststatus
      call mem_dealloc(tileget)

      write(DECinfo%output,*),""
      write(DECinfo%output,*),""
      write(DECinfo%output,*),"TESTING MPI_PUT"
      teststatus="SUCCESS"
      testint=2
      call get_tile_dim(j,test,testint)
      call mem_alloc(tileget,j)
      call random_number(tileget)
      !initiatilize with some weird number here 10 and after get compare the
      !norms of the tile and the local fortran array
      call array_print_tile_norm(test,2,normher)
      write(DECinfo%output,'("NORM OF TILE BEFORE PUT : ",f20.15)'),normher
      call print_norm(tileget,int(j,kind=8),ref)
      write(DECinfo%output,'("NORM OF FORT TO PUT     : ",f20.15)'),ref
      call array_put_tile(test,2,tileget,j)
      call array_print_tile_norm(test,2,normher)
      write(DECinfo%output,'("NORM OF TILE AFTER PUT  : ",f20.15)'),normher
      if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      write (DECinfo%output,'("PUT: NORM, TEST STATUS:  ",f20.15," : ",A7)'),normher,teststatus
     

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
      write(DECinfo%output,*),""
      write(DECinfo%output,*),""
      write(DECinfo%output,*),"TESTING MPI_ACCUMULATE"
      teststatus="SUCCESS"
      do i=1,j
        tileget(i)=tileget(i)+3.0E0_realk
      enddo
      call print_norm(tileget,int(j,kind=8),ref)
      write(DECinfo%output,'("NORM LOCAL ACCUMULATION : ",f20.15)'),ref
      !initialize the local tile with 3 and accumulate it --> compare norm
      tileget=3.0E0_realk
      call print_norm(tileget,int(j,kind=8),normher)
      write(DECinfo%output,'("NORM OF FORT TO ADD:      ",f20.15)'),normher
      call array_accumulate_tile(test,2,tileget,j)
      call array_print_tile_norm(test,2,normher)
      write(DECinfo%output,'("NORM REMOTE ACCUMULATION: ",f20.15)'),normher
      !use the tile with three in it, print its norm put and compare norms
      if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      write (DECinfo%output,'("ACC: NORM, TEST STATUS:  ",f20.15," : ",A7)'),normher,teststatus


      call array_free(test)
      test=array_init([nb,na,nv,no],4,TILED_DIST,MASTER_INIT,[0,0,0,0])
      write(DECinfo%output,*),""
      write(DECinfo%output,*),""
      write(DECinfo%output,*),"TESTING CONVERSION TO FORT"
      call array_convert(dummy1,test)
      teststatus="SUCCESS"
      call print_norm(dummy1,int(nb*na*nv*no,kind=8),ref)
      write(DECinfo%output,'("NORM OF DENSE ARRAY:      ",f20.15)'),ref
      call print_norm(dummy1,int(no*nv*na*nb,kind=8),normher)
      write(DECinfo%output,'("NORM OF PDM ARRAY :       ",f20.15)'),normher
      dummy2=1.0E13_realk
      call array_convert(test,dummy2)
      call print_norm(dummy2,int(no*nv*na*nb,kind=8),normher)
      write(DECinfo%output,'("NORM OF CONTRACTED ARRAY: ",f20.15)'),normher
      if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
      write (DECinfo%output,'("CTR: NORM, TEST STATUS:  ",f20.15," : ",A7)'),normher,teststatus
      teststatus="SUCCESS"
      do i=1,no*nv*na*nb
        if(abs(dummy1(i)-dummy2(i))>1.0E-12)then
          print *,"element",i,dummy1(i),dummy2(i)
          teststatus=" FAILED"
        endif
      enddo
      write (DECinfo%output,'("ORDER: TEST STATUS:                              ",A7)'),teststatus



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
       write (DECinfo%output,'("FIRST HALF ALLOCATION: ",A7)'),teststatus  
    endif
    !get the slaves into this routine
    if(master)then
      print *,"MASTER GETTING SLAVES"
      call ls_mpibcast(ARRAYTEST,infpar%master,infpar%lg_comm)
      write (DECinfo%output,*),""
      write (DECinfo%output,*),""
      write (DECinfo%output,*),"TESTING PARALLEL ACCESS TO THE SAME ROUTINES"
      write (DECinfo%output,*),""
    else
      print *,"SLAVE ARRIVED",infpar%lg_mynum
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!ALL OF THE SLAVES WILL BE HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    testint=2

    !initialize a matrix
    teststatus="SUCCESS"
    if(master) write (DECinfo%output,*),"ALLOC-DEALLOC TESTS"
    test=array_init([nb,nb+2,nb+3,nb+4],4,TILED_DIST,ALL_INIT,[nb,nb+2,40,2])
    test2=array_init([no+3,no+2,no+1,no],4,TILED_DIST,ALL_INIT,[no,40,40,10])
    call array_free(test)
    call array_free(test2)
    test2=array_init([nb,na,nv,no],4,TILED_DIST,ALL_INIT,[nb,nv-1,1,2])
    call array_free(test2)
    call array_print_mem_info(DECinfo%output,.true.,.true.,succ)
    if(succ/=0)teststatus=" FAILED"
    test2=array_init([nb,no,nv,no+1],4,TILED_DIST,ALL_INIT,[nb,2,3,4])
    if(master) write (DECinfo%output,'(" ALLOC-DEALLOC TESTS: ",A7)'),teststatus  
    if(master) write (DECinfo%output,*),"DONE -- NOW COMMUNICATION"
    if(master) write(DECinfo%output,*),""
    if(master) write(DECinfo%output,*),""
    !call lsmpi_barrier(infpar%lg_comm)

    !IF MY RANK IS THREE, PUT A MATRIX CONTAINING 10 in TILE 2 (ON THE
    !RESPECTIVE THREAD) 
    teststatus="SUCCESS"
    if(infpar%lg_mynum==3.or.master)then
      recver=3
      if(.not.master)then
        call get_tile_dim(j,test2,testint)
        call mem_alloc(tileget,j)
        tileget = 1.0E1_realk
        call array_put_tile(test2,2,tileget,j)
        call print_norm(tileget,int(j,kind=8),normher)
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        call mem_dealloc(tileget)
      else
        call ls_mpisendrecv(ref,infpar%lg_comm,recver,infpar%master)
        write(DECinfo%output,'("NORM PARALLEL 3LPN: ",f20.15)'),ref
      endif
    endif
    !BEFORE 2 CAN GET THE 
    call lsmpi_barrier(infpar%lg_comm)
    if(infpar%lg_mynum==2.or.master)then
      recver=2
      if(.not.master)then
        call get_tile_dim(j,test2,testint)
        call mem_alloc(tileget,j)
        call array_get_tile(test2,2,tileget,j)
        call print_norm(tileget,int(j,kind=8),normher)
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        do i=1,j
          tileget(i) = tileget(i) + 2.4E0_realk
        enddo
        call print_norm(tileget,int(j,kind=8),normher)
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        tileget = 2.4E0_realk
        call array_accumulate_tile(test2,[2,1,1,1],tileget,j)
        call mem_dealloc(tileget)
      else
        teststatus="SUCCESS"
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        write(DECinfo%output,'("NORM PARALLEL 2LGN: ",f20.15)'),normher
        if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
        write (DECinfo%output,'("PUT-GET: NORM, TEST STATUS: ",f19.10," : ",A7)'),normher,teststatus

        call ls_mpisendrecv(ref,infpar%lg_comm,recver,infpar%master)
        write(DECinfo%output,'("NORM PARALLEL 2LAC: ",f20.15)'),ref
      endif
    endif
    !BE CAREFUL ABOUT WHETER THE INFORMATION IS ALREADY TRANSMITTED --> AT
    !CRITICAL POINTS INSERT BARRIER STATEMENTS TO SYNCHONIZE THE NODES 
    call lsmpi_barrier(infpar%lg_comm)
    call array_print_tile_norm(test2,2,normher)
    call array_free(test2)
    if(master)then
       teststatus="SUCCESS"
       write(DECinfo%output,'("NORM PARALLEL WORK: ",f20.15)'),normher
       if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
       write (DECinfo%output,'("ACC2    : NORM, TEST STATUS: ",f19.10," : ",A7)'),ref,teststatus
    endif
    if(master) write (DECinfo%output,*),""

  
    !test extracting a tile with a different ordering than the dense matrix, put
    !that into pdm, get these tiles on each node and put them in reversed
    !reordering back into the full array, check norms and order
    teststatus="SUCCESS"
    call lsmpi_barrier(infpar%lg_comm)
    test2=array_init([no-4,nv+3,nv/7,no],4,TILED_DIST,ALL_INIT,[no-4,nv+3,5,2])
    test=array_init([nv/7,nv+3,no,no-4],4,TILED_DIST,ALL_INIT)
    call memory_allocate_array_dense(test)
    call random_number(test%elm1)
    call lsmpi_local_allreduce(test%elm1,test%nelms)
    if(infpar%lg_mynum==0)then
      write (msg,*),"local test norm master"
      call print_norm(test%elm1,test%nelms,msg)
    endif
    if(infpar%lg_mynum==1)then
      write (msg,*),"local test norm slave"
      call print_norm(test%elm1,test%nelms,msg)
    endif
    call print_norm(test%elm1,test%nelms,ref)
    call array_convert(test%elm1,test2,[4,2,1,3])
    call lsmpi_barrier(infpar%lg_comm)
    call array_mv_dense2tiled(test,.false.)
    call memory_allocate_array_dense(test2)
    call lsmpi_barrier(infpar%lg_comm)
    test2%elm1=0.0E0_realk
    do i=1,test2%ntiles
      call get_tile_dim(j,test2,i)
      call mem_alloc(tileget,j)
      call array_get_tile(test2,i,tileget,j)
      call put_tile_in_fort(tileget,i,test2%tdim,test2%elm1,test2%dims,4,[4,2,1,3])
      call mem_dealloc(tileget)
    enddo
    call lsmpi_barrier(infpar%lg_comm)
    call array_cp_tiled2dense(test,.true.)
    call lsmpi_barrier(infpar%lg_comm)
    if(infpar%lg_mynum==0)then
      write (msg,*),"local test 2 norm master"
      call print_norm(test2%elm1,test2%nelms,msg)
    endif
    if(infpar%lg_mynum==1)then
      write (msg,*),"local test 2 norm slave"
      call print_norm(test2%elm1,test2%nelms,msg)
    endif
    call print_norm(test2%elm1,test2%nelms,normher)
    do i=1,test%nelms
      if(abs(test%elm1(i)-test2%elm1(i))>1.0E-12)then
        teststatus=" FAILED"
      endif
    enddo
    call memory_deallocate_array_dense(test)
    call memory_deallocate_array_dense(test2)
    test%atype=TILED_DIST
    test2%atype=TILED_DIST
    call array_free(test)
    call array_free(test2)
    if(master)then
       write(DECinfo%output,'("PDM REORDERINGS: ",f20.15)'),normher
       if(abs(normher-ref)>1.0E-12_realk)teststatus=" FAILED"
       write (DECinfo%output,'("PDMR    : NORM, TEST STATUS: ",f19.10," : ",A7)'),ref,teststatus
    endif
    if(master) write (DECinfo%output,*),""
#endif

  end subroutine test_array_struct 

end module tensor_interface_module


#ifdef VAR_LSMPI
subroutine get_slaves_to_array_test()
  use precision
  use tensor_interface_module,only:test_array_struct
  implicit none
  call test_array_struct()
end subroutine get_slaves_to_array_test
#endif