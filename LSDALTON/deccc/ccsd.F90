!> @file
!> Subroutines related with construction of CC doubles residual
!> \author Patrick Ettenhuber, Marcin Ziolkowski, and Kasper Kristensen
module ccsd_module


  use precision
  use ptr_assoc_module!,only:ass_D4to1,ass_D2to1,ass_D1to3
  use lstiming!, only: lstimer
  use typedeftype!, only: lsitem,lssetting
  use matrix_module!, only:matrix
  use matrix_operations!, only: mat_init, mat_zero, mat_free
  use screen_mod!, only: DECscreenITEM
  use memory_handling!, only: mem_dealloc, mem_alloc
  use dec_typedef_module
  use BUILDAOBATCH!,only:build_batchesofaos,determine_maxbatchorbitalsize,&
!       & determine_MaxOrbitals
  use screen_mod!,only: free_decscreen, DECSCREENITEM
  use integralinterfaceDEC!,only: ii_precalc_decscreenmat, &
!       & ii_getbatchorbitalscreen, ii_get_decpacked4center_j_eri, &
!       & ii_getbatchorbitalscreenk, ii_get_decpacked4center_k_eri
  use integralinterfaceMod!, only: ii_get_h1, ii_get_h1_mixed_full,&
!       & ii_get_fock_mat_full
#ifdef VAR_LSMPI
  use infpar_module
  use lsmpi_type
#endif
  use integralparameters!, only: AORdefault
  use tensor_interface_module!, only: precondition_doubles_parallel
  use lspdm_tensor_operations_module!, only: array_init, array_change_atype_to_rep,&
  use tensor_basic_module!, only: DENSE,TILED,TILED_DIST,SCALAPACK,&
!         & NO_PDM,MASTER_INIT,REPLICATED,ALL_INIT,ass_1to3,ass_1to2,&
!         & ass_1to4,ass_2to1,&
!         & ass_4to1,ARR_MSG_LEN

    ! DEC DEPENDENCIES (within deccc directory)   
    ! *****************************************
#ifdef VAR_LSMPI
  use decmpi_module!, only: mpi_communicate_ccsd_calcdata,distribute_mpi_jobs
#endif
    use dec_fragment_utils
    use ri_simple_operations
    use array2_simple_operations!, only: array2_init, array2_add,&
!         & array2_transpose, array2_free, array2_add_to
    use array3_simple_operations!, only: array_reorder_3d
    use array4_simple_operations!, only: array4_init, operator(*),&
!         & array_reorder_4d, mat_transpose, array4_contract1,&
!         & array4_reorder, array4_free, array4_contract2_middle,&
!         & array4_read, array4_contract2, mat_transpose_p,mat_transpose_pl,&
!         & array4_alloc, array4_add_to, array4_scale, array4_contract3,&
!         & array4_read_file_type2, array4_write_file_type2,&
!         & array4_open_file, array4_read_file, array4_close_file,&
!         & array4_write_file
!         & array_change_atype_to_d,print_norm
    use ccintegrals!, only: get_gmo_simple,getL,dec_fock_transformation


    public :: getDoublesResidualMP2_simple, getDoublesResidualCCSD_simple, &
         & getDoublesResidual_explicite, getDoublesResidualCCSD_simple2, &
         & get_ccsd_residual_integral_driven_oldarray_wrapper, get_ccsd_residual_integral_driven, &
         & getFockCorrection, getInactiveFockFromRI,getInactiveFock_simple, &
         & precondition_singles, precondition_doubles,get_aot1fock, get_fock_matrix_for_dec, &
         & gett1transformation, getsinglesresidualccsd,fullmolecular_get_aot1fock,calculate_E2_and_permute, &
         & get_max_batch_sizes
    private

  interface Get_AOt1Fock
    module procedure Get_AOt1Fock_arraywrapper,Get_AOt1Fock_oa
  end interface Get_AOt1Fock
  interface get_fock_matrix_for_dec
    module procedure get_fock_matrix_for_dec_oa,get_fock_matrix_for_dec_arraywrapper
  end interface get_fock_matrix_for_dec

  interface precondition_singles
    module procedure precondition_singles_newarr,&
                    &precondition_singles_oldarr
  end interface precondition_singles

    interface precondition_doubles
      module procedure precondition_doubles_newarr,&
                      &precondition_doubles_oldarr
    end interface precondition_doubles
    

contains
  function precondition_doubles_newarr(omega2,ppfock,qqfock) result(prec)

    integer none
    type(array), intent(in) :: omega2
    type(array), intent(inout) :: ppfock, qqfock
    type(array) :: prec
    integer, dimension(4) :: dims
    integer :: a,i,b,j

    if(omega2%mode/=4.or.ppfock%mode/=2.or.qqfock%mode/=2)then
      call lsquit("ERROR(precondition_doubles_newarr):wrong number of modes&
      & for this operation",DECinfo%output)
    endif
    dims = omega2%dims

    if(omega2%atype==DENSE.and.&
    &(ppfock%atype==DENSE.or.ppfock%atype==REPLICATED).and.&
    &(qqfock%atype==DENSE.or.qqfock%atype==REPLICATED))then
      prec = array_init(dims,4)
     
      do j=1,dims(4)
        do i=1,dims(3)
          do b=1,dims(2)
            do a=1,dims(1)
     
             prec%elm4(a,b,i,j) = omega2%elm4(a,b,i,j) / &
                  ( ppfock%elm2(i,i) - qqfock%elm2(a,a) + &
                    ppfock%elm2(j,j) - qqfock%elm2(b,b) )
     
            end do
          end do
        end do
      end do

    elseif(omega2%atype==TILED_DIST.and.&
          &associated(ppfock%addr_p_arr).and.associated(qqfock%addr_p_arr))then
      prec = array_init(dims,4,TILED_DIST,MASTER_INIT,omega2%tdim)
      call array_change_atype_to_rep(ppfock)
      call array_change_atype_to_rep(qqfock)
      call precondition_doubles_parallel(omega2,ppfock,qqfock,prec)
      call array_change_atype_to_d(ppfock)
      call array_change_atype_to_d(qqfock)

    else
      call lsquit("ERROR(precondition_doubles_newarr):No preconditioning routine&
      & available for your choice of arrays",DECinfo%output)
    endif

  end function precondition_doubles_newarr


  !> \brief Preconditioning of doubles residual interface
  !> \author Kasper Kristensen
  !> \date October 2010
  function precondition_doubles_oldarr(omega2,ppfock,qqfock) result(prec)

    integer none
    type(array4), intent(in) :: omega2
    type(array2), intent(in) :: ppfock, qqfock
    type(array4) :: prec

    if(DECinfo%array4OnFile) then ! omega2 is written from file
       prec = precondition_doubles_file(omega2,ppfock,qqfock)
    else ! omega2 is stored in memory
       prec = precondition_doubles_memory(omega2,ppfock,qqfock)
    end if

  end function precondition_doubles_oldarr



  !> \brief Preconditioning of doucles residual
  !> when arrays are kept in memory.
  !> \author Marcin Ziolkowski
  function precondition_doubles_memory(omega2,ppfock,qqfock) result(prec)

    integer none
    type(array4), intent(in) :: omega2
    type(array2), intent(in) :: ppfock, qqfock
    type(array4) :: prec
    integer, dimension(4) :: dims
    integer :: a,i,b,j

    dims = omega2%dims
    prec = array4_init(dims)

    do a=1,dims(1)
       do i=1,dims(2)
          do b=1,dims(3)
             do j=1,dims(4)

                prec%val(a,i,b,j) = omega2%val(a,i,b,j) / &
                     ( ppfock%val(i,i) - qqfock%val(a,a) + &
                     ppfock%val(j,j) - qqfock%val(b,b) )

             end do
          end do
       end do
    end do

    return
  end function precondition_doubles_memory



  !> \brief Preconditioning of doubles residual
  !> when arrays are kept in memory.
  !> \author Kasper Kristensen
  function precondition_doubles_file(omega2,ppfock,qqfock) result(prec)

    integer none
    type(array4), intent(in) :: omega2
    type(array2), intent(in) :: ppfock, qqfock
    type(array4) :: prec
    integer :: a,i,b,j,dim1,dim2,dim3,dim4
    real(realk), pointer :: omega2val(:,:), precval(:,:)



    ! Sanity check
    ! ************

    ! Only implemented for storing type 2
    if(omega2%storing_type /= 2) then
       call lsquit('precondition_doubles_file: &
            & Only implemented for storing type 2!',-1)
    end if



    ! Initialize stuff
    ! ****************
    dim1 = omega2%dims(1)
    dim2 = omega2%dims(2)
    dim3 = omega2%dims(3)
    dim4 = omega2%dims(4)
    call mem_alloc(omega2val,dim1,dim2)
    call mem_alloc(precval,dim1,dim2)

    ! Use storing type 2 for preconditoned residual
    prec = array4_init(omega2%dims,2,.false.)

    ! Open file for residual and preconditioned residual
    call array4_open_file(omega2)
    call array4_open_file(prec)



    ! Carry out preconditioning
    ! *************************


    i_loop: do i=1,dim4
       a_loop: do a=1,dim3

          ! Read omega2(:,:,a,i) into omega2val
          call array4_read_file(omega2,a,i,omega2val,dim1,dim2)

          j_loop: do j=1,dim2
             b_loop: do b=1,dim1

                ! Calculate preconditioned residual: Element (b,j,a,i)
                precval(b,j) = omega2val(b,j) / &
                     & ( ppfock%val(i,i) - qqfock%val(a,a) + &
                     & ppfock%val(j,j) - qqfock%val(b,b) )

             end do b_loop
          end do j_loop

          ! At this point prec(:,:,a,i) has been calculated and stored in precval.
          ! Write these elements to the file referenced by prec using storing type 2.
          call array4_write_file(prec,a,i,precval,dim1,dim2)

       end do a_loop
    end do i_loop


    ! Free stuff
    call mem_dealloc(omega2val)
    call mem_dealloc(precval)
    call array4_close_file(omega2,'keep')
    call array4_close_file(prec,'keep')


  end function precondition_doubles_file





  !> \brief Doubles residual for MP2, transformed integrals enter this routine and don't change
  subroutine getDoublesResidualMP2_simple(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    type(array4), intent(in) :: gmo
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    type(array4) :: tmp
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k

    if(DECinfo%array4OnFile) then
       call getDoublesResidualMP2_OnFile(omega2,&
            &t2,gmo,pfock,qfock,nocc,nvirt)
    else
       call getDoublesResidualMP2_InMemory(omega2,&
            &t2,gmo,pfock,qfock,nocc,nvirt)
    end if

  end subroutine getDoublesResidualMP2_simple



  !> \brief Doubles residual for MP2, transformed integrals enter this routine and don't change
  !> Used when array4 elements are stored in memory.
  subroutine getDoublesResidualMP2_InMemory(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    type(array4), intent(in) :: gmo
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    type(array4) :: tmp
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k


    ! 1
    call array2_transpose(qfock)
    call array4_reorder(t2,[3,4,1,2])
    tmp = array4_init([nvirt,nocc,nvirt,nocc])
    call array4_contract1(t2,qfock,tmp,.true.)
    call array4_reorder(tmp,[3,4,1,2])
    call array4_add_to(omega2,1.0E0_realk,tmp)
    call array4_free(tmp)
    call array4_reorder(t2,[3,4,1,2])

    ! 2
    call array4_contract1(t2,qfock,omega2,.false.)
    call array2_transpose(qfock)

    ! 3
    call array4_reorder(t2,[4,3,2,1])
    tmp = array4_init([nocc,nvirt,nocc,nvirt])
    call array4_contract1(t2,pfock,tmp,.true.)
    call array4_reorder(t2,[4,3,2,1])
    call array4_reorder(tmp,[4,3,2,1])
    call array4_add_to(omega2,-1.0E0_realk,tmp)
    call array4_free(tmp)

    ! 4
    call array4_reorder(t2,[2,1,3,4])
    tmp = array4_init([nocc,nvirt,nvirt,nocc])
    call array4_contract1(t2,pfock,tmp,.true.)
    call array4_reorder(t2,[2,1,3,4])
    call array4_reorder(tmp,[2,1,3,4])
    call array4_add_to(omega2,-1.0E0_realk,tmp)
    call array4_free(tmp)

    ! 5
    call array4_add_to(omega2,1.0E0_realk,gmo)


    return
  end subroutine getDoublesResidualMP2_InMemory




  !> \brief Calculates doubles residual for MP2 when array4 elements are stored on file.
  !> Detailed equations are given inside subroutine.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine getDoublesResidualMP2_OnFile(omega2,t2,RHS,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(in) :: t2
    type(array4), intent(in) :: RHS
    type(array4), intent(inout) :: Omega2
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    type(array4) :: OmegaA, OmegaB
    integer :: a,b,i,j,dimsAB(4)
    real(realk), pointer :: t2tmp(:,:)
    real(realk),pointer :: OmegaA_ij(:,:), OmegaA_ji(:,:)
    real(realk),pointer :: OmegaB_ij(:,:), OmegaB_ji(:,:)
    real(realk),pointer :: OmegaAsym_ji(:,:,:), OmegaBsym_ji(:,:,:)
    real(realk), pointer :: RHStmp(:,:), Omega2tmp(:,:)
    real(realk), pointer :: OmegaAval(:,:,:), OmegaBval(:,:,:)
    real(realk) :: tcpu, twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)


    ! *******************************************************************
    !                    MP2 residual - main equations                  *
    ! *******************************************************************
    !
    ! Amplitudes and integrals satisfy:
    ! t_{aibj} = t_{bjai}
    ! RHS_{aibj} = RHS_{bjai}
    !
    !
    ! Residual Omega2:
    ! ''''''''''''''''
    ! Omega2_{bjai} = RHS_{bjai}
    !               + sum_{c} t_{bjci} F_{ca}       (OmegaA_{bjai})
    !               + sum_{c} t_{aicj} F_{cb}       (OmegaA_{aibj})
    !               - sum_{k} t_{bjak} F_{ki}       (OmegaB_{bjai})
    !               - sum_{k} t_{aibk} F_{kj}       (OmegaB_{aibj})
    !               = RHS_{bjai}
    !               + OmegaA_{bjai}                 ! OmegaAsym_{aibj} = OmegaA_{aibj}
    !               + OmegaA_{aibj}                 !                  + OmegaA_{bjai}
    !               - OmegaB_{bjai}                 ! Same definition for OmegaBsym
    !               - OmegaB_{aibj}
    !               = RHS_{bjai}
    !               + OmegaAsym_{bjai}
    !               - OmegaBsym_{bjai}              (*)
    !
    !  STRATEGY:
    !  1. Contruct OmegaA and OmegaB
    !  2. Symmetrize OmegaA and OmegaB and add to RHS.



    ! ****************************************************
    ! *       1. Construction of OmegaA and OmegaB       *
    ! ****************************************************


    ! Initialize stuff
    ! ****************

    ! Initialize memory for temporary arrays
    call mem_alloc(t2tmp,nvirt,nocc)
    call mem_alloc(OmegaAval,nvirt,nocc,nvirt)
    call mem_alloc(OmegaBval,nocc,nvirt,nvirt)

    ! Dimension for storing omegaA and omegaB
    dimsAB = [nvirt,nvirt,nocc,nocc]    ! order VVOO
    ! Initialize OmegaA and OmegaB using storing type 2
    OmegaA = array4_init(dimsAB,2,.false.)
    OmegaB = array4_init(dimsAB,2,.false.)

    ! Open amplitude file, and files for writing OmegaA and OmegaB
    call array4_open_file(OmegaA)
    call array4_open_file(OmegaB)
    call array4_open_file(t2)

    i_loop: do i=1,nocc

       a_loop: do a=1,nvirt

          ! Read t2(:,:,a,i) into t2tmp
          call array4_read_file_type2(t2,a,i,t2tmp,nvirt,nocc)


          ! Construct and store OmegaA
          ! **************************

          ! Matrix product, for each (ai) index:
          ! OmegaA_{bj} = sum_{c} Fbc * t_{cj}
          call dgemm('n','n',nvirt,nocc,nvirt,1E0_realk,qfock%val, &
               & nvirt,t2tmp,nvirt,0E0_realk,OmegaAval(:,:,a),nvirt)


          ! Construct and store OmegaB
          ! **************************

          ! Matrix product, for each ai index:
          ! OmegaB_{jb} = sum_{k} t_{bk} * F_{kj} = sum_{k} F_{jk} t_{kb}
          ! F is symmetric, t(virt,occ) is not. -> transpose t to (occ,virt).
          call dgemm('n','t',nocc,nvirt,nocc,1E0_realk,pfock%val, &
               & nocc, t2tmp,nvirt,0E0_realk,omegaBval(:,:,a),nocc)

       end do a_loop

       ! Now OmegaA(:,:,:,i) is stored in memory as OmegaAval(b,j,a) (for given "i")
       ! and OmegaB(:,:,:,i) is stored in memory as OmegaBval(j,b,a) (for given "i")
       ! We next store both of these to file in (j,i) chunks, i.e.
       ! Ordering: (b,a,j,i) for all a and b, and a particular (j,i) set.
       ! (Refering back to the amplitudes, b goes with j, and a goes with i).
       ! Same for both OmegaA and OmegaB!
       ! Very important that this particular storing is used to be able to symmetrize below.
       do j=1,nocc
          ! Write OmegaA(b,a) for given j and i, i.e. OmegaAval(:,j,:)
          call array4_write_file_type2(OmegaA,j,i,OmegaAval(1:nvirt,j,1:nvirt),&
               &dimsAB(1), dimsAB(2) )
          ! Write OmegaB(b,a) for given j and i, i.e. OmegaBval(j,:,:)
          call array4_write_file_type2(OmegaB,j,i,OmegaBval(j,1:nvirt,1:nvirt),&
               &dimsAB(1), dimsAB(2) )
       end do

    end do i_loop


    ! Free stuff and close t2 file
    call mem_dealloc(t2tmp)
    call mem_dealloc(OmegaAval)
    call mem_dealloc(OmegaBval)
    call array4_close_file(t2,'keep')


    ! At this point OmegaA and omegaB are stored on file using the ordering:
    ! (b,a,j,i) in (j,i) chunks using storing type 2.

    if (DECinfo%PL>1) call LSTIMER('RES: STEP 1',tcpu,twall,DECinfo%output)



    ! ***********************************************************************
    ! *      2. Symmetrize OmegaA and OmegaB, then calculate residual       *
    ! ***********************************************************************


    ! Initialize stuff
    ! ****************

    ! Allocate memory for temporary arrays
    call mem_alloc(OmegaA_ij,nvirt,nvirt)
    call mem_alloc(OmegaA_ji,nvirt,nvirt)
    call mem_alloc(OmegaB_ij,nvirt,nvirt)
    call mem_alloc(OmegaB_ji,nvirt,nvirt)
    call mem_alloc(OmegaAsym_ji,nvirt,nvirt,nocc)
    call mem_alloc(OmegaBsym_ji,nvirt,nvirt,nocc)
    call mem_alloc(RHStmp,nvirt,nocc)
    call mem_alloc(Omega2tmp,nvirt,nocc)

    ! Open files for final residual array omega2 and RHS array
    call array4_open_file(omega2)
    call array4_open_file(RHS)


    DoI : do i=1,nocc
       DoJ: do j=1,nocc

          ! Read in OmegaA(b,a,j,i) and OmegaB(b,a,j,i) for all (b,a) and a particular (j,i) set.
          call array4_read_file_type2(OmegaA,j,i,OmegaA_ji(1:nvirt,1:nvirt),&
               &OmegaA%dims(1),OmegaA%dims(2))
          call array4_read_file_type2(OmegaB,j,i,OmegaB_ji(1:nvirt,1:nvirt),&
               &OmegaB%dims(1),OmegaB%dims(2))
          ! Similarly for (i,j) (but only if non-diagonal)
          if(i /= j) then
             call array4_read_file_type2(OmegaA,i,j,OmegaA_ij(:,:),OmegaA%dims(1),OmegaA%dims(2))
             call array4_read_file_type2(OmegaB,i,j,OmegaB_ij(:,:),OmegaB%dims(1),OmegaB%dims(2))
          end if

          ! Now OmegaA(:,:,j,i) and OmegaA(:,:,i,j) are in memory
          ! and symmetrization can be performed, see Eq. (*).
          ! (Similarly for OmegaB.)
          symA_a: do a=1,nvirt
             symA_b: do b=1,nvirt

                if(i==j) then ! Diagonal
                   OmegaAsym_ji(b,a,j) = OmegaA_ji(b,a) + OmegaA_ji(a,b)
                   OmegaBsym_ji(b,a,j) = OmegaB_ji(b,a) + OmegaB_ji(a,b)
                else ! Off diagonal
                   OmegaAsym_ji(b,a,j) = OmegaA_ji(b,a) + OmegaA_ij(a,b)
                   OmegaBsym_ji(b,a,j) = OmegaB_ji(b,a) + OmegaB_ij(a,b)
                end if

             end do symA_b
          end do symA_a

       end do DoJ

       ! Now OmegaA(:,:,:,i) and OmegaB(:,:,:,i) are in memory,
       ! and we can construct the final residual, stored in (a,i) chunks.
       do a=1,nvirt
          call array4_read_file_type2(RHS,a,i,RHStmp,nvirt,nocc)

          ! Residual according to (*) above
          Omega2tmp(1:nvirt,1:nocc) = RHStmp(1:nvirt,1:nocc) &
               & +OmegaAsym_ji(1:nvirt,a,1:nocc) -OmegaBsym_ji(1:nvirt,a,1:nocc)

          call array4_write_file_type2(Omega2,a,i,Omega2tmp,nvirt,nocc)
       end do

    end do DoI


    ! Free stuff and close files
    call mem_dealloc(OmegaA_ij)
    call mem_dealloc(OmegaA_ji)
    call mem_dealloc(OmegaAsym_ji)
    call mem_dealloc(OmegaB_ij)
    call mem_dealloc(OmegaB_ji)
    call mem_dealloc(OmegaBsym_ji)
    call mem_dealloc(RHStmp)
    call mem_dealloc(Omega2tmp)
    call array4_close_file(OmegaA,'delete')
    call array4_free(OmegaA)
    call array4_close_file(OmegaB,'delete')
    call array4_free(OmegaB)
    call array4_close_file(Omega2,'keep')
    call array4_close_file(RHS,'keep')

    if (DECinfo%PL>1) call LSTIMER('RES: STEP 2',tcpu,twall,DECinfo%output)


  end subroutine getDoublesResidualMP2_OnFile



  !> \brief Simple double residual for CCSD
  subroutine getDoublesResidualCCSD_simple(omega2,t2,u,gao,aibj,iajb,nocc,nvirt, &
       ppfock,qqfock,xocc,xvirt,yocc,yvirt)

    implicit none
    real(realk) :: aStart,aEnd,bStart,bEnd,cStart,cEnd, &
         dStart,dEnd,eStart,eEnd
    type(array4), intent(inout) :: omega2,t2
    type(array4), intent(inout) :: u,gao,aibj,iajb
    type(array2), intent(inout) :: ppfock, qqfock
    integer, intent(in) :: nocc,nvirt
    type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt
    type(array4) :: abcd, tmp1, X
    type(array4) :: l1, l2, tmp
    type(array2) :: ppX,qqY,pptmp,qqtmp
    integer :: a,i,b,j,k,l,c,d

    aStart=0.0E0_realk; aEnd=0.0E0_realk
    bStart=0.0E0_realk; bEnd=0.0E0_realk
    cStart=0.0E0_realk; cEnd=0.0E0_realk
    dStart=0.0E0_realk; dEnd=0.0E0_realk
    eStart=0.0E0_realk; eEnd=0.0E0_realk


    ! -- A2
    call cpu_time(aStart)
    call array4_add_to(omega2,1.0E0_realk,aibj)

    abcd = get_gmo_simple(gao,xvirt,yvirt,xvirt,yvirt)
    call array4_reorder(t2,[1,3,2,4]) ! -> t2[ab,ij]
    call array4_reorder(abcd,[2,4,1,3])
    tmp1 = array4_init([nvirt,nvirt,nocc,nocc]) ! tmp1[ab,ij]
    call array4_contract2(abcd,t2,tmp1)
    call array4_reorder(tmp1,[1,3,2,4]) ! -> tmp1[ai,bj]
    call array4_add_to(omega2,1.0E0_realk,tmp1)
    call array4_free(tmp1)
    call array4_reorder(t2,[1,3,2,4]) ! -> t2[ai,bj]
    call array4_free(abcd)
    call cpu_time(aEnd)
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: A2 done, norm :',omega2*omega2

    ! -- B2
    call cpu_time(bStart)
    call array4_reorder(iajb,[2,4,1,3]) ! iajb[kc,ld] -> iajb[cd,kl] (iajb[ia,jb] -> iajb[ab,ij])
    call array4_reorder(t2,[1,3,2,4]) ! t2[ci,dj] -> t2[cd,ij] (t2[ai,bj] -> t2[ab,ij])
    tmp1 = array4_init([nocc,nocc,nocc,nocc])
    call array4_contract2(t2,iajb,tmp1) ! tmp1[ij,kl]
    X = get_gmo_simple(gao,xocc,yocc,xocc,yocc) ! X[ki,lj]
    call array4_reorder(X,[2,4,1,3]) ! X[ki,lj] -> X[ij,kl]
    call array4_add_to(X,1.0E0_realk,tmp1)
    call array4_free(tmp1)

    call array4_reorder(X,[3,4,1,2]) ! X[ij,kl] -> X[kl,ij]
    call array4_reorder(t2,[3,4,1,2]) ! t2[ab,kl]-> t2[kl,ab]
    tmp1 = array4_init([nvirt,nvirt,nocc,nocc])
    call array4_contract2(t2,X,tmp1) ! tmp1[ab,ij]
    call array4_reorder(tmp1,[1,3,2,4]) ! tmp1[ab,ij] -> tmp1[ai,bj]
    call array4_add_to(omega2,1.0E0_realk,tmp1)
    call array4_free(X)
    call array4_free(tmp1)
    call cpu_time(bEnd)
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: B2 done, norm :',omega2*omega2

    ! -- C2
    call cpu_time(cStart)

    ! std                                  ij,ab  ->    bi,aj
    call array4_reorder(t2,[4,1,3,2]) ! t2[li,ad] -> t2[dl,ai]
    call array4_reorder(iajb,[1,4,3,2]) ! iajb[dc,kl] -> iajb[dl,kc]
    X = get_gmo_simple(gao,xocc,yocc,xvirt,yvirt) ! X[ki,ac]
    tmp1 = array4_init([nvirt,nocc,nocc,nvirt]) ! tmp1[ai,kc]
    call array4_contract2(t2,iajb,tmp1)
    call array4_reorder(tmp1,[3,2,1,4]) ! -> tmp1[ki,ac]
    call array4_add_to(X,-0.5E0_realk,tmp1)
    call array4_free(tmp1)

    ! a
    call array4_reorder(X,[4,1,3,2]) ! X[ki,ac] -> X[ck,ai]
    tmp1 = array4_init([nvirt,nocc,nvirt,nocc])
    call array4_contract2(X,t2,tmp1) ! tmp1[ai,bj]
    call array4_add_to(omega2,-0.5E0_realk,tmp1)

    ! b
    call array4_reorder(tmp1,[1,4,3,2]) ! tmp1[aj,bi] -> tmp[ai,bj]
    call array4_add_to(omega2,-1.0E0_realk,tmp1)

    ! c
    call array4_reorder(tmp1,[3,2,1,4]) ! tmp1[] -> tmp1[]
    call array4_add_to(omega2,-0.5E0_realk,tmp1)

    ! d
    call array4_reorder(tmp1,[1,4,3,2]) ! tmp[] -> tmp1[]
    call array4_add_to(omega2,-1.0E0_realk,tmp1)

    call array4_free(tmp1)
    call array4_free(X)
    call array4_reorder(t2,[3,2,1,4]) ! t2[dl,ai] -> t2[al,di]
    call cpu_time(cEnd)
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: C2 done, norm :',omega2*omega2

    ! -- D2
    call cpu_time(dStart)
    l1 = getL(gao,xvirt,yocc,xocc,yvirt, &
         xvirt,yvirt,xocc,yocc)
    call array4_reorder(iajb,[3,1,2,4]) ! iajb[dl,kc] -> iajb[kd,lc]
    l2 = getL(iajb)
    X = array4_init([nvirt,nocc,nocc,nvirt]) ! X[ai,kc]
    call array4_reorder(u,[4,3,1,2]) ! u[ai,dl] -> u[ld,ai]
    call array4_contract2(u,l2,X) ! X[ai,kc]
    call array4_scale(X,0.5E0_realk)
    call array4_add_to(X,1.0E0_realk,l1)

    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: l1 norm : ',l1*l1
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: l1 norm : ',l2*l2

    call array4_free(l1)
    call array4_free(l2)

    ! a
    call array4_reorder(X,[3,4,1,2]) ! X[ai,kc] -> X[kc,ai]
    tmp = array4_init([nvirt,nocc,nvirt,nocc])
    call array4_contract2(X,u,tmp)
    call array4_add_to(omega2,0.5E0_realk,tmp)
    call array4_reorder(u,[3,4,2,1]) ! u[ld,ai] -> u[ai,dl]

    ! b
    call array4_reorder(tmp,[3,4,1,2]) ! tmp[bj,ai] -> tmp[ai,bj]
    call array4_add_to(omega2,0.5E0_realk,tmp)
    call array4_free(X)
    call array4_free(tmp)

    call cpu_time(dEnd)
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: D2 done, norm',omega2*omega2

    ! -- E2
    call cpu_time(eStart)
    qqY = array2_init([nvirt,nvirt]) ! qqY[b,c]
    qqtmp = array2_init([nvirt,nvirt])
    ppX = array2_init([nocc,nocc]) ! ppX[k,j]
    pptmp = array2_init([nocc,nocc])

    ! get qqY
    call array4_reorder(u,[4,3,2,1]) ! u[bk,dl] -> u[ldk,b]
    call array4_contract3(u,iajb,qqtmp)
    qqY = array2_add(1.0E0_realk,qqfock,-1.0E0_realk,qqtmp)
    call array2_free(qqtmp)

    ! get ppX
    call array4_reorder(u,[4,3,2,1]) ! u[jdl,c] -> u[cld,j]
    call array4_reorder(iajb,[4,3,2,1]) ! iajb[kdl,c] -> iajb[cld,k]
    call array4_contract3(iajb,u,pptmp)
    ppX = array2_add(1.0E0_realk,ppfock,1.0E0_realk,pptmp)
    call array4_reorder(iajb,[4,3,2,1])

    ! 1
    call array2_transpose(qqY)
    call array4_reorder(t2,[3,4,1,2])
    tmp = array4_init([nvirt,nocc,nvirt,nocc])
    call array4_contract1(t2,qqY,tmp,.true.)
    call array4_reorder(tmp,[3,4,1,2])
    call array4_add_to(omega2,1.0E0_realk,tmp)
    call array4_free(tmp)
    call array4_reorder(t2,[3,4,1,2])

    ! 2
    call array4_contract1(t2,qqY,omega2,.false.)
    call array2_transpose(qqY)

    ! 3
    call array4_reorder(t2,[4,3,2,1])
    tmp = array4_init([nocc,nvirt,nocc,nvirt])
    call array4_contract1(t2,ppX,tmp,.true.)
    call array4_reorder(t2,[4,3,2,1])
    call array4_reorder(tmp,[4,3,2,1])
    call array4_add_to(omega2,-1.0E0_realk,tmp)
    call array4_free(tmp)

    ! 4
    call array4_reorder(t2,[2,1,3,4])
    tmp = array4_init([nocc,nvirt,nvirt,nocc])
    call array4_contract1(t2,ppX,tmp,.true.)
    call array4_reorder(t2,[2,1,3,4])
    call array4_reorder(tmp,[2,1,3,4])
    call array4_add_to(omega2,-1.0E0_realk,tmp)
    call array4_free(tmp)

    call array2_free(ppX)
    call array2_free(qqY)
    call cpu_time(eEnd)

    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: E2 done, norm',omega2*omega2

  end subroutine getDoublesResidualCCSD_simple

  !> \brief Very simple debug version of CCSD residual
  subroutine getDoublesResidual_explicite(omega2,t2,u,gao,aibj,iajb,ppfock,qqfock, &
       xocc,xvirt,yocc,yvirt,nocc,nvirt,nbasis)

    implicit none
    type(array4), intent(inout) :: omega2,gao
    type(array4), intent(in) :: t2,u
    type(array4), intent(in) :: aibj,iajb
    type(array4) :: aikc,acki,tmp
    type(array2), intent(in) :: ppfock, qqfock
    type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt
    type(array2) :: pptmp, qqtmp
    integer, intent(in) :: nocc, nvirt, nbasis
    real(realk) :: starttime,endtime,aStart,aEnd,bStart,bEnd,cStart,cEnd, &
         dStart,dEnd,eStart,eEnd
    integer :: a,b,c,d,i,j,k,l

    type(array4) :: abcd,X,l1,l2

    aStart=0.0E0_realk; aEnd=0.0E0_realk
    bStart=0.0E0_realk; bEnd=0.0E0_realk
    cStart=0.0E0_realk; cEnd=0.0E0_realk
    dStart=0.0E0_realk; dEnd=0.0E0_realk
    eStart=0.0E0_realk; eEnd=0.0E0_realk
    starttime=0.0E0_realk
    endtime=0.0E0_realk

    ! -- A2
    call cpu_time(aStart)
    write(DECinfo%output,'(a)') ' debug explicit :: A2'
    omega2%val = omega2%val + aibj%val
    abcd = get_gmo_simple(gao,xvirt,yvirt,xvirt,yvirt)

    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt

                do d=1,nvirt
                   do c=1,nvirt
                      omega2%val(a,i,b,j) = omega2%val(a,i,b,j) + &
                           t2%val(c,i,d,j)*abcd%val(a,c,b,d)
                   end do
                end do

             end do
          end do
       end do
    end do

    call array4_free(abcd)
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') ' debug explicit :: A2 done, norm',omega2*omega2
    call cpu_time(aEnd)

    ! -- B2
    call cpu_time(bStart)
    write(DECinfo%output,'(a)') ' debug explicit :: B2'
    X = get_gmo_simple(gao,xocc,yocc,xocc,yocc)

    do j=1,nocc
       do l=1,nocc
          do i=1,nocc
             do k=1,nocc

                do d=1,nvirt
                   do c=1,nvirt
                      X%val(k,i,l,j) = X%val(k,i,l,j) + t2%val(c,i,d,j)*iajb%val(k,c,l,d)
                   end do
                end do

             end do
          end do
       end do
    end do

    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt

                do l=1,nocc
                   do k=1,nocc
                      omega2%val(a,i,b,j) = omega2%val(a,i,b,j) + &
                           t2%val(a,k,b,l)*X%val(k,i,l,j)
                   end do
                end do

             end do
          end do
       end do
    end do

    call array4_free(X)
    call cpu_time(bEnd)
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') ' debug explicit :: B2 done, norm',omega2*omega2

    ! -- C2
    call cpu_time(cStart)
    write(DECinfo%output,'(a)') ' debug explicit :: C2'
    X = get_gmo_simple(gao,xocc,yocc,xvirt,yvirt)

    do c=1,nvirt
       do a=1,nvirt
          do i=1,nocc
             do k=1,nocc

                do d=1,nvirt
                   do l=1,nocc
                      X%val(k,i,a,c) = X%val(k,i,a,c) - &
                           0.5E0_realk * t2%val(a,l,d,i)*iajb%val(k,d,l,c)
                   end do
                end do

             end do
          end do
       end do
    end do

    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt

                do c=1,nvirt
                   do k=1,nocc
                      omega2%val(a,i,b,j) = omega2%val(a,i,b,j) - &
                           0.5E0_realk*t2%val(b,k,c,j) * X%val(k,i,a,c)
                   end do
                end do

             end do
          end do
       end do
    end do

    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt

                do c=1,nvirt
                   do k=1,nocc
                      omega2%val(a,i,b,j) = omega2%val(a,i,b,j) - &
                           t2%val(b,k,c,i) * X%val(k,j,a,c)
                   end do
                end do

             end do
          end do
       end do
    end do

    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt

                do c=1,nvirt
                   do k=1,nocc
                      omega2%val(a,i,b,j) = omega2%val(a,i,b,j) - &
                           0.5E0_realk*t2%val(a,k,c,i) * X%val(k,j,b,c)
                   end do
                end do

             end do
          end do
       end do
    end do

    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt

                do c=1,nvirt
                   do k=1,nocc
                      omega2%val(a,i,b,j) = omega2%val(a,i,b,j) - &
                           t2%val(a,k,c,j) * X%val(k,i,b,c)
                   end do
                end do

             end do
          end do
       end do
    end do

    call array4_free(X)

    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') ' debug explicit :: C2 done, norm',omega2*omega2
    call cpu_time(cEnd)

    ! -- D2
    call cpu_time(dStart)
    write(DECinfo%output,'(a)') ' debug explicit :: D2'

    l2 = array4_init([nocc,nvirt,nocc,nvirt])
    do c=1,nvirt
       do k=1,nocc
          do d=1,nvirt
             do l=1,nocc
                l2%val(l,d,k,c) = 2.0E0_realk*iajb%val(l,d,k,c) - iajb%val(l,c,k,d)
             end do
          end do
       end do
    end do

    l1 = array4_init([nvirt,nocc,nocc,nvirt])
    aikc = get_gmo_simple(gao,xvirt,yocc,xocc,yvirt)
    acki = get_gmo_simple(gao,xvirt,yvirt,xocc,yocc)

    do c=1,nvirt
       do k=1,nocc
          do i=1,nocc
             do a=1,nvirt
                l1%val(a,i,k,c) = 2.0E0_realk*aikc%val(a,i,k,c) - acki%val(a,c,k,i)
             end do
          end do
       end do
    end do

    tmp = array4_init([nvirt,nocc,nocc,nvirt])
    do c=1,nvirt
       do k=1,nocc
          do i=1,nocc
             do a=1,nvirt

                do d=1,nvirt
                   do l=1,nocc
                      tmp%val(a,i,k,c) = tmp%val(a,i,k,c) + u%val(a,i,d,l)*l2%val(l,d,k,c)
                   end do
                end do

             end do
          end do
       end do
    end do

    call array4_scale(tmp,0.5E0_realk)
    call array4_add_to(tmp,1.0E0_realk,l1)

    ! standard
    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt

                do c=1,nvirt
                   do k=1,nocc
                      omega2%val(a,i,b,j) = omega2%val(a,i,b,j) + 0.5E0_realk*u%val(b,j,c,k)*tmp%val(a,i,k,c)
                   end do
                end do

             end do
          end do
       end do
    end do

    ! permuted
    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt

                do c=1,nvirt
                   do k=1,nocc
                      omega2%val(a,i,b,j) = omega2%val(a,i,b,j) + 0.5E0_realk*u%val(a,i,c,k)*tmp%val(b,j,k,c)
                   end do
                end do

             end do
          end do
       end do
    end do

    call array4_free(tmp)
    call array4_free(l1)
    call array4_free(l2)

    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') ' debug explicit :: D2 done, norm',omega2*omega2
    call cpu_time(dEnd)

    ! -- E2
    call cpu_time(eStart)
    write(DECinfo%output,'(a)') ' debug explicit :: E2'

    pptmp = array2_init([nocc,nocc])
    qqtmp = array2_init([nvirt,nvirt])

    call array2_add_to(pptmp,1.0E0_realk,ppfock)
    call array2_add_to(qqtmp,1.0E0_realk,qqfock)

    ! do qqtmp
    do b=1,nvirt
       do c=1,nvirt

          do d=1,nvirt
             do k=1,nocc
                do l=1,nocc
                   qqtmp%val(b,c) = qqtmp%val(b,c) - u%val(b,k,d,l)*iajb%val(l,d,k,c)
                end do
             end do
          end do

       end do
    end do

    ! do pptmp
    do k=1,nocc
       do j=1,nocc

          do c=1,nvirt
             do d=1,nvirt
                do l=1,nocc
                   pptmp%val(k,j) = pptmp%val(k,j) + u%val(c,l,d,j)*iajb%val(k,d,l,c)
                end do
             end do
          end do

       end do
    end do

    ! second part of E2
    do a=1,nvirt
       do i=1,nocc
          do b=1,nvirt
             do j=1,nocc

                do c=1,nvirt
                   omega2%val(a,i,b,j) = omega2%val(a,i,b,j) + t2%val(a,i,c,j)*qqtmp%val(b,c)
                end do

             end do
          end do
       end do
    end do

    do a=1,nvirt
       do i=1,nocc
          do b=1,nvirt
             do j=1,nocc

                do c=1,nvirt
                   omega2%val(a,i,b,j) = omega2%val(a,i,b,j) + t2%val(c,i,b,j)*qqtmp%val(a,c)
                end do

             end do
          end do
       end do
    end do

    do a=1,nvirt
       do i=1,nocc
          do b=1,nvirt
             do j=1,nocc

                do k=1,nocc
                   omega2%val(a,i,b,j) = omega2%val(a,i,b,j) - t2%val(a,k,b,j)*pptmp%val(k,i)
                end do

             end do
          end do
       end do
    end do

    do a=1,nvirt
       do i=1,nocc
          do b=1,nvirt
             do j=1,nocc

                do k=1,nocc
                   omega2%val(a,i,b,j) = omega2%val(a,i,b,j) - t2%val(a,i,b,k)*pptmp%val(k,j)
                end do

             end do
          end do
       end do
    end do

    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') ' debug explicit :: C2 done, norm',omega2*omega2
    call cpu_time(eEnd)

  end subroutine getDoublesResidual_explicite


  !> \brief Double residual for CCSD (devel version)
  subroutine getDoublesResidualCCSD_simple2(omega2,t2,u,gao,aibj,iajb,nocc,nvirt, &
       ppfock,qqfock,xocc,xvirt,yocc,yvirt,nbas)

    implicit none
    real(realk) :: starttime,endtime,aStart,aEnd,bStart,bEnd,cStart,cEnd, &
         dStart,dEnd,eStart,eEnd,partial_start,partial_end
    type(array4), intent(inout) :: omega2,t2
    type(array4), intent(inout) :: u,gao,aibj,iajb
    type(array2), intent(inout) :: ppfock, qqfock
    integer, intent(in) :: nocc,nvirt,nbas
    type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt
    type(array4) :: abcd, tmp1, X
    type(array4) :: l1, l2, tmp
    type(array2) :: ppX,qqY,pptmp,qqtmp
    integer :: a,i,b,j,k,l,c,d
    type(array4) :: t2_half_ao,t2_ao,X_half_ao,X_ao, &
         ABe_IJ,IJ_AlBe,Inu_AlBe

    aStart=0.0E0_realk; aEnd=0.0E0_realk
    bStart=0.0E0_realk; bEnd=0.0E0_realk
    cStart=0.0E0_realk; cEnd=0.0E0_realk
    dStart=0.0E0_realk; dEnd=0.0E0_realk
    eStart=0.0E0_realk; eEnd=0.0E0_realk
    starttime=0.0E0_realk; endtime=0.0E0_realk
    partial_start=0.0E0_realk; partial_end=0.0E0_realk

    call cpu_time(starttime)

    ! -- A2
    call cpu_time(aStart)

    ! transform double amplitudes to AO in place of virtual indices
    call cpu_time(partial_start)
    t2_half_ao = array4_init([nbas,nocc,nvirt,nocc]) ! half transformed t2
    call array2_transpose(yvirt)
    call array4_contract1(t2,yvirt,t2_half_ao,.true.) ! t2[ai,bj] -> t2[mu i,bj]
    call array4_reorder(t2_half_ao,[3,1,2,4]) ! t2[mu i,bj] -> t2[b mu,ij]
    t2_ao = array4_init([nbas,nbas,nocc,nocc])
    call array4_contract1(t2_half_ao,yvirt,t2_ao,.true.) ! t2[b mu,ij] -> t2[nu mu,ij]
    call array4_free(t2_half_ao)
    call array2_transpose(yvirt)
    call array4_reorder(t2_ao,[2,1,3,4]) ! t2[nu mu,ij] -> t2[mu nu,ij]
    call cpu_time(partial_end)
    if(DECinfo%PL>1) &
         write(DECinfo%output,'(/,a,f16.3,a)') ' time :: t2[ai,bj] -> t2[mu nu,ij] : ', &
         partial_end-partial_start,' s'

    ! reorder two-electron integrals
    call cpu_time(partial_start)
    call array4_read(gao) ! gao[al mu, be nu]
    !call array4_reorder(gao,[2,4,1,3]) ! gao[al mu,be nu] -> gao[mu nu, al be]
    call array4_reorder(gao,[1,3,2,4])
    call cpu_time(partial_end)
    if(DECinfo%PL>1) &
         write(DECinfo%output,'(a,f16.3,a)') ' time :: sort GAO                  : ', &
         partial_end-partial_start,' s'

    ! contract aplitudes with two-electron integrals in AO basis
    call cpu_time(partial_start)
    X_ao = array4_init([nbas,nbas,nocc,nocc])
    !call array4_contract2(gao,t2_ao,X_ao)
    call array4_contract2_middle(gao,t2_ao,X_ao)
    call array4_free(t2_ao)
    call cpu_time(partial_end)
    if(DECinfo%PL>1) &
         write(DECinfo%output,'(a,f16.3,a)') ' time :: contract with GAO         : ', &
         partial_end-partial_start,' s'

    ! get g[ai,bj] and add it directly to omega2
    call cpu_time(partial_start)
    Inu_AlBe = array4_init([nocc,nbas,nbas,nbas]) ! gao[MuNu,AlBe] -> gao[INu,AlBe]
    call array4_contract1(gao,yocc,Inu_AlBe,.true.)
    call array4_alloc(gao)
    call array4_reorder(Inu_AlBe,[2,1,3,4]) ! Inu_AlBe[INu,AlBe] -> Inu_AlBe[NuI,AlBe]

    IJ_AlBe = array4_init([nocc,nocc,nbas,nbas])
    call array4_contract1(Inu_AlBe,yocc,IJ_AlBe,.true.) ! Inu_AlBe[NuI,AlBe] -> IJ_AlBe[JI,AlBe]
    call array4_free(Inu_AlBe)

    call array4_reorder(IJ_AlBe,[4,3,2,1]) ! IJ_AlBe[JI,AlBe] -> IJ_AlBe[BeAl,IJ]
    ABe_IJ = array4_init([nvirt,nbas,nocc,nocc])
    call array4_contract1(IJ_AlBe,xvirt,ABe_IJ,.true.) ! IJ_AlBe[BeAl,IJ] -> ABe_IJ[ABe,IJ]
    call array4_free(IJ_AlBe)

    call array4_reorder(ABe_IJ,[2,3,1,4]) ! ABe_IJ[BAl,IJ] -> ABe_IJ[AlI,BJ]
    call array4_contract1(ABe_IJ,xvirt,omega2,.false.) ! ABe_IJ[AlI,BJ] -> ABe_IJ[AI,BJ]
    call array4_free(ABe_IJ)

    call cpu_time(partial_end)
    if(DECinfo%PL>1) &
         write(DECinfo%output,'(a,f16.3,a)') ' time :: omega2 += g[ai,bj]        : ', &
         partial_end-partial_start,' s'

    ! transform to mo basis
    call cpu_time(partial_start)
    X_half_ao = array4_init([nvirt,nbas,nocc,nocc])
    call array4_contract1(X_ao,xvirt,X_half_ao,.true.)
    call array4_free(X_ao)
    call array4_reorder(X_half_ao,[2,1,3,4])
    X = array4_init([nvirt,nvirt,nocc,nocc])
    call array4_contract1(X_half_ao,xvirt,X,.true.)
    call array4_free(X_half_ao)
    call array4_reorder(X,[2,3,1,4])
    call array4_add_to(omega2,1.0E0_realk,X)
    call array4_free(X)
    call cpu_time(partial_end)
    if(DECinfo%PL>1) &
         write(DECinfo%output,'(a,f16.3,a)') ' time :: transform result to MO    : ', &
         partial_end-partial_start,' s'

    call cpu_time(aEnd)
    if(DECinfo%PL>1) &
         write(DECinfo%output,'(a,f16.3,a,/)') ' time :: total A term              : ', &
         aEnd-aStart,' s'

    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: A2 done, norm :',omega2*omega2

    ! -- B2
    call cpu_time(bStart)
    call array4_reorder(iajb,[2,4,1,3]) ! iajb[kc,ld] -> iajb[cd,kl] (iajb[ia,jb] -> iajb[ab,ij])
    call array4_reorder(t2,[1,3,2,4]) ! t2[ci,dj] -> t2[cd,ij] (t2[ai,bj] -> t2[ab,ij])
    tmp1 = array4_init([nocc,nocc,nocc,nocc])
    call array4_contract2(t2,iajb,tmp1) ! tmp1[ij,kl]
    X = get_gmo_simple(gao,xocc,yocc,xocc,yocc) ! X[ki,lj]
    call array4_reorder(X,[2,4,1,3]) ! X[ki,lj] -> X[ij,kl]
    call array4_add_to(X,1.0E0_realk,tmp1)
    call array4_free(tmp1)

    call array4_reorder(X,[3,4,1,2]) ! X[ij,kl] -> X[kl,ij]
    call array4_reorder(t2,[3,4,1,2]) ! t2[ab,kl]-> t2[kl,ab]
    tmp1 = array4_init([nvirt,nvirt,nocc,nocc])
    call array4_contract2(t2,X,tmp1) ! tmp1[ab,ij]
    call array4_reorder(tmp1,[1,3,2,4]) ! tmp1[ab,ij] -> tmp1[ai,bj]
    call array4_add_to(omega2,1.0E0_realk,tmp1)
    call array4_free(X)
    call array4_free(tmp1)
    call cpu_time(bEnd)
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: B2 done, norm :',omega2*omega2

    ! -- C2
    call cpu_time(cStart)

    ! std                                  ij,ab  ->    bi,aj
    call array4_reorder(t2,[4,1,3,2]) ! t2[li,ad] -> t2[dl,ai]
    call array4_reorder(iajb,[1,4,3,2]) ! iajb[dc,kl] -> iajb[dl,kc]
    X = get_gmo_simple(gao,xocc,yocc,xvirt,yvirt) ! X[ki,ac]
    tmp1 = array4_init([nvirt,nocc,nocc,nvirt]) ! tmp1[ai,kc]
    call array4_contract2(t2,iajb,tmp1)
    call array4_reorder(tmp1,[3,2,1,4]) ! -> tmp1[ki,ac]
    call array4_add_to(X,-0.5E0_realk,tmp1)
    call array4_free(tmp1)

    ! a
    call array4_reorder(X,[4,1,3,2]) ! X[ki,ac] -> X[ck,ai]
    tmp1 = array4_init([nvirt,nocc,nvirt,nocc])
    call array4_contract2(X,t2,tmp1) ! tmp1[ai,bj]
    call array4_add_to(omega2,-0.5E0_realk,tmp1)

    ! b
    call array4_reorder(tmp1,[1,4,3,2]) ! tmp1[aj,bi] -> tmp[ai,bj]
    call array4_add_to(omega2,-1.0E0_realk,tmp1)

    ! c
    call array4_reorder(tmp1,[3,2,1,4]) ! tmp1[] -> tmp1[]
    call array4_add_to(omega2,-0.5E0_realk,tmp1)

    ! d
    call array4_reorder(tmp1,[1,4,3,2]) ! tmp[] -> tmp1[]
    call array4_add_to(omega2,-1.0E0_realk,tmp1)

    call array4_free(tmp1)
    call array4_free(X)
    call array4_reorder(t2,[3,2,1,4]) ! t2[dl,ai] -> t2[al,di]
    call cpu_time(cEnd)
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: C2 done, norm :',omega2*omega2

    ! -- D2
    call cpu_time(dStart)
    l1 = getL(gao,xvirt,yocc,xocc,yvirt, &
         xvirt,yvirt,xocc,yocc)
    call array4_reorder(iajb,[3,1,2,4]) ! iajb[dl,kc] -> iajb[kd,lc]
    l2 = getL(iajb)
    X = array4_init([nvirt,nocc,nocc,nvirt]) ! X[ai,kc]
    call array4_reorder(u,[4,3,1,2]) ! u[ai,dl] -> u[ld,ai]
    call array4_contract2(u,l2,X) ! X[ai,kc]
    call array4_scale(X,0.5E0_realk)
    call array4_add_to(X,1.0E0_realk,l1)

    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: l1 norm : ',l1*l1
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: l1 norm : ',l2*l2

    call array4_free(l1)
    call array4_free(l2)

    ! a
    call array4_reorder(X,[3,4,1,2]) ! X[ai,kc] -> X[kc,ai]
    tmp = array4_init([nvirt,nocc,nvirt,nocc])
    call array4_contract2(X,u,tmp)
    call array4_add_to(omega2,0.5E0_realk,tmp)
    call array4_reorder(u,[3,4,2,1]) ! u[ld,ai] -> u[ai,dl]

    ! b
    call array4_reorder(tmp,[3,4,1,2]) ! tmp[bj,ai] -> tmp[ai,bj]
    call array4_add_to(omega2,0.5E0_realk,tmp)
    call array4_free(X)
    call array4_free(tmp)

    call cpu_time(dEnd)
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: D2 done, norm',omega2*omega2

    ! -- E2
    call cpu_time(eStart)
    qqY = array2_init([nvirt,nvirt]) ! qqY[b,c]
    qqtmp = array2_init([nvirt,nvirt])
    ppX = array2_init([nocc,nocc]) ! ppX[k,j]
    pptmp = array2_init([nocc,nocc])

    ! get qqY
    call array4_reorder(u,[4,3,2,1]) ! u[bk,dl] -> u[ldk,b]
    call array4_contract3(u,iajb,qqtmp)
    qqY = array2_add(1.0E0_realk,qqfock,-1.0E0_realk,qqtmp)
    call array2_free(qqtmp)

    ! get ppX
    call array4_reorder(u,[4,3,2,1]) ! u[jdl,c] -> u[cld,j]
    call array4_reorder(iajb,[4,3,2,1]) ! iajb[kdl,c] -> iajb[cld,k]
    call array4_contract3(iajb,u,pptmp)
    ppX = array2_add(1.0E0_realk,ppfock,1.0E0_realk,pptmp)
    call array4_reorder(iajb,[4,3,2,1])

    ! 1
    call array2_transpose(qqY)
    call array4_reorder(t2,[3,4,1,2])
    tmp = array4_init([nvirt,nocc,nvirt,nocc])
    call array4_contract1(t2,qqY,tmp,.true.)
    call array4_reorder(tmp,[3,4,1,2])
    call array4_add_to(omega2,1.0E0_realk,tmp)
    call array4_free(tmp)
    call array4_reorder(t2,[3,4,1,2])

    ! 2
    call array4_contract1(t2,qqY,omega2,.false.)
    call array2_transpose(qqY)

    ! 3
    call array4_reorder(t2,[4,3,2,1])
    tmp = array4_init([nocc,nvirt,nocc,nvirt])
    call array4_contract1(t2,ppX,tmp,.true.)
    call array4_reorder(t2,[4,3,2,1])
    call array4_reorder(tmp,[4,3,2,1])
    call array4_add_to(omega2,-1.0E0_realk,tmp)
    call array4_free(tmp)

    ! 4
    call array4_reorder(t2,[2,1,3,4])
    tmp = array4_init([nocc,nvirt,nvirt,nocc])
    call array4_contract1(t2,ppX,tmp,.true.)
    call array4_reorder(t2,[2,1,3,4])
    call array4_reorder(tmp,[2,1,3,4])
    call array4_add_to(omega2,-1.0E0_realk,tmp)
    call array4_free(tmp)

    call array2_free(ppX)
    call array2_free(qqY)
    call cpu_time(eEnd)

    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: E2 done, norm',omega2*omega2

    call cpu_time(endtime)
    if(DECinfo%PL>1) then
       write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD A2 : ',aEnd-aStart,' s'
       write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD B2 : ',bEnd-bStart,' s'
       write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD C2 : ',cEnd-cStart,' s'
       write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD D2 : ',dEnd-dStart,' s'
       write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD E2 : ',eEnd-eStart,' s'
       write(DECinfo%output,'(a,f16.3,a)')  &
            ' time :: CCSD doubles : ',endtime-starttime,' s'
    end if

    return
  end subroutine getDoublesResidualCCSD_simple2


  subroutine get_ccsd_residual_integral_driven_oldarray_wrapper(deltafock,omega2,t2,fock,govov,nocc,nvirt,&
       & ppfock,qqfock,pqfock,qpfock,xocc,xvirt,yocc,yvirt,nbas,MyLsItem, omega1,iter)
       implicit none
       type(array2),intent(in) :: deltafock
       type(array4), intent(inout) :: omega2,t2
       type(array2), intent(inout) :: ppfock, qqfock,pqfock,qpfock,omega1,fock
       type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt
       type(lsitem), intent(inout) :: MyLsItem
       integer,intent(in) :: nbas
       integer,intent(in) :: nocc
       integer,intent(in) :: nvirt
       integer,intent(in) :: iter
       real(realk),target,intent(inout) :: govov(nvirt*nocc*nvirt*nocc)
       real(realk),pointer :: t2_p(:),xo_p(:),xv_p(:),yo_p(:),yv_p(:)
       type(array) :: t2a,ga,o2

      ! get t2 and omega2 into the correct order
      call array4_reorder(t2,[1,3,2,4])
      call ass_D4to1(t2%val,t2_p,t2%dims)
      t2a=array_init(t2%dims,4)
      call memory_deallocate_array_dense(t2a)
      call ass_D4to1(t2%val,t2a%elm1,t2%dims)
      call assoc_ptr_arr(t2a)

      call array4_reorder(omega2,[1,3,2,4])
      o2=array_init(omega2%dims,4)
      call memory_deallocate_array_dense(o2)
      call ass_D4to1(omega2%val,o2%elm1,omega2%dims)
      call assoc_ptr_arr(o2)
  
      ga=array_init([nocc,nvirt,nocc,nvirt],4)
      call memory_deallocate_array_dense(ga)
      ga%elm1 => govov
      call assoc_ptr_arr(ga)

      call ass_D2to1(xocc%val,xo_p,xocc%dims)
      call ass_D2to1(xvirt%val,xv_p,xvirt%dims)
      call ass_D2to1(yocc%val,yo_p,yocc%dims)
      call ass_D2to1(yvirt%val,yv_p,yvirt%dims)
      !call get_ccsd_residual_integral_driven(deltafock%val,omega2%val,t2_p,fock%val,govov,nocc,nvirt,&
      ! & ppfock%val,qqfock%val,pqfock%val,qpfock%val,xo_p,xv_p,yo_p,yv_p,nbas,MyLsItem,&
      ! & omega1%val,iter)
      call get_ccsd_residual_integral_driven(deltafock%val,o2,t2a,fock%val,ga,nocc,nvirt,&
       & ppfock%val,qqfock%val,pqfock%val,qpfock%val,xo_p,xv_p,yo_p,yv_p,nbas,MyLsItem,&
       & omega1%val,iter)

      nullify(t2a%elm1)
      nullify(t2a%elm4)
      call array_free(t2a)
      nullify(o2%elm1)
      nullify(o2%elm4)
      call array_free(o2)
      nullify(ga%elm1)
      nullify(ga%elm4)
      call array_free(ga)

      !reorder ampitudes and residuals for use within the solver
      call array4_reorder(t2,[1,3,2,4])
      call array4_reorder(omega2,[1,3,2,4])
      nullify(xo_p)
      nullify(yo_p)
      nullify(xv_p)
      nullify(yv_p)
  end subroutine get_ccsd_residual_integral_driven_oldarray_wrapper


  !> \brief Get CCSD residual in an integral direct fashion.
  !> \author Patrick Ettenhuber
  !> \date December 2012
  subroutine get_ccsd_residual_integral_driven(deltafock,omega2,t2,fock,govov,no,nv,&
       ppfock,qqfock,pqfock,qpfock,xo,xv,yo,yv,nb,MyLsItem, omega1,iter,rest)
#ifdef VAR_OMP
    use omp_lib,only:omp_get_wtime
#endif
    implicit none
    !> Number of basis functions
    integer,intent(in) :: nb
    !> Number of occupied orbitals
    integer,intent(in) :: no
    !> Number of virtual orbitals
    integer,intent(in) :: nv

    ! derived types needed for the calculation
    !> Long-range correction to Fock matrix
    real(realk),intent(in)    :: deltafock(nb*nb)
    !> the zeroed doubles residual vector on input, on output this contains the
    !full doubles residual
    !real(realk),intent(inout) :: omega2(nv*nv*no*no)
    type(array),intent(inout) :: omega2
    !> the current guess amplitudes
    !real(realk),pointer,intent(in) ::t2(:)
    type(array),intent(inout) ::t2
    !> on output this contains the occupied-occupied block of the t1-fock matrix
    real(realk),intent(inout) :: ppfock(no*no)
    !> on output this contains the virtual-virtual block of the t1-fock matrix
    real(realk),intent(inout) :: qqfock(nv*nv)
    !> on output this contains the occupied-virtual block of the t1-fock matrix
    real(realk),intent(inout) :: pqfock(no*nv)
    !> on output this contains the virtual-occupied block of the t1-fock matrix
    real(realk),intent(inout) :: qpfock(nv*no)
    !> the ao-fock-matrix
    real(realk),intent(inout) :: fock(nb*nb)
    !> zeroed on input, on output this contains the singles residual
    real(realk),intent(inout) :: omega1(no*nv)
    !> on input this contains the transformation matrices for t1-transformed
    !integrals
    real(realk),pointer :: xo(:),xv(:),yo(:),yv(:)

    !> LS item with information needed for integrals
    type(lsitem), intent(inout) :: MyLsItem

    !govov is passed from the dec ccsd driver and is returned,
    !only it is used in another shape here in the subroutine
    !real(realk), dimension(nvirt*nocc*nvirt*nocc) :: govov
    ! for the master iter is the iteration, mdimg = 0
    ! for the slaves iter contains maximum allowed dim alpha
    ! and mdimg maximum allowed dim gamma
    integer,intent(in) :: iter
    !real(realk),intent(inout) :: govov(nv*no*nv*no)
    type(array),intent(inout) :: govov
    !logical that specifies whether the amplitudes were read
    logical, optional, intent(inout) :: rest

    ! elementary types needed for the calculation
    real(realk), pointer :: w0(:), w1(:), w2(:), w3(:),Had(:),t2_d(:,:,:,:),&
                            &Gbi(:),uigcj(:), tGammadij(:),tpl(:),tmi(:),sio4(:),&
                            &gvvoo(:),gvoov(:),gvvoo_r(:),gvoov_r(:)

    type(matrix) :: Dens,iFock
    ! Variables for mpi
    logical :: master
    integer :: fintel,nintel,fe,ne
    integer(kind=ls_mpik) :: nnod
    real(realk) :: startt, stopp
    
    type(array) :: u2
    type(array) :: gvoova,gvvooa
    !special arrays for scheme=1
    type(array) :: t2jabi,u2kcjb
    integer,pointer :: mpi_task_distribution(:)
    integer(kind=ls_mpik) :: win_in_g
#ifdef VAR_LSMPI
    ! stuff for direct communication
    type(c_ptr) :: gvvoo_c,gvoov_c,sio4_c,gvvoo_p,gvoov_p
    integer(kind=ls_mpik) :: gvvoo_w, gvoov_w, sio4_w
    integer(kind=ls_mpik) :: ierr, hstatus, nctr
    integer :: rcnt(infpar%lg_nodtot),dsp(infpar%lg_nodtot)
    character*(MPI_MAX_PROCESSOR_NAME) :: hname
    real(realk),pointer :: mpi_stuff(:)
    type(c_ptr) :: mpi_ctasks
    !integer(kind=ls_mpik),pointer :: win_in_g(:)
#endif

    ! CHECKING and MEASURING variables
    integer(kind=long) :: maxsize64,dummy64
    integer :: myload,double_2G_nel,nelms,n4
    real(realk) :: tcpu, twall,tcpu1,twall1,tcpu2,twall2, deb1,deb2,MemFree,ActuallyUsed
    real(realk) :: tcpu_start,twall_start, tcpu_end,twall_end,time_a, time_c, time_d,time_singles
    real(realk) :: time_doubles,time_start,timewall_start,wait_time,max_wait_time
    integer     :: scheme

    ! variables used for BATCH construction and INTEGRAL calculation
    integer :: alphaB,gammaB,dimAlpha,dimGamma
    integer :: dim1,dim2,dim3,K,MinAObatch
    integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
    integer :: iorb,nthreads,idx
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    Character(80)        :: FilenameCS,FilenamePS
    Character(80),pointer:: BatchfilenamesCS(:,:)
    Character(80),pointer:: BatchfilenamesPS(:,:)
    Character            :: INTSPEC(5)
    logical :: FoundInMem,fullRHS, doscreen
    integer :: MaxAllowedDimAlpha,MaxActualDimAlpha,nbatchesAlpha,nbatches
    integer :: MaxAllowedDimGamma,MaxActualDimGamma,nbatchesGamma
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
    integer :: a,b,i,j,l,m,n,c,d,fa,fg,la,lg,worksize
    integer :: nb2,nb3,nv2,no2,b2v,o2v,v2o,no3
    integer(kind=8) :: nb4,o2v2,no4
    integer :: tlen,tred,nor,nvr,goffs,aoffs
    integer :: prev_alphaB,mpi_buf
    logical :: jobtodo,first_round,dynamic_load,restart,print_debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !TEST AND DEVELOPMENT VARIABLES!!!!!
    real(realk) :: op_start,op_stop
    integer :: testmode(4)
    integer(kind=long) :: xyz,zyx1,zyx2
    logical :: debug
    character(ARR_MSG_LEN) :: msg
    integer :: def_atype
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef VAR_OMP
    integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
    TYPE(DECscreenITEM)    :: DecScreen
    restart=.false.
    if(present(rest))restart=rest
    master = .true.
    def_atype = DENSE
    scheme = 0
    dynamic_load = DECinfo%dyn_load
    startt=0.0E0_realk
    stopp=0.0E0_realk
    !double_2G_nel=250000000
    double_2G_nel=170000000
    print_debug = (DECinfo%PL>1)
    call print_norm(t2)



!HACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print_debug = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




#ifdef VAR_DEBUG
    double_2G_nel=20
    print_debug = .true.
#endif
    

    ! Set integral info
    ! *****************
    INTSPEC(1)='R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)='R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)='R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)='R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)='C' !C = Coulomb operator
#ifdef VAR_LSMPI
    master=(infpar%lg_mynum == 0)
    nnod=infpar%lg_nodtot
    call get_int_dist_info(int(no*no*nv*nv,kind=8),fintel,nintel)
    def_atype=TILED_DIST
#endif

    ! Some timings
    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu_start,twall_start,DECinfo%output)
    call LSTIMER('START',time_start,timewall_start,DECinfo%output)

    doscreen = MyLsItem%setting%scheme%cs_screen.OR.MyLsItem%setting%scheme%ps_screen
    debug = .false.

    ! Set some shorthand notations
    nb2 = nb*nb
    nb3 = nb2*nb
    nb4 = int(nb3*nb,kind=long)
    nv2 = nv*nv
    no2 = no*no
    no3 = no2*no
    no4 = int(no3*no,kind=long)
    b2v = nb2*nv
    o2v = no2*nv
    v2o = nv2*no
    o2v2 = int(nv2*no2,kind=long)
    nor=no*(no+1)/2
    nvr=nv*(nv+1)/2
!print*,"HACK:random t"
!if(master)call random_number(t2%elm1)

    ! Initialize stuff
    nullify(orb2batchAlpha)
    nullify(batchdimAlpha)
    nullify(batchsizeAlpha)
    nullify(batch2orbAlpha)
    nullify(batchindexAlpha)
    nullify(orb2batchGamma)
    nullify(batchdimGamma)
    nullify(batchsizeGamma)
    nullify(batch2orbGamma)
    nullify(batchindexGamma)
    nullify(w0)
    nullify(w1)
    nullify(w2)
    nullify(w3)
    nullify(Had)
    nullify(Gbi)
    nullify(uigcj)
    nullify(tGammadij)
    nullify(tpl)
    nullify(tmi)
    nullify(sio4)
    nullify(gvoov)
    nullify(gvvoo)
    nullify(gvoov_r)
    nullify(gvvoo_r)
#ifdef VAR_LSMPI
    nullify(mpi_task_distribution)
#endif

    if(master) then
    !==================================================
    !                  Batch construction             !
    !==================================================


    ! Get free memory and determine maximum batch sizes
    ! -------------------------------------------------
      call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch)
      call get_currently_available_memory(MemFree)
      call get_max_batch_sizes(scheme,nb,nv,no,MaxAllowedDimAlpha,MaxAllowedDimGamma,&
                                     &MinAObatch,DECinfo%manual_batchsizes,iter,MemFree,.true.)

      !SOME WORDS ABOUT THE CHOSEN SCHEME:
      ! Depending on the availability of memory on the nodes a certain scheme
      ! for the calculations is chosen, hereby the schemes 4-1 are the MPI
      ! schemes with decreasing memory reqirements. Hereby scheme 4 and 0 can be
      ! used without MPI. Scheme 0 should never be chosen with MPI and
      ! schemes 1-3 should never be chosen without MPI

      ! scheme 4: everything is treated locally, only the main integral driven
      !           loop is MPI-parallel
      ! scheme 3: treat govov, gvoov and gvvoo in the main part in PDM
      ! scheme 2: additionally to 3 also the amplitudes, u, the residual are
      !           treated in PDM, the strategy is to only use one V^2O^2 in 
      !           local mem
      ! scheme 1: is the high scaling scheme with the same constraints as in 2,
      !           this reduces the communication (probably not even worth 
      !           implementing)
      ! scheme 0: the "tradional" high scaling scheme for non-MPI calculations

#ifndef VAR_LSMPI
      !scheme 1 is a pure mpi-scheme, it is not selected by get_max_batch_size,
      !but may be chosen by the user (for debugging reasons). Therefore we have
      !to quit here
      if(scheme==3.or.scheme==2.or.scheme==1) call lsquit("ERROR(ccsd_residual_integral_driven):wrong choice of scheme",-1)
#else
      !allocate the dense part of the arrays if all can be kept in local memory.
      !do that for master only, as the slaves recieve the data via StartUpSlaves
      if((scheme==4.or.scheme==0).and.govov%atype/=DENSE)then
        if(iter==1) call memory_allocate_array_dense(govov)
        if(iter/=1) call array_cp_tiled2dense(govov,.false.)
      endif 
#endif
    endif


!all communication for MPI prior to the loop
#ifdef VAR_LSMPI
    StartUpSlaves: if(master .and. infpar%lg_nodtot>1) then
      call ls_mpibcast(CCSDDATA,infpar%master,infpar%lg_comm)
      call mpi_communicate_ccsd_calcdata(omega2,t2,govov,xo,xv,yo,&
      &yv,MyLsItem,nb,nv,no,iter,scheme,DECinfo%solver_par)
    endif StartUpSlaves
      
    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(scheme,infpar%master)
    call ls_mpi_buffer(dynamic_load,infpar%master)
    call ls_mpi_buffer(restart,infpar%master)
    call ls_mpi_buffer(MaxAllowedDimAlpha,infpar%master)
    call ls_mpi_buffer(MaxAllowedDimGamma,infpar%master)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

    hstatus = 80
    CALL MPI_GET_PROCESSOR_NAME(hname,hstatus,ierr)

    !dense part was allocated in the communicate subroutine
    if(scheme==4.or.scheme==0)then
      govov%atype    = DENSE
    endif
    govov%init_type  = ALL_INIT
    omega2%init_type = ALL_INIT
    t2%init_type     = ALL_INIT
#endif

    !if the residual is handeled as dense, allocate and zero it, adjust the
    !access parameters to the data
    if(omega2%atype/=DENSE.and.(scheme==0.or.scheme==3.or.scheme==4))then
      call memory_allocate_array_dense(omega2)
      omega2%atype=DENSE
    endif
    call array_zero(omega2)

    !ZERO the integral matrix if  first iteration
    if(iter==1.and.(scheme==1.or.scheme==2.or.scheme==3))then
      do i=1,govov%nlti
        govov%ti(i)%t=0.0E0_realk
      enddo
    elseif(iter==1.and.(scheme==0.or.scheme==4))then
      govov%elm1=0.0E0_realk
    endif

    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************

    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchGamma,nb)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
         & nb,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,&
         &nbatchesGamma,orb2BatchGamma)
    if(master)write(DECinfo%output,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma,&
                                       & 'with maximum size',MaxActualDimGamma

    ! Translate batchindex to orbital index
    ! -------------------------------------
    call mem_alloc(batch2orbGamma,nbatchesGamma)
    do idx=1,nbatchesGamma
       call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx))
       batch2orbGamma(idx)%orbindex = 0
       batch2orbGamma(idx)%norbindex = 0
    end do
    do iorb=1,nb
       idx = orb2batchGamma(iorb)
       batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
       K = batch2orbGamma(idx)%norbindex
       batch2orbGamma(idx)%orbindex(K) = iorb
    end do


    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************

    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchAlpha,nb)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
         & nb,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha)
    if(master)write(DECinfo%output,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha&
                                      &, 'with maximum size',MaxActualDimAlpha

    ! Translate batchindex to orbital index
    ! -------------------------------------
    call mem_alloc(batch2orbAlpha,nbatchesAlpha)
    do idx=1,nbatchesAlpha
       call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
       batch2orbAlpha(idx)%orbindex = 0
       batch2orbAlpha(idx)%norbindex = 0
    end do
    do iorb=1,nb
       idx = orb2batchAlpha(iorb)
       batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
       K = batch2orbAlpha(idx)%norbindex
       batch2orbAlpha(idx)%orbindex(K) = iorb
    end do



    ! ************************************************
    ! *  Allocate matrices used in the batched loop  *
    ! ************************************************


    !PRINT some information
    if(master) then
      if(scheme==4) write(DECinfo%output,'("Using memory intensive scheme (NON-PDM)")')
      if(scheme==3) write(DECinfo%output,'("Using memory intensive scheme with direct updates")')
      if(scheme==2) write(DECinfo%output,'("Using memory intensive scheme only 1x V^2O^2")')
      if(scheme==1) write(DECinfo%output,'("Using memory saving scheme with direct updates")')
      if(scheme==0) write(DECinfo%output,'("Using memory saving scheme (NON-PDM), this is only for non-MPI")')
      ActuallyUsed=get_min_mem_req(no,nv,nb,MaxActualDimAlpha,MaxActualDimGamma,3,scheme,.false.)
      write(DECinfo%output,'("Using",1f8.4,"% of available Memory in part B on master")')ActuallyUsed/MemFree*100
      ActuallyUsed=get_min_mem_req(no,nv,nb,MaxActualDimAlpha,MaxActualDimGamma,2,scheme,.false.)
      write(DECinfo%output,'("Using",1f8.4,"% of available Memory in part C on master")')ActuallyUsed/MemFree*100
    endif
    

    !get everything that needs to be gained before t2 is put into pdm
    call mem_alloc(tpl,int(nor*nvr,kind=long))
    call mem_alloc(tmi,int(nor*nvr,kind=long))
    call get_tpl_and_tmi(t2%elm1,nv,no,tpl,tmi)
    !call mem_alloc(u2,no2*nv2)
    if(scheme==2.or.scheme==1)then
      u2=array_init([nv,nv,no,no],4,TILED_DIST,ALL_INIT)
      call array_zero(u2)
      if(master)then 
        call array_add(u2,2.0E0_realk,t2%elm1,[2,1,3,4])
        call array_add(u2,-1.0E0_realk,t2%elm1,[2,1,4,3])
      endif
      if(scheme==1)then
        !fill the reordered amplitude, t2  abij -> jabi [4,1,2,3]
        t2jabi=array_init([no,nv,nv,no],4,TILED_DIST,ALL_INIT)
        call array_convert(t2%elm1,t2jabi,[4,1,2,3])
      endif
      call array_mv_dense2tiled(t2,.true.)
      if(scheme==1)then
        !fill the reordered u2 array used in d2 direct term
        u2kcjb=array_init([no,nv,no,nv],4,TILED_DIST,ALL_INIT)
        call array_cp_tiled2dense(u2,.false.)
        !dcij -> jdic [4,1,3,2]
        call array_convert(u2%elm1,u2kcjb,[4,1,3,2])
        call memory_deallocate_array_dense(u2)
      endif
    else
      u2=array_init([nv,nv,no,no],4)
      !calculate u matrix: t[c d i j] -> t[d c i j], 2t[d c i j] - t[d c j i] = u [d c i j]
      call array_reorder_4d(2.0E0_realk,t2%elm1,nv,nv,no,no,[2,1,3,4],0.0E0_realk,u2%elm1)
      call array_reorder_4d(-1.0E0_realk,t2%elm1,nv,nv,no,no,[2,1,4,3],1.0E0_realk,u2%elm1)
    endif

    call mem_alloc(Had,nv*nb)
    call mem_alloc(Gbi,nb*no)
    if(DECinfo%ccModel>2)then
#ifdef VAR_LSMPI
      call mem_alloc(sio4,sio4_c,int(nor*no2,kind=long))
      call lsmpi_win_create(sio4,sio4_w,int(nor*no2,kind=long),infpar%lg_comm)
#else
      call mem_alloc(sio4,nor*no2)
#endif
      if(scheme==4)then
#ifdef VAR_LSMPI
        call mem_alloc(gvoov,gvoov_c,int(no2*nv2,kind=long))
        call mem_alloc(gvvoo,gvvoo_c,int(no2*nv2,kind=long))
        call lsmpi_win_create(gvoov,gvoov_w,int(no2*nv2,kind=long),infpar%lg_comm)
        call lsmpi_win_create(gvvoo,gvvoo_w,int(no2*nv2,kind=long),infpar%lg_comm)
#else
        call mem_alloc(gvoov,no2*nv2)
        call mem_alloc(gvvoo,no2*nv2)
#endif
        gvoov=0.0E0_realk
        gvvoo=0.0E0_realk
      elseif(scheme==3)then
#ifdef VAR_LSMPI
        call mem_alloc(gvoov_r,gvoov_p,int(nintel,kind=long))
        call mem_alloc(gvvoo_r,gvvoo_p,int(nintel,kind=long))
        gvvoo_r=0.0E0_realk
        gvoov_r=0.0E0_realk
        call lsmpi_win_create(gvoov_r,gvoov_w,int(nintel,kind=long),infpar%lg_comm)
        call lsmpi_win_create(gvvoo_r,gvvoo_w,int(nintel,kind=long),infpar%lg_comm)
#else
        call mem_alloc(gvoov,int(no2*nv2,kind=long))
        call mem_alloc(gvvoo,int(no2*nv2,kind=long))
        gvvoo=0.0E0_realk
        gvoov=0.0E0_realk
#endif
      elseif(scheme==2)then
        gvvooa=array_init([nv,no,no,nv],4,TILED_DIST,ALL_INIT)
        gvoova=array_init([nv,no,nv,no],4,TILED_DIST,ALL_INIT)
        call array_zero(gvvooa)
        call array_zero(gvoova)
      endif
    endif
    !zero the matrix
    Had=0.0E0_realk
    Gbi=0.0E0_realk
   
    ! allocate working arrays depending on the batch sizes
    maxsize64 = nb2*MaxActualDimAlpha*MaxActualDimGamma
    call mem_alloc(w0,maxsize64)

    maxsize64 = max(int(nb2*MaxActualDimAlpha*MaxActualDimGamma,kind=8),int(v2o*MaxActualDimAlpha,kind=8))
    maxsize64 = max(maxsize64,int(o2v*MaxActualDimGamma,kind=8))
    if(scheme==4.or.scheme==3) maxsize64 = max(maxsize64,int(o2v*MaxActualDimAlpha,kind=8))
    call mem_alloc(w1,maxsize64)

    maxsize64 = max(int(nb*nb*MaxActualDimAlpha*MaxActualDimGamma,kind=8),int(nv2*no2,kind=8))
    maxsize64 = max(maxsize64,int(nor*no2,kind=8))
    if(scheme==0) maxsize64=max(maxsize64,int(o2v*MaxActualDimAlpha,kind=8))
    call mem_alloc(w2,maxsize64)

    maxsize64 = max(int(nv*no*MaxActualDimAlpha*MaxActualDimGamma,kind=8),int(no2*MaxActualDimAlpha*MaxActualDimGamma,kind=8))
    maxsize64 = max(maxsize64,int(o2v*MaxActualDimAlpha,kind=8))
    maxsize64 = max(maxsize64,int(2*nor*MaxActualDimAlpha*MaxActualDimGamma,kind=8)) 
    maxsize64 = max(maxsize64,int(nor*nv*MaxActualDimAlpha,kind=8)) 
    maxsize64 = max(maxsize64,int(nor*nv*MaxActualDimGamma,kind=8)) 
    maxsize64 = max(maxsize64,int(no*nor*MaxActualDimAlpha,kind=8)) 
    maxsize64 = max(maxsize64,int(no*nor*MaxActualDimGamma,kind=8)) 
    call mem_alloc(w3,maxsize64)

    !allocate semi-permanent storage arrays for loop
    !print *,"allocing help things:",o2v*MaxActualDimGamma*2,&
    !      &(8.0E0_realk*o2v*MaxActualDimGamma*2)/(1024.0E0_realk*1024.0E0_realk*1024.0E0_realk)
    call mem_alloc(uigcj,o2v*MaxActualDimGamma)
    if(DECinfo%ccModel>2.and.(scheme==0.or.scheme==1)) call mem_alloc(tGammaDij,o2v*MaxActualDimGamma)

    w1=0.0E0_realk
    w2=0.0E0_realk
    w3=0.0E0_realk
    w0=0.0E0_realk
    uigcj=0.0E0_realk
    if(DECinfo%ccModel>2)then
      if(scheme==0.or.scheme==1)then
        tGammadij=0.0E0_realk
      endif
      sio4=0.0E0_realk
    endif


    ! This subroutine builds the full screening matrix.
    call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
         & nbatches,nbatchesAlpha,nbatchesGamma,INTSPEC)
    IF(mylsitem%setting%scheme%cs_screen .OR. &
         & mylsitem%setting%scheme%ps_screen)THEN
       call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
            & nb,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,DECinfo%output,DECinfo%output)
       call II_getBatchOrbitalScreenK(DecScreen,mylsitem%setting,&
            & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
            & batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,DECinfo%output,DECinfo%output)
    ENDIF
    !setup LHS screening - the full AO basis is used so we can use the
    !                      full matrices:        FilenameCS and FilenamePS
    !Note that it is faster to calculate the integrals in the form
    !(dimAlpha,dimGamma,nbasis,nbasis) so the full AO basis is used on the RHS
    !but the integrals is stored and returned in (nbasis,nbasis,dimAlpha,dimGamma)

#ifdef VAR_OMP
    nthreads=OMP_GET_MAX_THREADS()
    if(master)write(DECinfo%output,*) 'Starting DEC-CCSD residuals - OMP. Number of threads: ', OMP_GET_MAX_THREADS()
#else
    nthreads=1
    if(master)write(DECinfo%output,*) 'Starting DEC-CCSD integral/amplitudes - NO OMP!'
#endif

#ifdef VAR_LSMPI
    if(.not.dynamic_load)then
      ! Calculate the batches for a good load balance
      call mem_alloc(mpi_task_distribution,nbatchesAlpha*nbatchesGamma)
      mpi_task_distribution = 0
      myload = 0
      
      call distribute_mpi_jobs(mpi_task_distribution,nbatchesAlpha,nbatchesGamma,batchdimAlpha,&
           &batchdimGamma,myload,scheme,no,nv,nb,batch2orbAlpha,batch2orbGamma)
    else
      !call mem_alloc(win_in_g,nbatchesGamma)
      !call mem_alloc(mpi_task_distribution,mpi_ctasks,nbatchesGamma) 
      call mem_alloc(mpi_stuff,mpi_ctasks,nbatchesGamma) 
      mpi_stuff=0.0E0_realk
      if(master) mpi_stuff(1) = float(infpar%lg_nodtot)
      call lsmpi_win_create(mpi_stuff,win_in_g,nbatchesGamma,infpar%lg_comm)
    endif
    !startt = omp_get_wtime()
#endif
    myload = 0



    if(master)call LSTIMER('CCSD part A',time_start,timewall_start,DECinfo%output)

    fullRHS=(nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)


    !**********************************
    ! Begin the loop over gamma batches
    !**********************************

    first_round=.false.
    if(dynamic_load)first_round=.true.

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
       dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
       GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
       GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
       !short hand notation
       fg = GammaStart
       lg = dimGamma

       !Lambda^h [gamma d] u[d c i j] = u [gamma c i j]
       if(scheme==2.or.scheme==1)then
         !call array_cp_tiled2dense(u2,.false.)
         call array_convert(u2,w2)
         call dgemm('n','n',lg,o2v,nv,1.0E0_realk,yv(fg),nb,w2,nv,0.0E0_realk,w1,lg)
         !call memory_deallocate_array_dense(u2)
       else
         call dgemm('n','n',lg,o2v,nv,1.0E0_realk,yv(fg),nb,u2%elm1,nv,0.0E0_realk,w1,lg)
       endif
       !u [gamma c i j ] -> u [i gamma c j]
       call array_reorder_4d(1.0E0_realk,w1,lg,nv,no,no,[3,1,2,4],0.0E0_realk,uigcj)
       call print_norm(uigcj,int(lg*nv*no*no,kind=8))

       !Lambda^h [gamma c] * t [c d i j] = t [gamma d i j]
       if(DECinfo%ccModel>2.and.(scheme==0.or.scheme==1)) then
         if(scheme==1)then
           call array_convert(t2,w2)
           call dgemm('n','n',lg,o2v,nv,1.0E0_realk,yv(fg),nb,w2,nv,0.0E0_realk,tGammadij,lg)
         else
           call dgemm('n','n',lg,o2v,nv,1.0E0_realk,yv(fg),nb,t2%elm1,nv,0.0E0_realk,tGammadij,lg)
         endif
       endif

       alphaB=0
       
    !**********************************
    ! Begin the loop over alpha batches
    !**********************************

    !BatchAlpha: do alphaB = 1,nbatchesAlpha    ! AO batches
    BatchAlpha: do while(alphaB<=nbatchesAlpha) ! AO batches
      
      !check if the current job is to be done by current node
      call check_job(scheme,first_round,dynamic_load,alphaB,gammaB,nbatchesAlpha,&
        &nbatchesGamma,mpi_task_distribution,win_in_g,print_debug)
       !break the loop if alpha become too large, necessary to account for all
       !of the mpi and non mpi schemes, this is accounted for, because static,
       !and dynamic load balancing are enabled
       if(alphaB>nbatchesAlpha)exit

       dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
       AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
       AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch

       !short hand notation
       fa = AlphaStart
       la = dimAlpha
       myload=myload+la*lg
       !u[k gamma  c j] * Lambda^p [alpha j] ^T = u [k gamma c alpha]
       call dgemm('n','t',no*nv*lg,la,no,1.0E0_realk,uigcj,no*nv*lg,xo(fa),nb,0.0E0_realk,w1,nv*no*lg)
       call lsmpi_poke()
       !Transpose u[k gamma c alpha]^T -> u[c alpha k gamma]
       !call mat_transpose(w1,no*lg, nv*la,w3)
       call array_reorder_4d(1.0E0_realk,w1,no,lg, nv,la,[3,4,1,2],0.0E0_realk,w3)
       call lsmpi_poke()

       !print*,"GAMMA:",fg,nbatchesGamma,"ALPHA:",fa,nbatchesAlpha
       !print*,"--------------------------------------------------"

       !setup RHS screening - here we only have a set of AO basisfunctions
       !                      so we use the batchscreening matrices.
       !                      like BatchfilenamesCS(alphaB,gammaB)
       !Note that it is faster to calculate the integrals in the form
       !(dimAlpha,dimGamma,nbasis,nbasis) so the subset of the AO basis is used on the LHS
       !but the integrals is stored and returned in (nbasis,nbasis,dimAlpha,dimGamma)
       IF(doscreen) Mylsitem%setting%LST_GAB_RHS => DECSCREEN%masterGabLHS
       IF(doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p
       ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
       ! ************************************************************************************
       dim1 = nb*nb*dimAlpha*dimGamma   ! dimension for integral array
       ! Store integral in tmp1(1:dim1) array in (beta,delta,alphaB,gammaB) order
       call LSTIMER('START',tcpu1,twall1,DECinfo%output)
       !Mylsitem%setting%scheme%intprint=6
       call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, Mylsitem%setting, w1,batchindexAlpha(alphaB),&
            &batchindexGamma(gammaB),&
            &batchsizeAlpha(alphaB),batchsizeGamma(gammaB),nb,nb,dimAlpha,dimGamma,fullRHS,nbatches,INTSPEC)
       call lsmpi_poke()
       !Mylsitem%setting%scheme%intprint=0
       call LSTIMER('START',tcpu2,twall2,DECinfo%output)
      
       !if(master)call LSTIMER('INTEGRAL1',time_start,timewall_start,DECinfo%output)
       call array_reorder_4d(1.0E0_realk,w1,nb,nb,la,lg,[4,2,3,1],0.0E0_realk,w0)
       call lsmpi_poke()

       ! I [gamma delta alpha beta] * Lambda^p [beta l] = I[gamma delta alpha l]
       call dgemm('n','n',lg*la*nb,no,nb,1.0E0_realk,w0,lg*nb*la,xo,nb,0.0E0_realk,w2,lg*nb*la)
       call lsmpi_poke()
       !Transpose I [gamma delta alpha l]^T -> I [alpha l gamma delta]
       !call mat_transpose(w2,lg*nb,la*no,w1)
       call array_reorder_4d(1.0E0_realk,w2,lg,nb,la,no,[3,4,1,2],0.0E0_realk,w1)
       call lsmpi_poke()


       !u [b alpha k gamma] * I [alpha k gamma delta] =+ Had [a delta]
       call dgemm('n','n',nv,nb,lg*la*no,1.0E0_realk,w3,nv,w1,lg*la*no,1.0E0_realk,Had,nv)
       call lsmpi_poke()
       call print_norm(Had,int(nb*nv,kind=8))

       !VVOO
       if (DECinfo%ccModel>2.and.(scheme==4.or.scheme==3.or.scheme==2)) then
        !I [alpha  i gamma delta] * Lambda^h [delta j]          = I [alpha i gamma j]
        call dgemm('n','n',la*no*lg,no,nb,1.0E0_realk,w1,la*no*lg,yo,nb,0.0E0_realk,w3,la*no*lg)
        call lsmpi_poke()
        ! gvvoo = (vv|oo) constructed from w2                 = I [alpha i j  gamma]
        call array_reorder_4d(1.0E0_realk,w3,la,no,lg,no,[1,2,4,3],0.0E0_realk,w2)
        call lsmpi_poke()
        !I [alpha  i j gamma] * Lambda^h [gamma b]            = I [alpha i j b]
        call dgemm('n','n',la*no2,nv,lg,1.0E0_realk,w2,la*no2,yv(fg),nb,0.0E0_realk,w3,la*no2)
        call lsmpi_poke()
        !Lambda^p [alpha a]^T * I [alpha i j b]             =+ gvvoo [a i j b]
        if(scheme==4)then
          call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w3,la,1.0E0_realk,gvvoo,nv)
        elseif(scheme==3)then
#ifdef VAR_LSMPI
          call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w3,la,0.0E0_realk,w2,nv)
          call dist_int_contributions(w2,int(no2*nv2,kind=long),gvvoo_w)
#else
          call lsquit("ERROR(ccsd_integral_driven):this should never appear gvvoo_w",-1)
#endif
        elseif(scheme==2)then
          call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w3,la,0.0E0_realk,w2,nv)
          call array_add(gvvooa,1.0E0_realk,w2)
        endif
        call lsmpi_poke()
       endif
       call print_norm(gvvoo,o2v2)

       ! I [alpha l gamma delta] * Lambda^h [delta c] = I[alpha l gamma c]
       call dgemm('n','n',lg*la*no,nv,nb,1.0E0_realk,w1,la*no*lg,yv,nb,0.0E0_realk,w3,la*no*lg)
       call lsmpi_poke()
       !I [alpha l gamma c] * u [l gamma c j]  =+ Gbi [alpha j]
       call dgemm('n','n',la,no,nv*no*lg,1.0E0_realk,w3,la,uigcj,nv*no*lg,1.0E0_realk,Gbi(fa),nb)
       call lsmpi_poke()
       call print_norm(Gbi,int(nb*no,kind=8))
       
       
       !CALCULATE govov FOR ENERGY
       !Reorder I [alpha j gamma b]                      -> I [alpha j b gamma]
       if(iter==1.or.(scheme==4.or.scheme==3.or.scheme==2))&
       &call array_reorder_4d(1.0E0_realk,w3,la,no,lg,nv,[1,2,4,3],0.0E0_realk,w2)
       
       if(iter==1)then
         !I [alpha  j b gamma] * Lambda^h [gamma a]          = I [alpha j b a]
         call dgemm('n','n',la*no*nv,nv,lg,1.0E0_realk,w2,la*no*nv,yv(fg),nb,0.0E0_realk,w1,la*no*nv)
         call lsmpi_poke()
         !Lambda^p [alpha i]^T * I [alpha j b a]             =+ govov [i j b a]
         if(scheme==4)then
           call dgemm('t','n',no,v2o,la,1.0E0_realk,xo(fa),nb,w1,la,1.0E0_realk,govov%elm1,no)
           !call array_add(govov,1.0E0_realk,w2,no2*nv2)
         else
           ! i a j b
           call dgemm('t','n',no,v2o,la,1.0E0_realk,xo(fa),nb,w1,la,0.0E0_realk,w2,no)
           call array_add(govov,1.0E0_realk,w2,[1,4,2,3])
         endif
         call lsmpi_poke()
       endif
       call print_norm(govov)

       !VOOV
       if((restart.and.iter==1).and..not.scheme==4)&
         &call array_reorder_4d(1.0E0_realk,w3,la,no,lg,nv,[1,2,4,3],0.0E0_realk,w2)
       if (DECinfo%ccModel>2.and.(scheme==4.or.scheme==3.or.scheme==2).and.(iter/=1.or.restart)) then
        ! gvoov = (vo|ov) constructed from w2               = I [alpha j b  gamma]
        !I [alpha  j b gamma] * Lambda^h [gamma i]          = I [alpha j b i]
        call dgemm('n','n',la*no*nv,no,lg,1.0E0_realk,w2,la*no*nv,yo(fg),nb,0.0E0_realk,w1,la*no*nv)
        call lsmpi_poke()
        !Lambda^p [alpha a]^T * I [alpha j b i]             =+ gvoov [a j b i]
        if(scheme==4)then
          call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w1,la,1.0E0_realk,gvoov,nv)
        elseif(scheme==3)then
#ifdef VAR_LSMPI
          call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w1,la,0.0E0_realk,w2,nv)
          call dist_int_contributions(w2,int(no2*nv2,kind=long),gvoov_w)
#else
          call lsquit("ERROR(ccsd_integral_driven):this should never appear gvoov_w",-1)
#endif
        elseif(scheme==2)then
          call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w1,la,0.0E0_realk,w2,nv)
          call array_add(gvoova,1.0E0_realk,w2)
        endif
        call lsmpi_poke()
       endif
       call print_norm(gvoov,o2v2)

       !prepare w0 to contain L
       if(DECinfo%ccModel>2.and.(scheme==0.or.scheme==1)) call dscal(int(lg*la*nb2),2.0E0_realk,w0,1)

       IF(doscreen)Mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGabKLHS(alphaB)%p
       IF(doscreen)Mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGabKRHS(gammaB)%p

       !Mylsitem%setting%scheme%intprint=6
       call II_GET_DECPACKED4CENTER_K_ERI(DECinfo%output,DECinfo%output, &
            & Mylsitem%setting,w1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
            & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),dimAlpha,nb,dimGamma,nb,nbatches,INTSPEC,fullRHS)
       call lsmpi_poke()
       !Mylsitem%setting%scheme%intprint=0


      if(DECinfo%ccmodel>2)then
        if(scheme==0.or.scheme==1)then
          call get_d_term_int_direct(w0,w1,w2,w3,no,nv,nb,fa,fg,la,lg,&
          &xo,yo,xv,yv,u2,uigcj,omega2,u2kcjb,scheme)
          call lsmpi_poke()
        endif

        if(fa<=fg+lg-1)then
        !CHECK WHETHER THE TERM HAS TO BE DONE AT ALL, i.e. when the first
        !element in the alpha batch has a smaller index as the last element in
        !the gamma batch, chose the trafolength as minimum of alpha batch-length
        !and the difference between first element of alpha batch and last element
        !of gamma batch
        call get_a22_and_prepb22_terms_ex(w0,w1,w2,w3,tpl,tmi,no,nv,nb,fa,fg,la,lg,&
             &xo,yo,xv,yv,omega2,sio4,scheme)       
        call lsmpi_poke()

        endif
      endif
      call print_norm(omega2)
      

      !(w0):I[ delta gamma alpha beta] <- (w1):I[ alpha beta gamma delta ]
      call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,3,1,4],0.0E0_realk,w0)
      call lsmpi_poke()
      ! (w3):I[i gamma alpha beta] = Lambda^h[delta i] I[delta gamma alpha beta]
      call dgemm('t','n',no,lg*la*nb,nb,1.0E0_realk,yo,nb,w0,nb,0.0E0_realk,w2,no)
      call lsmpi_poke()
      ! (w0):I[i gamma alpha j] = (w3):I[i gamma alpha beta] Lambda^h[beta j]
      call dgemm('n','n',no*lg*la,no,nb,1.0E0_realk,w2,no*lg*la,yo,nb,0.0E0_realk,w0,no*lg*la)
      call lsmpi_poke()
      ! (w3):I[alpha gamma i j] <- (w0):I[i gamma alpha j]
      if(DECinfo%ccModel>2)call add_int_to_sio4(w0,w2,w3,no,nv,nb,fa,fg,la,lg,xo,sio4)
      call lsmpi_poke()
      call print_norm(sio4,int(no*no*nor,kind=8))


      ! (w2):I[gamma i j alpha] <- (w0):I[i gamma alpha j]
      call array_reorder_4d(1.0E0_realk,w0,no,lg,la,no,[2,1,4,3],0.0E0_realk,w2)
      call lsmpi_poke()
      ! (w3):I[b i j alpha] = Lamda^p[gamma b] (w2):I[gamma i j alpha]
      call dgemm('t','n',nv,no2*la,lg,1.0E0_realk,xv(fg),nb,w2,lg,0.0E0_realk,w3,nv)
      call lsmpi_poke()
      ! Omega += Lambda^p[alpha a]^T (w3):I[b i j alpha]^T
      if(scheme==1.or.scheme==2)then
        call dgemm('t','t',nv,o2v,la,0.5E0_realk,xv(fa),nb,w3,o2v,0.0E0_realk,w2,nv)
        call array_add(omega2,1.0E0_realk,w2)
      else
        call dgemm('t','t',nv,o2v,la,0.5E0_realk,xv(fa),nb,w3,o2v,1.0E0_realk,omega2%elm1,nv)
      endif
      call lsmpi_poke()
      call print_norm(omega2)

      if(DECinfo%ccmodel>2.and.(scheme==0.or.scheme==1))then
        call get_c_term_int_direct(w0,w1,w2,w3,no,nv,nb,fa,fg,la,lg,xo,&
             &yo,xv,yv,t2,tGammadij,omega2,t2jabi,scheme)
        call lsmpi_poke()
      endif

    end do BatchAlpha
    end do BatchGamma

    ! Free integral stuff
    ! *******************
    nullify(Mylsitem%setting%LST_GAB_LHS)
    nullify(Mylsitem%setting%LST_GAB_RHS)
    call free_decscreen(DECSCREEN)

    ! Free gamma stuff
    call mem_dealloc(orb2batchGamma)
    call mem_dealloc(batchdimGamma)
    call mem_dealloc(batchsizeGamma)
    call mem_dealloc(batchindexGamma)
    do idx=1,nbatchesGamma
       call mem_dealloc(batch2orbGamma(idx)%orbindex)
       batch2orbGamma(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbGamma)

    ! Free alpha stuff
    call mem_dealloc(orb2batchAlpha)
    call mem_dealloc(batchdimAlpha)
    call mem_dealloc(batchsizeAlpha)
    call mem_dealloc(batchindexAlpha)
    do idx=1,nbatchesAlpha
       call mem_dealloc(batch2orbAlpha(idx)%orbindex)
       batch2orbAlpha(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbAlpha)

    ! free arrays only needed in the batched loops
    call mem_dealloc(w0)
    call mem_dealloc(uigcj)
    call mem_dealloc(tpl)
    call mem_dealloc(tmi)

    if(scheme==1)then
      call array_free(t2jabi)
    endif
    !if(DECinfo%ccModel>2 .and. on )call mem_dealloc(tGammadij)
    if(DECinfo%ccModel>2.and.scheme==0) call mem_dealloc(tGammadij)

    ! free working matrices and adapt to new requirements
    call mem_dealloc(w1)
    call mem_dealloc(w2)
    call mem_dealloc(w3)
    
#ifdef VAR_LSMPI
    if(scheme==3)then
      call mem_alloc(gvvoo,gvvoo_c,int(no2*nv2,kind=long))
      call mem_alloc(gvoov,gvoov_c,int(no2*nv2,kind=long))
    endif
    ! Finish the MPI part of the Residual calculation
    startt=MPI_wtime()
    call lsmpi_barrier(infpar%lg_comm)
    stopp=MPI_wtime()
    wait_time = stopp - startt
    max_wait_time = wait_time


#ifdef VAR_DEBUG
    write(*,'("--rank",I2,", load: ",I5,", w-time:",f15.4)'),infpar%mynum,myload,wait_time
    call lsmpi_local_reduction(wait_time,infpar%master)
    call lsmpi_local_max(max_wait_time,infpar%master)
    if(master)then
      write(*,'("----------------------------------------------------------")')
      write(*,'("sum: ",f15.4," 0: ",f15.4," Max: ",f15.4)'),wait_time,wait_time/(infpar%nodtot*1.0E0_realk),max_wait_time
    endif
#endif
    if (master) call LSTIMER('CCSD part B',time_start,timewall_start,DECinfo%output)
    startt=MPI_wtime()
    if(infpar%lg_nodtot>1.or.scheme==3) then
       if(iter==1.and.scheme==4.or.scheme==0)then
         call lsmpi_local_allreduce_chunks(govov%elm1,o2v2,double_2G_nel)
         call lsmpi_barrier(infpar%lg_comm)
       elseif(scheme==3)then
         call array_cp_tiled2dense(govov,.false.)
       endif
       ! The following block is structured like this due to performance reasons
       !***********************************************************************
       if(DECinfo%ccModel>2)then
         call lsmpi_barrier(infpar%lg_comm)
         !print *,"calling allreduce",nor*no*no,double_2G_nel
         dummy64 = int(nor*no*no,kind=long)
         call lsmpi_local_allreduce_chunks(sio4,dummy64,double_2G_nel)
         !print *,"second allred done"
         call lsmpi_barrier(infpar%lg_comm)
         if(scheme==4)then
           call lsmpi_local_allreduce_chunks(gvvoo,o2v2,double_2G_nel)
           call lsmpi_local_allreduce_chunks(gvoov,o2v2,double_2G_nel)
         elseif(scheme==3)then
#ifdef VAR_MPIWIN
           dummy64 = int(no2*nv2,kind=long)
           call lsmpi_win_fence(gvoov_w,.true.)
           call collect_int_contributions_f(gvoov,dummy64,gvoov_w)
           call lsmpi_win_fence(gvoov_w,.false.)
           call lsmpi_win_fence(gvvoo_w,.true.)
           call collect_int_contributions_f(gvvoo,dummy64,gvvoo_w)
           call lsmpi_win_fence(gvvoo_w,.false.)
#else
           do nctr = 0,infpar%lg_nodtot-1
             call get_int_dist_info(o2v2,fe,ne,nctr)
             rcnt(nctr+1) = ne
             dsp(nctr+1)  = fe - 1
           enddo
           call lsmpi_local_allgatherv(gvvoo_r,gvvoo,rcnt,dsp)
           call lsmpi_local_allgatherv(gvoov_r,gvoov,rcnt,dsp)
#endif
         endif
       endif
    end if
    if(.not.dynamic_load)then
      call mem_dealloc(mpi_task_distribution)
    else
      call lsmpi_win_free(win_in_g)
      call mem_dealloc(mpi_stuff,mpi_ctasks)
    endif
    stopp=MPI_wtime()
#ifdef VAR_DEBUG
    if(master) print*,"MPI part of the calculation finished, comm-time",stopp-startt
#endif    
    !free windows and deallocate partial int matrices in scheme 1
    if(DECinfo%ccModel>2.and.(scheme==3.or.scheme==4))then
      call lsmpi_win_free(gvoov_w)
      call lsmpi_win_free(gvvoo_w)
      if(scheme==3)then
        call mem_dealloc(gvvoo_r,gvvoo_p)
        call mem_dealloc(gvoov_r,gvoov_p)
      endif
    endif
#endif


    ! Reallocate 1 temporary array
    maxsize64 = max(int(nv2*no2,kind=8),int(nb2,kind=8))
    maxsize64 = max(maxsize64,int(nv2*nor,kind=8))
    call mem_alloc(w1,maxsize64)
    call dcopy(int(o2v2),omega2%elm1,1,w1,1)
    call lsmpi_local_reduction(w1,o2v2,infpar%master,double_2G_nel)
    if(master)call print_norm(w1,o2v2)
    w1=0.0E0_realk

    !reorder integral for use within the solver and the c and d terms
    if(iter==1.and.(scheme==4.or.scheme==0))then
      call array_reorder_4d(1.0E0_realk,govov%elm1,no,no,nv,nv,[1,4,2,3],0.0E0_realk,w1)
      call dcopy(no2*nv2,w1,1,govov%elm1,1)
#ifdef VAR_LSMPI
      if(DECinfo%solver_par)then
        govov%atype     = TILED_DIST
        call array_convert(govov%elm1,govov)
        if(scheme==4)govov%atype = DENSE
      endif
#endif
    endif

    if(DECinfo%ccModel>2)then

      !get B2.2 contributions
      !**********************
      call get_B22_contrib_mo(sio4,t2,w1,w2,no,nv,nb,omega2,scheme)
#ifdef VAR_LSMPI
      call lsmpi_win_free(sio4_w)
      call mem_dealloc(sio4,sio4_c)
#else
      call mem_dealloc(sio4)
#endif
      if(scheme==4.or.scheme==3.or.scheme==2)then

        !Get the C2 and D2 terms if enough memory is available
        !*****************************************************
#ifdef VAR_OMP
        startt=omp_get_wtime()
#elif VAR_LSMPI
        startt=MPI_wtime()
#endif
        call get_cnd_terms_mo(w1,w2,w3,t2,u2,govov,gvoov,gvvoo,no,nv,omega2,gvvooa,gvoova,scheme)
        if(scheme==4.or.scheme==3)then
#ifdef VAR_LSMPI
          call mem_dealloc(gvoov,gvoov_c)
          call mem_dealloc(gvvoo,gvvoo_c)
#else
          call mem_dealloc(gvoov)
          call mem_dealloc(gvvoo)
#endif
        elseif(scheme==1.or.scheme==2)then
          call array_free(gvoova)
          call array_free(gvvooa)
        endif
#ifdef VAR_OMP
        stopp=omp_get_wtime()
#elif VAR_LSMPI
        stopp=MPI_wtime()
#endif

!OUTPUT
#ifdef VAR_LSMPI
        write(*,'(I3,"C and D   :",f15.4)'),infpar%lg_mynum,stopp-startt
#else
        write(DECinfo%output,'("C and D   :",f15.4)')stopp-startt
#endif
      endif
    endif


    !IN CASE OF MPI (AND CORRECT SCHEEME) REDUCE TO MASTER
    !*****************************************************
#ifdef VAR_LSMPI
    if(infpar%lg_nodtot>1) then
      if(scheme==4.or.scheme==3.or.scheme==0)&
       &call lsmpi_local_reduction(omega2%elm1,o2v2,infpar%master,double_2G_nel)
      call lsmpi_local_reduction(Gbi,nb*no,infpar%master,double_2G_nel)
      call lsmpi_local_reduction(Had,nb*nv,infpar%master,double_2G_nel)
    endif
    !convert stuff
    !set for correct access again, save as i a j b
    if(DECinfo%solver_par)then
      if((master.and..not.(scheme==2.or.scheme==1)).or.scheme==3)&
      &call memory_deallocate_array_dense(govov)
      govov%atype      = TILED_DIST
    endif
    govov%init_type  = MASTER_INIT
    omega2%init_type = MASTER_INIT
    t2%init_type     = MASTER_INIT
    if(scheme==1.or.scheme==2) u2%init_type     = MASTER_INIT
#endif
    

 
    ! slaves should exit the subroutine after the main work is done
    if(.not. master) then
      call mem_dealloc(w1)
      call mem_dealloc(Had)
      call mem_dealloc(Gbi)
      if(scheme==4.or.scheme==3.or.scheme==0)call array_free(u2)
      return
    endif
    call print_norm(omega2%elm1,o2v2)


    !allocate the density matrix
    call mat_init(iFock,nb,nb)
    call mat_init(Dens,nb,nb)

    !calculate inactive fock matrix in ao basis
    call dgemm('n','t',nb,nb,no,1.0E0_realk,yo,nb,xo,nb,0.0E0_realk,Dens%elms,nb)
    call mat_zero(iFock)
    call dec_fock_transformation(iFock,Dens,MyLsItem,.false.)
    

    !THIS IS NOT YET IMPLEMENTED -- as soon as it is, do not use type(matrix)
    !anymore
    !call II_get_fock_mat_full(DECinfo%output,DECinfo%output,MyLsItem%setting,nb,&
    !& Dens%elms,.false.,iFock%elms)
    !use dens as temporay array 


    call ii_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsItem%setting,&
         & Dens%elms,nb,nb,AORdefault,AORdefault)
    ! Add one- and two-electron contributions to Fock matrix
    call daxpy(nb2,1.0E0_realk,Dens%elms,1,iFock%elms,1)
    !Free the density matrix
    call mat_free(Dens)

    ! KK: Add long-range Fock correction
    call daxpy(nb2,1.0E0_realk,deltafock,1,iFock%elms,1)
#ifdef VAR_OMP
    startt=omp_get_wtime()
#elif VAR_LSMPI
    startt=MPI_wtime()
#endif
    !Transform inactive Fock matrix into the different mo subspaces
    if (DECinfo%ccModel>2) then
      ! -> Foo
      call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,iFock%elms,nb,0.0E0_realk,w1,no)
      call dgemm('n','n',no,no,nb,1.0E0_realk,w1,no,yo,nb,0.0E0_realk,ppfock,no)
      ! -> Fov
      call dgemm('n','n',no,nv,nb,1.0E0_realk,w1,no,yv,nb,0.0E0_realk,pqfock,no)
      ! -> Fvo
      call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,iFock%elms,nb,0.0E0_realk,w1,nv)
      call dgemm('n','n',nv,no,nb,1.0E0_realk,w1,nv,yo,nb,0.0E0_realk,qpfock,nv)
      ! -> Fvv
      call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1,nv,yv,nb,0.0E0_realk,qqfock,nv)
    else
      ! -> Foo
      call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,fock,nb,0.0E0_realk,w1,no)
      call dgemm('n','n',no,no,nb,1.0E0_realk,w1,no,yo,nb,0.0E0_realk,ppfock,no)
      ! -> Fov
      call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,iFock%elms,nb,0.0E0_realk,w1,no)
      call dgemm('n','n',no,nv,nb,1.0E0_realk,w1,no,yv,nb,0.0E0_realk,pqfock,no)
      ! -> Fvo
      call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,iFock%elms,nb,0.0E0_realk,w1,nv)
      call dgemm('n','n',nv,no,nb,1.0E0_realk,w1,nv,yo,nb,0.0E0_realk,qpfock,nv)
      ! -> Fvv
      call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,fock,nb,0.0E0_realk,w1,nv)
      call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1,nv,yv,nb,0.0E0_realk,qqfock,nv)
    endif

    !Free the AO fock matrix
    call mat_free(iFock)


#ifdef VAR_OMP
    stopp=omp_get_wtime()
#elif VAR_LSMPI
    stopp=MPI_wtime()
#endif
    write(DECinfo%output,'("Fock trafo:",f15.4)')stopp-startt
#ifdef VAR_OMP
    startt=omp_get_wtime()
#elif VAR_LSMPI
    startt=MPI_wtime()
#endif


    !GET SINGLES CONTRIBUTIONS
    !*************************

    !calculate singles J term
    ! F [a i] = Omega [a i]
    call dcopy(no*nv,qpfock,1,omega1,1)

    !calculate singles I term
    ! Reorder u [c a i k] -> u [a i c k]
    if(scheme==4.or.scheme==3.or.scheme==0)then
      call array_reorder_4d(1.0E0_realk,u2%elm1,nv,nv,no,no,[2,3,4,1],0.0E0_realk,w1)
    elseif(scheme==2.or.scheme==1)then
      call array_convert(u2,w1,[4,1,2,3])
    endif
    ! u [a i k c] * F[k c] =+ Omega [a i]
    call dgemv('n',nv*no,nv*no,1.0E0_realk,w1,nv*no,pqfock,1,1.0E0_realk,omega1,1)

    !calculate singles G term
    ! Lambda^p [alpha a]^T Gbi [alpha i] =+ Omega [a i]
    call dgemm('t','n',nv,no,nb,1.0E0_realk,xv,nb,Gbi,nb,1.0E0_realk,omega1,nv)
    !calculate singles H term
    ! (-1) Had [a delta] * Lambda^h [delta i] =+ Omega[a i]
    call dgemm('n','n',nv,no,nb,-1.0E0_realk,Had,nv,yo,nb,1.0E0_realk,omega1,nv)
    

    !GET DOUBLES E2 TERM - AND INTRODUCE PERMUTATIONAL SYMMMETRY
    !***********************************************************
    call calculate_E2_and_permute(ppfock,qqfock,w1,t2,xo,yv,Gbi,Had,no,nv,nb,omega2,scheme)

    call mem_dealloc(Had)
    call mem_dealloc(Gbi)

#ifdef VAR_OMP
    stopp=omp_get_wtime()
#elif VAR_LSMPI
    stopp=MPI_wtime()
#endif
    write(DECinfo%output,'("S and E   :",f15.4)')stopp-startt


#ifdef VAR_LSMPI
    if(DECinfo%solver_par.and.(scheme==4.or.scheme==3.or.scheme==0))then
      call array_mv_dense2tiled(omega2,.true.)
      call array_mv_dense2tiled(t2,.true.)
    endif
#endif
    call mem_dealloc(w1)
    call array_free(u2)

    call LSTIMER('CCSD part C',time_start,timewall_start,DECinfo%output)
    call LSTIMER('CCSD RESIDUAL',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu_end,twall_end,DECinfo%output)

  end subroutine get_ccsd_residual_integral_driven

  subroutine calculate_E2_and_permute(ppf,qqf,w1,t2,xo,yv,Gbi,Had,no,nv,nb,omega2,s)
    implicit none
    real(realk),intent(inout)::ppf(:)
    real(realk),intent(inout)::qqf(:)
    real(realk),pointer:: w1(:)
    type(array),intent(inout) :: t2
    real(realk),pointer:: xo(:)
    real(realk),pointer:: yv(:)
    real(realk),pointer:: Gbi(:)
    real(realk),pointer:: Had(:)
    integer, intent(in) :: no,nv,nb
    type(array),intent(inout) :: omega2
    integer, intent(in) :: s
    integer :: no2,nv2,v2o,o2v
    logical :: master 
    real(realk),pointer :: w2(:),w3(:)
    integer :: i,ml
    integer(kind=ls_mpik) :: me,nnod,nod
    integer :: ml1,fai1,l1,tl1,lai1
    integer :: ml2,fai2,l2,tl2,lai2
    integer :: fri,tri
    character(ARR_MSG_LEN) :: msg
    real(realk) :: nrm
    master = .true.
    nrm=0.0E0_realk

    no2=no*no
    nv2=nv*nv
    v2o=nv*nv*no
    o2v=no*no*nv
    
#ifdef VAR_LSMPI
    master=(infpar%lg_mynum==infpar%master)
    if((s==2.or.s==1).and.master)then
      call share_E2_with_slaves(ppf,qqf,t2,xo,yv,Gbi,Had,no,nv,nb,omega2,s)
    endif
#endif

    if(s==4.or.s==3.or.s==0)then
      !calculate first part of doubles E term and its permutation
      ! F [k j] + Lambda^p [alpha k]^T * Gbi [alpha j] = G' [k j]
      call dcopy(no2,ppf,1,w1,1)
      if (DECinfo%ccModel>2) call dgemm('t','n',no,no,nb,1.0E0_realk,xo,nb,Gbi,nb,1.0E0_realk,w1,no)
      ! (-1) t [a b i k] * G' [k j] =+ Omega [a b i j]
      call dgemm('n','n',v2o,no,no,-1.0E0_realk,t2%elm1,v2o,w1,no,1.0E0_realk,omega2%elm1,v2o)
     
      !calculate second part of doubles E term
      ! F [b c] - Had [a delta] * Lambda^h [delta c] = H' [b c]
      call dcopy(nv2,qqf,1,w1,1)
      if (DECinfo%ccModel>2) call dgemm('n','n',nv,nv,nb,-1.0E0_realk,Had,nv,yv,nb,1.0E0_realk,w1,nv)
      ! H'[a c] * t [c b i j] =+ Omega [a b i j]
      call dgemm('n','n',nv,o2v,nv,1.0E0_realk,w1,nv,t2%elm1,nv,1.0E0_realk,omega2%elm1,nv)
     
      !INTRODUCE PERMUTATION
      call dcopy(nv2*no2,omega2%elm1,1,w1,1)
      call array_reorder_4d(1.0E0_realk,w1,nv,nv,no,no,[2,1,4,3],1.0E0_realk,omega2%elm1)
    elseif(s==1.or.s==2)then
#ifdef VAR_LSMPI
       omega2%init_type=ALL_INIT
       t2%init_type=ALL_INIT
       me   = 0
       nnod = 1
       nnod = infpar%lg_nodtot
       me   = infpar%lg_mynum
      
      !Setting transformation variables for each rank
      !**********************************************
      call mo_work_dist(nv*nv*no,fai1,tl1)
      call mo_work_dist(nv*no*no,fai2,tl2)

      call mem_alloc(w3,max(tl1*no,tl2*nv))
      call mem_alloc(w2,max(nv2,no2))

      !DO ALL THINGS DEPENDING ON 1

      !calculate first part of doubles E term and its permutation
      ! F [k j] + Lambda^p [alpha k]^T * Gbi [alpha j] = G' [k j]
      call dcopy(no2,ppf,1,w2,1)
      if (DECinfo%ccModel>2) call dgemm('t','n',no,no,nb,1.0E0_realk,xo,nb,Gbi,nb,1.0E0_realk,w2,no)
      ! (-1) t [a b i k] * G' [k j] =+ Omega [a b i j]
      !if(me==0) call array_convert(t2,w1,t2%nelms)
      call array_gather_tilesinfort(t2,w1,int(nv*nv*no*no,kind=long),infpar%master)
      do nod=1,nnod-1
        call mo_work_dist(nv*nv*no,fri,tri,nod)
        if(me==0)then
          do i=1,no
            call dcopy(tri,w1(fri+(i-1)*no*nv*nv),1,w3(1+(i-1)*tri),1)
          enddo
        endif
        if(me==0.or.me==nod)then
          call ls_mpisendrecv(w3,int(no*tri,kind=long),infpar%lg_comm,infpar%master,nod)
        endif
      enddo
      if(me==0)then
        do i=1,no
          call dcopy(tl1,w1(fai1+(i-1)*no*nv*nv),1,w3(1+(i-1)*tl1),1)
        enddo
      endif
      w1=0.0E0_realk
      call dgemm('n','n',tl1,no,no,-1.0E0_realk,w3,tl1,w2,no,0.0E0_realk,w1(fai1),v2o)
      call lsmpi_local_reduction(w1,int(nv*nv*no*no,kind=long),infpar%master)
      call array_scatteradd_densetotiled(omega2,1.0E0_realk,w1,int(nv*nv*no*no,kind=long),infpar%master)
      !if(me==0)then
      !  call array_add(omega2,1.0E0_realk,w1,omega2%nelms)
      !endif


      !DO ALL THINGS DEPENDING ON 2

      !calculate second part of doubles E term
      ! F [b c] - Had [a delta] * Lambda^h [delta c] = H' [b c]
      call dcopy(nv2,qqf,1,w2,1)
      if (DECinfo%ccModel>2) call dgemm('n','n',nv,nv,nb,-1.0E0_realk,Had,nv,yv,nb,1.0E0_realk,w2,nv)
      ! H'[a c] * t [c b i j] =+ Omega [a b i j]
      !if(me==0) call array_convert(t2,w1,t2%nelms)
      call array_gather_tilesinfort(t2,w1,int(nv*nv*no*no,kind=long),infpar%master)
      do nod=1,nnod-1
        call mo_work_dist(nv*no*no,fri,tri,nod)
        if(me==0)then
          do i=1,tri
            call dcopy(nv,w1(1+(fri+i-2)*nv),1,w3(1+(i-1)*nv),1)
          enddo
        endif
        if(me==0.or.me==nod)then
          call ls_mpisendrecv(w3,int(nv*tri,kind=long),infpar%lg_comm,infpar%master,nod)
        endif
      enddo
      if(me==0)then
        do i=1,tl2
          call dcopy(nv,w1(1+(fai2+i-2)*nv),1,w3(1+(i-1)*nv),1)
        enddo
      endif
      w1=0.0E0_realk
      call dgemm('n','n',nv,tl2,nv,1.0E0_realk,w2,nv,w3,nv,0.0E0_realk,w1(1+(fai2-1)*nv),nv)
      call lsmpi_local_reduction(w1,int(nv*nv*no*no,kind=long),infpar%master)
      call array_scatteradd_densetotiled(omega2,1.0E0_realk,w1,int(nv*nv*no*no,kind=long),infpar%master)
      !if(me==0)then
      !  call array_add(omega2,1.0E0_realk,w1,omega2%nelms)
      !endif
      
      call mem_dealloc(w2)
      call mem_dealloc(w3)
      
      !INTRODUCE PERMUTATION
      !please note: the barrier statements are crucial because one sided
      !communcation is used
      omega2%init_type=MASTER_INIT
      t2%init_type=MASTER_INIT
      !call lsmpi_barrier(infpar%lg_comm)
      !if(me==0)then
      !  call array_convert(omega2,w1,omega2%nelms)
      !  call array_add(w1,1.0E0_realk,omega2,omega2%nelms,[2,1,4,3])
      !call lsmpi_barrier(infpar%lg_comm)
      !  call array_convert(w1,omega2,omega2%nelms)
      !endif
      !call lsmpi_barrier(infpar%lg_comm)
      call array_gather_tilesinfort(omega2,w1,int(nv*nv*no*no,kind=long),infpar%master)
      call array_gatheradd_tilestofort(omega2,1.0E0_realk,w1,int(nv*nv*no*no,kind=long),&
                          &infpar%master,[2,1,4,3])
      call array_scatter_densetotiled(omega2,w1,int(no*no*nv*nv,kind=long),infpar%master)
#endif
    endif
    
  end subroutine calculate_E2_and_permute


  subroutine check_job(s,fr,dyn,a,g,na,ng,static,dynamic,prnt)
    implicit none
    logical,intent(in) :: dyn,prnt
    logical,intent(inout) :: fr
    integer,intent(in) :: s,g,na,ng
    integer,intent(inout) :: a
    integer :: static(:)
    !integer(kind=ls_mpik) :: dynamic(:)
    integer(kind=ls_mpik) :: dynamic
    real(realk) :: mpi_buf,el 
    integer(kind=ls_mpik) :: i, job
#ifdef VAR_LSMPI
       el=0
       !ugly construction to get both schemes in
       if(.not.dyn)then
         a=a+1
         do while(a<=na)
           if(static((a-1)*ng+g)/=infpar%lg_mynum)then
             a=a+1
           else
             a=a-1
             exit
           endif
         enddo
       else
         !BE CAREFUL THIS IS HIGHLY EXPERIMENTAL CODE
         if(fr)then
           a=infpar%lg_mynum
         else
           mpi_buf=1.0E0_realk
           call lsmpi_win_lock(infpar%master,dynamic,'e')
           call lsmpi_get(el,g,infpar%master,dynamic)
           call lsmpi_acc(mpi_buf,g,infpar%master,dynamic)
           call lsmpi_win_unlock(infpar%master,dynamic)
           if(s==4.or.s==1)then
             do i=1,infpar%lg_nodtot-1
               call lsmpi_win_lock(i,dynamic,'s')
               call lsmpi_win_unlock(i,dynamic)
             enddo
           endif
           a=int(el)
         endif
       endif
       if(fr) fr=.false.
#endif
       a=a+1
       if(a>na)return
#ifdef VAR_LSMPI
       if(prnt)write (*, '("Rank ",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")'),infpar%mynum,&
       &a,na,g,ng
#else
       write (*, '("starting job (",I3,"/",I3,",",I3,"/",I3,")")')a,&
       &na,g,ng
#endif
       call lsmpi_poke()
  end subroutine check_job
  
  subroutine mo_work_dist(m,fai,tl,nod)
    implicit none
    integer,intent(in) :: m
    integer,intent(inout)::fai
    integer,intent(inout)::tl
    integer(kind=ls_mpik),optional,intent(inout)::nod
    integer :: l,ml,me,nnod
    
    me   = 0
    nnod = 1
#ifdef VAR_LSMPI
    nnod = infpar%lg_nodtot
    me   = infpar%lg_mynum
#endif
      
    if(present(nod))me=nod

    !Setting transformation variables for each rank
    !**********************************************
    l   = (m) / nnod
    ml  = mod(m,nnod)
    fai = me * l + 1
    tl  = l

    if(ml>0)then
      if(me<ml)then
        fai = fai + me
        tl  = l + 1
      else
        fai = fai + ml
        tl  = l
      endif
    endif

  end subroutine mo_work_dist

  !> \brief Routine to get the c and the d terms from t1 tranformed integrals
  !using a simple mpi-parallelization
  !> \author Patrick Ettenhuber
  !> \Date January 2013 
  subroutine get_cnd_terms_mo(w1,w2,w3,t2,u2,govov,gvoov,gvvoo,no,nv,omega2,gvvooa,gvoova,s)
    implicit none
    !> input some empty workspace of zise v^2*o^2 
    real(realk), intent(inout) :: w1(:)
    real(realk),pointer :: w2(:),w3(:)
    !> the t1-transformed integrals
    real(realk), intent(in) :: gvvoo(:),gvoov(:)
    type(array), intent(inout) :: govov,gvvooa,gvoova
    !> number of occupied orbitals 
    integer, intent(in) :: no
    !> nuber of virtual orbitals
    integer, intent(in) :: nv
    !> ampitudes on input ordered as abij
    !real(realk), intent(in) :: t2(:)
    type(array), intent(inout) :: t2
    !> u on input u{aibj}=2t{aibj}-t{ajbi} ordered as abij
    !real(realk), intent(in) :: u2(:)
    type(array), intent(inout) :: u2
    !> the residual to add the contribution
    !real(realk), intent(inout) :: omega2(:)
    type(array), intent(inout) :: omega2
    !> integer specifying the scheme
    integer, intent(in) :: s
    integer :: tl,fai,lai,i,faif,lead
    integer :: l,ml
    integer(kind=ls_mpik) :: nod,me,nnod
    real(realk) :: nrm1,nrm2,nrm3,nrm4
    integer :: a,b,j,fri,tri
    character(ARR_MSG_LEN) :: msg



      me   = 0
      nnod = 1
#ifdef VAR_LSMPI
      nnod = infpar%lg_nodtot
      me   = infpar%lg_mynum
#endif
      
      !Setting transformation variables for each rank
      !**********************************************
      call mo_work_dist(nv*no,fai,tl)
     
      if(s==4)then
        faif = fai
        lead = no * nv
        call mem_alloc(w2,int(no*no*nv*nv,kind=long))
        call mem_alloc(w3,int(no*no*nv*nv,kind=long))
      elseif(s==3.or.s==2)then
        faif = 1
        lead = tl
        call mem_alloc(w2,int(tl*no*nv,kind=long))
        call mem_alloc(w3,int(tl*no*nv,kind=long))
      else
        call lsquit("ERROR(get_cnd_terms_mo):no valid scheme",-1)
      endif


      !calculate doubles C term
      !*************************

      if(s==2.or.s==3.or.s==4)then

        !Reorder gvvoo [a c k i] -> goovv [a i c k]
        if(s==4)then
          call array_reorder_4d(1.0E0_realk,gvvoo,nv,no,no,nv,[1,3,4,2],0.0E0_realk,w2)
        elseif(s==3)then
          call array_reorder_4d(1.0E0_realk,gvvoo,nv,no,no,nv,[1,3,4,2],0.0E0_realk,w1)
          do i=1,tl
            call dcopy(no*nv,w1(fai+i-1),no*nv,w2(i),tl)
          enddo
        elseif(s==2)then
#ifdef VAR_LSMPI
          !for a conversion from tiled to dense the reverse order has to be
          !given
          !if(me==0) call array_convert(gvvooa,w1,gvvooa%nelms,[1,4,2,3])
          call array_gather_tilesinfort(gvvooa,w1,int(no*no*nv*nv,kind=long),infpar%master,[1,4,2,3])
          do nod=1,nnod-1
            call mo_work_dist(no*nv,fri,tri,nod)
            if(me==0)then
              do i=1,tri
                call dcopy(no*nv,w1(fri+i-1),no*nv,w2(i),tri)
              enddo
            endif
            if(me==0.or.me==nod)then
              call ls_mpisendrecv(w2,int(no*nv*tri,kind=long),infpar%lg_comm,infpar%master,nod)
            endif
          enddo
          if(me==0)then
            do i=1,tl
              call dcopy(no*nv,w1(fai+i-1),no*nv,w2(i),tl)
            enddo
          endif
#endif
        endif


        !Reorder t [a d l i] -> t [a i d l]
        if(s==4)then
          call array_reorder_4d(1.0E0_realk,t2%elm1,nv,nv,no,no,[1,4,2,3],0.0E0_realk,w3)
        elseif(s==3)then
          call array_reorder_4d(1.0E0_realk,t2%elm1,nv,nv,no,no,[1,4,2,3],0.0E0_realk,w1)
          do i=1,tl
            call dcopy(no*nv,w1(fai+i-1),no*nv,w3(i),tl)
          enddo
        elseif(s==2)then
#ifdef VAR_LSMPI
          !if(me==0) call array_convert(t2,w1,t2%nelms,[1,3,4,2])
          call array_gather_tilesinfort(t2,w1,int(no*no*nv*nv,kind=long),infpar%master,[1,3,4,2])
          do nod=1,nnod-1
            call mo_work_dist(no*nv,fri,tri,nod)
            if(me==0)then
              do i=1,tri
                call dcopy(no*nv,w1(fri+i-1),no*nv,w3(i),tri)
              enddo
            endif
            if(me==0.or.me==nod)then
              call ls_mpisendrecv(w3,int(no*nv*tri,kind=long),infpar%lg_comm,infpar%master,nod)
            endif
          enddo
          if(me==0)then
            do i=1,tl
              call dcopy(no*nv,w1(fai+i-1),no*nv,w3(i),tl)
            enddo
          endif
#endif
        endif


        !Reorder govov [k d l c] -> govov [d l c k]
        if(s==3.or.s==4)then
          call array_reorder_4d(1.0E0_realk,govov%elm1,no,nv,no,nv,[2,3,4,1],0.0E0_realk,w1)
        elseif(s==2)then
#ifdef VAR_LSMPI
          !if(me==0)call array_convert(govov,w1,govov%nelms,[4,1,2,3])
          call array_gather_tilesinfort(govov,w1,int(no*no*nv*nv,kind=long),infpar%master,[4,1,2,3])
          call ls_mpibcast(w1,int(no*no*nv*nv,kind=long),infpar%master,infpar%lg_comm)
#endif
        endif
       
        !(-0.5) * t [a i d l] * govov [d l c k] + goovv [a i c k] = C [a i c k]
        call dgemm('n','n',tl,no*nv,no*nv,-0.5E0_realk,w3(faif),lead,w1,no*nv,1.0E0_realk,w2(faif),lead)
        !(-1) * C [a i c k] * t [c k b j] = preOmC [a i b j]
        if(s==4)then
          w1=0.0E0_realk
          call dgemm('n','t',tl,no*nv,no*nv,-1.0E0_realk,w2(faif),lead,w3,no*nv,0.0E0_realk,w1(fai),no*nv)
        elseif(s==3)then
          call array_reorder_4d(1.0E0_realk,t2%elm1,nv,nv,no,no,[1,4,2,3],0.0E0_realk,w1)
          call dgemm('n','t',tl,no*nv,no*nv,-1.0E0_realk,w2(faif),lead,w1,no*nv,0.0E0_realk,w3,lead)
          w1=0.0E0_realk
          do i=1,tl
            call dcopy(no*nv,w3(i),tl,w1(fai+i-1),no*nv)
          enddo
        elseif(s==2)then
#ifdef VAR_LSMPI
          !call array_convert(t2,w1,t2%nelms,[1,3,4,2])
          call array_gather_tilesinfort(t2,w1,int(no*no*nv*nv,kind=long),infpar%master,[1,3,4,2])
          call ls_mpibcast(w1,int(no*no*nv*nv,kind=long),infpar%master,infpar%lg_comm)
          call dgemm('n','t',tl,no*nv,no*nv,-1.0E0_realk,w2(faif),lead,w1,no*nv,0.0E0_realk,w3,lead)
          w1=0.0E0_realk
          do i=1,tl
            call dcopy(no*nv,w3(i),tl,w1(fai+i-1),no*nv)
          enddo
#endif
        endif
        if(s==3.or.s==4)then
          !contribution 1: 0.5*preOmC [a i b j] -> =+ Omega [a b i j]
          call array_reorder_4d(0.5E0_realk,w1,nv,no,nv,no,[1,3,2,4],1.0E0_realk,omega2%elm1)
          !contribution 3: preOmC [a j b i] -> =+ Omega [a b i j]
          call array_reorder_4d(1.0E0_realk,w1,nv,no,nv,no,[1,3,4,2],1.0E0_realk,omega2%elm1)
        elseif(s==2)then
#ifdef VAR_LSMPI
          call lsmpi_local_reduction(w1,nv*nv*no*no,infpar%master)
          call array_scatteradd_densetotiled(omega2,0.5E0_realk,w1,int(no*no*nv*nv,kind=long),infpar%master,[1,3,2,4])
          call array_scatteradd_densetotiled(omega2,1.0E0_realk,w1,int(no*no*nv*nv,kind=long),infpar%master,[1,3,4,2])
          !if(me==0)then
          !  call array_add(omega2,0.5E0_realk,w1,omega2%nelms,[1,3,2,4])
          !  call array_add(omega2,1.0E0_realk,w1,omega2%nelms,[1,3,4,2])
          !endif
#endif
        endif
      else
        call lsquit("ERROR(get_cnd_terms_mo):no valid scheme",-1)
      endif


      !calculate doubles D term
      !************************
      if(s==4.or.s==3.or.s==2)then
        !(-1) * gvvoo [a c k i] -> + 2*gvoov[a i k c] = L [a i k c]
        !write (msg,*),"gvoov"

        if(s==4)then
          !call print_norm(gvoov,no*no*nv*nv,msg)
          call array_reorder_4d(2.0E0_realk,gvoov,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w2)
          call array_reorder_4d(-1.0E0_realk,gvvoo,nv,no,no,nv,[1,3,2,4],1.0E0_realk,w2)
          !call print_norm(w2,no*no*nv*nv)
        elseif(s==3)then
          call array_reorder_4d(2.0E0_realk,gvoov,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w1)
          call array_reorder_4d(-1.0E0_realk,gvvoo,nv,no,no,nv,[1,3,2,4],1.0E0_realk,w1)
          do i=1,tl
            call dcopy(no*nv,w1(fai+i-1),no*nv,w2(i),tl)
          enddo
          !write (msg,*),infpar%lg_mynum,"has L1:"
          !call print_norm(w2,lead*no*nv,msg)
        elseif(s==2)then
#ifdef VAR_LSMPI
          !call print_norm(gvoova,msg)
          call array_gather_tilesinfort(gvoova,w1,int(no*no*nv*nv,kind=long),infpar%master,[1,3,4,2])
          if(me==0)then
            !call array_convert(gvoova,w1,gvoova%nelms,[1,3,4,2])
            call dscal(int(gvoova%nelms),2.0E0_realk,w1,1)
            !call array_add(w1,-1.0E0_realk,gvvooa,gvvooa%nelms,[1,3,2,4])
          endif
          call array_gatheradd_tilestofort(gvvooa,-1.0E0_realk,w1,int(no*no*nv*nv,kind=long),infpar%master,[1,3,2,4])
          do nod=1,nnod-1
            call mo_work_dist(no*nv,fri,tri,nod)
            if(me==0)then
              do i=1,tri
                call dcopy(no*nv,w1(fri+i-1),no*nv,w2(i),tri)
              enddo
            endif
            if(me==0.or.me==nod)then
              call ls_mpisendrecv(w2,int(no*nv*tri,kind=long),infpar%lg_comm,infpar%master,nod)
            endif
          enddo
          if(me==0)then
            do i=1,tl
              call dcopy(no*nv,w1(fai+i-1),no*nv,w2(i),tl)
            enddo
          endif
#endif
        endif

        !Transpose u [d a i l] -> u [a i l d]
        if(s==4)then
          !call mat_transpose(u2%elm1,nv,no*no*nv,w3)
          call array_reorder_4d(1.0E0_realk,u2%elm1,nv,nv,no,no,[2,3,4,1],0.0E0_realk,w3)
        elseif(s==3)then
          !call mat_transpose(u2%elm1,nv,no*no*nv,w1)
          call array_reorder_4d(1.0E0_realk,u2%elm1,nv,nv,no,no,[2,3,4,1],0.0E0_realk,w1)
          do i=1,tl
            call dcopy(no*nv,w1(fai+i-1),no*nv,w3(i),tl)
          enddo
        elseif(s==2)then
#ifdef VAR_LSMPI
          !if(me==0)call array_convert(u2,w1,u2%nelms,[4,1,2,3])
          call array_gather_tilesinfort(u2,w1,int(no*no*nv*nv,kind=long),infpar%master,[4,1,2,3])
          do nod=1,nnod-1
            call mo_work_dist(no*nv,fri,tri,nod)
            if(me==0)then
              do i=1,tri
                call dcopy(no*nv,w1(fri+i-1),no*nv,w3(i),tri)
              enddo
            endif
            if(me==0.or.me==nod)then
              call ls_mpisendrecv(w3,int(no*nv*tri,kind=long),infpar%lg_comm,infpar%master,nod)
            endif
          enddo
          if(me==0)then
            do i=1,tl
              call dcopy(no*nv,w1(fai+i-1),no*nv,w3(i),tl)
            enddo
          endif
#endif
        endif
        !call print_norm(w3,no*no*nv*nv)

        !(-1) * govov [l c k d] + 2*govov[l d k c] = L [l d k c]
        if(s==3.or.s==4)then
          call array_reorder_4d(2.0E0_realk,govov%elm1,no,nv,no,nv,[1,2,3,4],0.0E0_realk,w1)
          call array_reorder_4d(-1.0E0_realk,govov%elm1,no,nv,no,nv,[1,4,3,2],1.0E0_realk,w1)
        elseif(s==2)then
#ifdef VAR_LSMPI
          call array_gather_tilesinfort(govov,w1,int(nv*nv*no*no,kind=long),infpar%master)
          if(me==0)then
            !call array_convert(govov,w1,govov%nelms)
            call dscal(int(govov%nelms),2.0E0_realk,w1,1)
            !call array_add(w1,-1.0E0_realk,govov,govov%nelms,[1,4,3,2])
          endif
          call array_gatheradd_tilestofort(govov,-1.0E0_realk,w1,int(no*no*nv*nv,kind=long),infpar%master,[1,4,3,2])
          call ls_mpibcast(w1,int(nv*nv*no*no,kind=long),infpar%master,infpar%lg_comm)
#endif
        endif

       
        ! (0.5) * u [a i l d] * L [l d k c] + L [a i k c] = D [a i k c]
        call dgemm('n','n',tl,nv*no,nv*no,0.5E0_realk,w3(faif),lead,w1,nv*no,1.0E0_realk,w2(faif),lead)
        ! (0.5)*D[a i k c] * u [b j k c]^T  = preOmD [a i b j]
        if(s==4)then
          w1=0.0E0_realk
          call dgemm('n','t',tl,nv*no,nv*no,0.5E0_realk,w2(faif),lead,w3,nv*no,0.0E0_realk,w1(fai),nv*no)
          !call print_norm(w1,lead*no*nv)
        elseif(s==3)then
          !call mat_transpose(u2%elm1,nv,no*no*nv,w1)
          call array_reorder_4d(1.0E0_realk,u2%elm1,nv,nv,no,no,[2,3,4,1],0.0E0_realk,w1)
          call dgemm('n','t',tl,nv*no,nv*no,0.5E0_realk,w2(faif),lead,w1,nv*no,0.0E0_realk,w3,lead)
          !write (msg,*),infpar%lg_mynum,"has 2:"
          !call print_norm(w3,lead*no*nv,msg)
          w1=0.0E0_realk
          do i=1,tl
            call dcopy(no*nv,w3(i),tl,w1(fai+i-1),no*nv)
          enddo
        elseif(s==2)then
#ifdef VAR_LSMPI
          call array_gather_tilesinfort(u2,w1,int(no*no*nv*nv,kind=long),infpar%master,[4,1,2,3])
          call ls_mpibcast(w1,int(nv*nv*no*no,kind=long),infpar%master,infpar%lg_comm)
          call dgemm('n','t',tl,nv*no,nv*no,0.5E0_realk,w2(faif),lead,w1,nv*no,0.0E0_realk,w3,lead)
          w1=0.0E0_realk
          do i=1,tl
            call dcopy(no*nv,w3(i),tl,w1(fai+i-1),no*nv)
          enddo
#endif
        endif
       
        ! preOmD [a i b j] -> =+ Omega [a b i j]
        if(s==4.or.s==3)then
          call array_reorder_4d(1.0E0_realk,w1,nv,no,nv,no,[1,3,2,4],1.0E0_realk,omega2%elm1)
        elseif(s==2)then
#ifdef VAR_LSMPI
          call lsmpi_local_reduction(w1,nv*nv*no*no,infpar%master)
          call array_scatteradd_densetotiled(omega2,1.0E0_realk,w1,int(no*no*nv*nv,kind=long),infpar%master,[1,3,2,4])
          !if(me==0)then
          !  call array_add(omega2,1.0E0_realk,w1,omega2%nelms,[1,3,2,4])
          !endif
#endif
        endif
      else
        call lsquit("ERROR(get_cnd_terms_mo):no valid scheme",-1)
      endif

      call mem_dealloc(w2)
      call mem_dealloc(w3)

  end subroutine get_cnd_terms_mo

  !> \brief Get the b2.2 contribution constructed in the kobayashi scheme after
  !the loop to avoid steep scaling ste  !> \author Patrick Ettenhuber
  !> \date December 2012
  subroutine get_B22_contrib_mo(sio4,t2,w1,w2,no,nv,nb,om2,s)
    implicit none
    !> the sio4 matrix from the kobayashi terms on input
    real(realk), intent(in) :: sio4(:)
    !> amplitudes
    !real(realk), intent(in) :: t2(*)
    type(array), intent(inout) :: t2
    !> some workspave
    real(realk), intent(inout) :: w1(:)
    real(realk), pointer :: w2(:)
    !> number of occupied, virutal and ao indices
    integer, intent(in) :: no,nv,nb
    !> residual to be updated
    !real(realk), intent(inout) :: om2(*)
    type(array), intent(inout) :: om2
    !> integer specifying the calc-scheme
    integer, intent(in) :: s
    integer :: nor,i,j,pos
    integer :: ml,l,tl,fai,lai
    integer :: tri,fri
    integer(kind=ls_mpik) :: nod,me,nnod,massa
    real(realk) :: nrm1,nrm2,nrm3,nrm4
    integer :: pos1, pos2, mv((nv*nv)/2),st
    me   = 0
    massa = 0
    nnod = 1
#ifdef VAR_LSMPI
    massa = infpar%master
    nnod = infpar%lg_nodtot
    me   = infpar%lg_mynum
#endif
      
    !Setting transformation variables for each rank
    !**********************************************
    call mo_work_dist(nv*nv,fai,tl)
    
    nor=no*(no+1)/2

    ! do contraction
    if(s==4.or.s==3.or.s==0)then
      w1=0.0E0_realk
      call dgemm('n','n',tl,nor,no*no,0.5E0_realk,t2%elm1(fai),nv*nv,sio4,no*no,0.0E0_realk,w1(fai),nv*nv)
    elseif(s==2.or.s==1)then
#ifdef VAR_LSMPI
      call mem_alloc(w2,tl*no*no)
      !if(me==0)call array_convert(t2,w1,t2%nelms)
      call array_gather_tilesinfort(t2,w1,int(nv*nv*no*no,kind=long),massa)
      do nod=1,nnod-1
        call mo_work_dist(nv*nv,fri,tri,nod)
        if(me==0)then
          do i=1,tri
            call dcopy(no*no,w1(fri+i-1),nv*nv,w2(i),tri)
          enddo
        endif
        if(me==0.or.me==nod)then
          call ls_mpisendrecv(w2,int(no*no*tri,kind=long),infpar%lg_comm,infpar%master,nod)
        endif
      enddo
      if(me==0)then
        do i=1,tl
          call dcopy(no*no,w1(fai+i-1),nv*nv,w2(i),tl)
        enddo
      endif
      w1=0.0E0_realk
      call dgemm('n','n',tl,nor,no*no,0.5E0_realk,w2,tl,sio4,no*no,0.0E0_realk,w1(fai),nv*nv)
      call mem_dealloc(w2)
#endif
    endif
   

    do j=no,1,-1
      do i=j,1,-1
        pos1=1+((i+j*(j-1)/2)-1)*nv*nv
        pos2=1+(i-1)*nv*nv+(j-1)*no*nv*nv
        if(j/=1) call dcopy(nv*nv,w1(pos1),1,w1(pos2),1)
      enddo
    enddo
    do j=no,1,-1
      do i=j,1,-1
          pos1=1+(i-1)*nv*nv+(j-1)*no*nv*nv
          pos2=1+(j-1)*nv*nv+(i-1)*no*nv*nv
          if(i/=j) call dcopy(nv*nv,w1(pos1),1,w1(pos2),1)
          call alg513(w1(pos1:nv*nv+pos1-1),nv,nv,nv*nv,mv,(nv*nv)/2,st)
      enddo
    enddo
    if(s==4.or.s==3.or.s==0)then
      call daxpy(int(no*no*nv*nv),1.0E0_realk,w1,1,om2%elm1,1)
    elseif(s==2.or.s==1)then
#ifdef VAR_LSMPI
      call lsmpi_local_reduction(w1,int(nv*nv*no*no,kind=long),infpar%master)
      call array_scatteradd_densetotiled(om2,1.0E0_realk,w1,int(no*no*nv*nv,kind=long),infpar%master)
      !if(me==0)then
      !  call array_add(om2,1.0E0_realk,w1,om2%nelms)
      !endif
#endif
    endif
  end subroutine get_B22_contrib_mo

  !> \brief subroutine to add contributions to the sio4 matrix which enters the
  !B2.2 term in the "non"-parallel region
  !> \author Patrick Ettenhuber
  !> \Date September 2012
  subroutine add_int_to_sio4(w0,w2,w3,no,nv,nb,fa,fg,la,lg,xo,sio4)
    implicit none
    !> workspace containing the pariteially transformed integrals ordered as I(i
    !gamma alpha j)
    real(realk),pointer :: w0(:)
    !> arbitrary workspace of correct size
    real(realk),pointer :: w2(:),w3(:)
    !> number of occupied, virutal and ao indices
    integer, intent(in) :: no,nv,nb
    !> first alpha and first gamma indices of the current loop
    integer, intent(in) :: fa,fg
    !> lengths of the alpha ang gamma batches in the currnet loop
    integer, intent(in) :: la,lg
    !> transformation matrix for t1 transformed integrals "Lambda p"
    real(realk),intent(in) :: xo(nb*no)
    !> sio4 storage space to update during the batched loops
    real(realk),pointer :: sio4(:)
    integer :: nor,pos,i,j

    nor=no*(no+1)/2

    ! (w3):I[alpha gamma i j] <- (w0):I[i gamma alpha j]
    call array_reorder_4d(1.0E0_realk,w0,no,lg,la,no,[2,3,1,4],0.0E0_realk,w2)
    ! (w2):I[alpha gamma i <= j] <- (w3):I[alpha gamma i j]
    do j=1,no
      do i=1,j
        pos=1+((i+j*(j-1)/2)-1)*la*lg
        call dcopy(la*lg,w2(1+(i-1)*la*lg+(j-1)*la*lg*no),1,w3(pos),1)
      enddo
    enddo
    ! (w3):I[ gamma i <= j alpha] <- (w2):I[alpha gamma i <= j]
    call array_reorder_3d(1.0E0_realk,w3,lg,la,nor,[2,3,1],0.0E0_realk,w2)
    ! (w2):I[ l i <= j alpha] <- (w3):Lambda^p [gamma l ]^T I[gamma i <= j alpha]
    call dgemm('t','n',no,nor*lg,la,1.0E0_realk,xo(fa),nb,w2,la,0.0E0_realk,w3,no)
    ! (sio4):I[ k l i <= j] <-+ (w2):Lambda^p [alpha k ]^T I[l i <= j alpha]^T
    call dgemm('t','t',no,nor*no,lg,1.0E0_realk,xo(fg),nb,w3,nor*no,1.0E0_realk,sio4,no)

  end subroutine add_int_to_sio4

  !> \brief calculate d term integral direct form echange and coulomb integrals
  !> \author Patrick Ettenhuber
  !> \date December 2012
  subroutine get_c_term_int_direct(w0,w1,w2,w3,no,nv,nb,fa,fg,la,lg,xo,&
             &yo,xv,yv,t2,tGammadij,omega2,t2jabi,s)
    implicit none
    !> empty workspace
    real(realk),intent(inout) :: w0(:),w2(:),w3(:)
    !> workpace containing the exchange ao integrals 
    real(realk),intent(inout) :: w1(:)
    !> partially transformed amplitudes
    real(realk),intent(inout) :: tGammadij(:)
    !> number of occupied, virutal and ao indices
    integer, intent(in) :: no,nv,nb
    !> first alpha and first gamma indices of the current loop
    integer, intent(in) :: fa,fg
    !> lengths of the alpha ang gamma batches in the currnet loop
    integer, intent(in) :: la,lg
    !> integer specifying the scheme
    integer, intent(in) :: s
    !> amplitudes
    !real(realk), intent(in) :: t2(:)
    type(array), intent(inout) :: t2
    !> residual to be updated
    !real(realk), intent(inout) :: omega2(:)
    type(array), intent(inout) :: omega2
    ! t2 reordered
    type(array), intent(inout) :: t2jabi
    !> Lambda transformation matrices
    real(realk), intent(in) :: xo(:),yo(:),xv(:),yv(:)
    integer :: nb2
    nb2=nb*nb
    !CALC C-TERM 
    !Get the integrals in the correct shape
    call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[1,3,2,4],0.0E0_realk,w0)
    call dcopy(la*lg*nb2,w0,1,w1,1)
    !(w2): I[l alpha gamma delta] = Lambda^p [beta l]^T * (w0):I[alpha gamma delta, beta]^T
    call dgemm('t','t',no,la*lg*nb,nb,1.0E0_realk,xo,nb,w0,la*lg*nb,0.0E0_realk,w2,no)
    !(w0): I[d l alpha gamma] = Lambda^h[delta d]^T * (w2):I[l alpha gamma, delta]^T
    call dgemm('t','t',nv,no*la*lg,nb,1.0E0_realk,yv,nb,w2,no*la*lg,0.0E0_realk,w0,nv)
    !(w3): t[ j a d l] <- t[a d l, j]^T
    if(s==1)then
      call array_convert(t2,w2,[2,3,4,1])
      !call array_convert(t2jabi,w2,t2jabi%nelms)
    else
      call array_reorder_4d(1.0E0_realk,t2%elm1,nv,nv,no,no,[4,1,2,3],0.0E0_realk,w2)
    endif
!   op_stop=omp_get_wtime()
    !(w2):C'[j a alpha gamma] = -0.5 * t[j a d l] * (w0):I[d l alpha gamma]
    call dgemm('n','n',no*nv,lg*la,nv*no,-0.5E0_realk,w2,no*nv,w0,nv*no,0.0E0_realk,w3,no*nv)
    !(w0):I[a gamma alpha delta] = Lambda^p[beta a]^T * I[alpha gamma delta, beta]^T 
    call dgemm('t','t',nv,la*lg*nb,nb,1.0E0_realk,xv,nb,w1,la*lg*nb,0.0E0_realk,w0,nv)
    !(w2):C'[j a alpha gamma] += Lamda:^h[delta j]^T * (w0):I[a gamma alpha, delta]^T
    call dgemm('t','t',no,nv*la*lg,nb,1.0E0_realk,yo,nb,w0,nv*la*lg,1.0E0_realk,w3,no)
    !(w3)t[alpha gamma b i] = Lambda^p[alpha k] * t[gamma b i, k]^T
    call dgemm('n','t',la,lg*nv*no,no,1.0E0_realk,xo(fa),nb,tGammadij,lg*nv*no,0.0E0_realk,w1,la)
    !(w0):C[j a b i] = -1* (w2):C'[j a alpha gamma] * (w3):t[alpha gamma b i]
    call dgemm('n','n',no*nv,no*nv,la*lg,-1.0E0_realk,w3,no*nv,w1,la*lg,0.0E0_realk,w2,no*nv)
    
    if(s==1)then
      call array_add(omega2,1.0E0_realk,w2,[2,3,4,1])
      call array_add(omega2,0.5E0_realk,w2,[2,3,1,4])
    else 
      call array_reorder_4d(1.0E0_realk,w2,no,nv,nv,no,[2,3,4,1],1.0E0_realk,omega2%elm1)
      call array_reorder_4d(0.5E0_realk,w2,no,nv,nv,no,[2,3,1,4],1.0E0_realk,omega2%elm1)
    endif
  end subroutine get_c_term_int_direct

  !> \brief calculate d term integral direct form echange and coulomb integrals
  !> \author Patrick Ettenhuber
  !> \date December 2012
  subroutine get_d_term_int_direct(w0,w1,w2,w3,no,nv,nb,fa,fg,la,lg,xo,&
             &yo,xv,yv,u2,uigcj,omega2,u2kcjb,s)
    implicit none
    !> workspace containing exchange ao integrals
    real(realk), intent(inout) :: w1(:)
    !> workspace containing coulomb ao integrals
    real(realk), intent(inout) :: w0(:)
    !> empty workspace
    real(realk), intent(inout) :: w2(:),w3(:)
    !> partially transformed u matrix
    real(realk), intent(inout) :: uigcj(:)
    !> u matrix
    !real(realk), intent(inout) :: u2(:)
    type(array), intent(inout) :: u2
    !> number of occupied, virutal and ao indices
    integer, intent(in) :: no,nv,nb
    !> first alpha and first gamma indices of the current loop
    integer, intent(in) :: fa,fg
    !> lengths of the alpha ang gamma batches in the currnet loop
    integer, intent(in) :: la,lg
    !> the residual to update
    !real(realk), intent(inout) :: omega2(:)
    type(array), intent(inout) :: omega2
    !> reordered u2 array, only saved in scheme=1
    type(array), intent(inout) :: u2kcjb
    !> scheme
    integer, intent(in) :: s
    !> Lambda transformation matrices
    real(realk), intent(in) :: xo(:),yo(:),xv(:),yv(:)
    integer :: o2v
    o2v=no*no*nv

    call array_reorder_4d(-1.0E0_realk,w1,la,nb,lg,nb,[3,2,1,4],1.0E0_realk,w0)
       
    !Construct L gamma delta alpha beta from g alpha beta gamma delta with
    !(w0) L[gamma delta alpha beta] = 2* I^C[gamma delta alpha beta] - I^E[gamma delta alpha beta]
    
    !(w2):L[gamma delta alpha k] = (w0):L[gamma delta alpha beta] * Lambda^p[beta k]
    call dgemm('n','n',lg*la*nb,no,nb,1.0E0_realk,w0,lg*nb*la,xo,nb,0.0E0_realk,w2,lg*nb*la)
    !(w0):L[gamma alpha k delta] <- (w2):L[gamma delta alpha k]
    call array_reorder_4d(1.0E0_realk,w2,lg,nb,la,no,[1,3,4,2],0.0E0_realk,w0)
    !(w3):L[gamma alpha k c] = L[gamma alpha k delta] * Lambda^h[delta c]
    call dgemm('n','n',lg*la*no,nv,nb,1.0E0_realk,w0,lg*no*la,yv,nb,0.0E0_realk,w3,lg*no*la)
    !(w2):u[k c j b] <- u[b c j k]
    if(s==1)then
      !call array_convert(u2kcjb,w2,u2kcjb%nelms)
      call array_convert(u2,w2,[2,4,3,1])
    else
      call array_reorder_4d(1.0E0_realk,u2%elm1,nv,nv,no,no,[4,1,3,2],0.0E0_realk,w2)
    endif
    !(w0)D[gamma alpha j b] = 0.5 * (w3):L[gamma alpha k c] * (w2):u[k c j b]
    call dgemm('n','n',lg*la,no*nv,no*nv,0.5E0_realk,w3,lg*la,w2,no*nv,0.0E0_realk,w0,lg*la)
    !(w2):D'[i alpha j b] = Lambda^h[gamma i]^T * (w0):D[gamma alpha j b]
    call dgemm('t','n',no,la*no*nv,lg,1.0E0_realk,yo(fg),nb,w0,lg,0.0E0_realk,w2,no) 
    !(w3):D'[alpha b i j] <- (w2):D[i alpha j b]
    call array_reorder_4d(1.0E0_realk,w2,no,la,no,nv,[2,4,1,3],0.0E0_realk,w3)
    !Omega^D += Lambda^p[alpha a]^T * (w3):D'[alpha b i j]    <---------TO COME
    if(s==1)then
      call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w3,la,0.0E0_realk,w2,nv) 
      call array_add(omega2,1.0E0_realk,w2)
    else
      call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w3,la,1.0E0_realk,omega2%elm1,nv) 
    endif
    
    !(w2)u[i gamma a alpha] = u[i gamma a l] * Lambda^p[alpha l]^T
    call dgemm('n','t',no*nv*lg,la,no,1.0E0_realk,uigcj,no*nv*lg,xo(fa),nb,0.0E0_realk,w2,no*nv*lg)
    !(w3):u[i a gamma alpha] <- (w2):u[i gamma a alpha]
    call array_reorder_4d(1.0E0_realk,w2,no,lg,nv,la,[1,3,2,4],0.0E0_realk,w3)
    !(w2)D''[iajb] = 0.5 * (w3):u[i a gamma alpha] * (w0):D[gamma alpha j b]
    call dgemm('n','n',no*nv,no*nv,la*lg,0.5E0_realk,w3,no*nv,w0,la*lg,0.0E0_realk,w2,no*nv)
    if(s==1)then
      call array_add(omega2,1.0E0_realk,w2,[2,4,1,3])
    else
      call array_reorder_4d(1.0E0_realk,w2,no,nv,no,nv,[2,4,1,3],1.0E0_realk,omega2%elm1)
    endif
  end subroutine get_d_term_int_direct

  !> \brief calculate a and b terms in a kobayashi fashion
  !> \author Patrick Ettenhuber
  !> \date December 2012
  subroutine get_a22_and_prepb22_terms_ex(w0,w1,w2,w3,tpl,tmi,no,nv,nb,fa,fg,la,lg,&
  &xo,yo,xv,yv,om2,sio4,s)       
    implicit none
    !> workspace with exchange integrals
    real(realk),intent(inout) :: w0(:)
    !> empty workspace of correct sizes
    real(realk),intent(inout) :: w1(:),w2(:),w3(:)
    !> the t+ and t- combinations with a value of the amplitudes with the
    !diagonal elements divided by two
    real(realk),intent(inout) :: tpl(:),tmi(:)
    !> number of occupied, virutal and ao indices
    integer, intent(in) :: no,nv,nb
    !> first alpha and first gamma indices of the current loop
    integer, intent(in) :: fa,fg
    !> lengths of the alpha ang gamma batches in the currnet loop
    integer, intent(in) :: la,lg
    !> the doubles residual to update
    !real(realk), intent(inout) :: om2(:)
    type(array), intent(inout) :: om2
    !> the lambda transformation matrices
    real(realk), intent(in)    :: xo(:),yo(:),xv(:),yv(:)
    !> the sio4 matrix to calculate the b2.2 contribution
    real(realk),intent(inout) ::sio4(:)
    !> scheme
    integer,intent(in) :: s
    integer :: goffs,aoffs,tlen,tred,nor,nvr

    nor=no*(no+1)/2
    nvr=nv*(nv+1)/2
    
    !Determine the offsets in the alpha and gamma indices, which arise from
    !non uniform batch distributions, e.g. consider:
    !Case 1:   gamma _____
    !               |\    |
    !               |_\___|
    !  alpha        |  \  |
    !               |___\_| here an offset in gamma for the second gamma batch is
    !               needed, since the elements to be considered do not
    !               begin with the first element in that batch
    goffs=0
    if(fa-fg>0)goffs=fa-fg
    !Case 2:   gamma ______
    !               |\ |   |
    !               | \|   |
    !  alpha        |  |\  |
    !               |__|_\_| here an offset in alpha for the second alpha batch is
    !               needed, since the elements to be considered do not
    !               begin with the first element in that batch
    aoffs=0
    if(fg-fa>0)aoffs=fg-fa
    
    !Determine the dimension of the triangular part in the current batch
    !and the total number of elements (tred) in the triangular plus
    !rectangular parts
    tred=0
    tlen=min(min(min(la,lg),fg+lg-fa),fa+la-fg)
    
    !calculate amount of triangular elements in the current batch
    if(fa+la-1>=fg)tred= tlen*(tlen+1)/2
    !add the rectangular contribution if lg is larger than tlen 
    if(fa>=fg.and.fg+lg-fa-tlen>0)tred=tred+(fg+lg-fa-tlen)*la
    if(fa<fg)tred=tred+lg*aoffs
    if(fa<fg.and.fa+la>fg) tred=tred+(lg-tlen)*(la-aoffs)
    !if only rectangular
    if(tlen<=0)then
      tlen=0
      aoffs=0
      goffs=0
      tred=la*lg
    endif
    !!SYMMETRIC COMBINATION
    !(w0):I+ [delta alpha<=gamma beta] <= (w1):I [alpha beta gamma delta] + (w1):I[alpha delta gamma beta]
    call get_I_plusminus_le(w0,w1,w2,'p',fa,fg,la,lg,nb,tlen,tred,goffs)
    call lsmpi_poke()
    !(w2):I+ [delta alpha<=gamma c] = (w0):I+ [delta alpha<=gamma beta] * Lambda^h[beta c]
    call dgemm('n','n',nb*tred,nv,nb,1.0E0_realk,w0,nb*tred,yv,nb,0.0E0_realk,w2,nb*tred)
    call lsmpi_poke()
    !(w0):I+ [alpha<=gamma c d] = (w2):I+ [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
    call dgemm('t','n',tred*nv,nv,nb,1.0E0_realk,w2,nb,yv,nb,0.0E0_realk,w0,nv*tred)
    call lsmpi_poke()
    !(w2):I+ [alpha<=gamma c>=d] <= (w0):I+ [alpha<=gamma c d] 
    call get_I_cged(w2,w0,tred,nv)
    call lsmpi_poke()
    !(w3.1):sigma+ [alpha<=gamma i>=j] = (w2):I+ [alpha<=gamma c>=d] * (w0):t+ [c>=d i>=j]
    call dgemm('n','n',tred,nor,nvr,0.5E0_realk,w2,tred,tpl,nvr,0.0E0_realk,w3,tred)
    call lsmpi_poke()
    
    
    !!ANTI-SYMMETRIC COMBINATION
    !(w0):I- [delta alpha<=gamma beta] <= (w1):I [alpha beta gamma delta] + (w1):I[alpha delta gamma beta]
    call get_I_plusminus_le(w0,w1,w2,'m',fa,fg,la,lg,nb,tlen,tred,goffs)
    call lsmpi_poke()
    !(w2):I- [delta alpha<=gamma c] = (w0):I- [delta alpha<=gamma beta] * Lambda^h[beta c]
    call dgemm('n','n',nb*tred,nv,nb,1.0E0_realk,w0,nb*tred,yv,nb,0.0E0_realk,w2,nb*tred)
    call lsmpi_poke()
    !(w0):I- [alpha<=gamma c d] = (w2):I- [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
    call dgemm('t','n',tred*nv,nv,nb,1.0E0_realk,w2,nb,yv,nb,0.0E0_realk,w0,nv*tred)
    call lsmpi_poke()
    !(w2):I- [alpha<=gamma c<=d] <= (w0):I- [alpha<=gamma c d] 
    call get_I_cged(w2,w0,tred,nv)
    call lsmpi_poke()
    !(w3.2):sigma- [alpha<=gamma i<=j] = (w2):I- [alpha<=gamma c>=d] * (w0):t- [c>=d i>=j]
    call dgemm('n','n',tred,nor,nvr,0.5E0_realk,w2,tred,tmi,nvr,0.0E0_realk,w3(tred*nor+1),tred)
    call lsmpi_poke()
    
    !COMBINE THE TWO SIGMAS OF W3 IN W2
    !(w2):sigma[alpha<=gamma i<=j]=0.5*(w3.1):sigma+ [alpha<=gamma i<=j] + 0.5*(w3.2):sigma- [alpha <=gamm i<=j]
    !(w2):sigma[alpha>=gamma i<=j]=0.5*(w3.1):sigma+ [alpha<=gamma i<=j] - 0.5*(w3.2):sigma- [alpha <=gamm i<=j]
    call combine_and_transform_sigma(om2,w0,w2,w3,xv,xo,sio4,nor,&
    &tlen,tred,fa,fg,la,lg,no,nv,nb,goffs,aoffs,s)  
  end subroutine get_a22_and_prepb22_terms_ex

  !> \brief Combine sigma matrixes in symmetric and antisymmetric combinations 
  !> \author Patrick Ettenhuber
  !> \date October 2012
  subroutine combine_and_transform_sigma(omega,w0,w2,w3,xvirt,xocc,sio4,nor,&
  &tlen,tred,fa,fg,la,lg,no,nv,nb,goffs,aoffs,s)
    implicit none
    !\> omega should be the residual matrix which contains the second parts
    !of the A2 and B2 term
    !real(realk),intent(inout) :: omega(nv*nv*no*no)
    type(array),intent(inout) :: omega
    !> w0 is just some workspace on input
    real(realk),intent(inout) :: w0(:)
    !> w2 is just some workspace on input
    real(realk),intent(inout) :: w2(:)
    !> w3 contains the symmetric and antisymmetric combinations 
    real(realk),intent(inout) :: w3(:)
    !> sio4 are the reduced o4 integrals whic are used to calculate the B2.2
    !contribution after the loop, update them in the loops
    real(realk),intent(inout) :: sio4(:)
    !> Lambda p virutal part
    real(realk),intent(in) :: xvirt(:)
    !> Lambda p occupied part
    real(realk),intent(in) :: xocc(:)
    !> number of reduced occupied indices 
    integer,intent(in)::nor
    !> total number of upper triangular elements in the batch
    integer,intent(in)::tlen
    !> first element of alpha and gamma in the current batch
    integer,intent(in)::fa,fg
    !> length of alpha and gamma in the current batch
    integer,intent(in)::la,lg
    !> number of triangular elements in the batch
    integer,intent(in)::tred
    !>number of occupied, virtual and ao indices
    integer,intent(in)::no,nv,nb
    !>offsets in the currnt alpha and gamma batches to get to the first upper
    !triangular element
    integer,intent(in)::goffs,aoffs
    !> scheme
    integer,intent(in)::s
    !> the doubles amplitudes
    !real(realk),intent(in) :: amps(nv*nv*no*no)
    !type(array),intent(in) :: amps
    real(realk) :: scaleitby
    integer ::occ,gamm,alpha,pos,pos2,pos21,nel2cp,case_sel,full1,full2,offset1,offset2,ttri
    integer :: l1,l2,i,j,lsa,lsg,gamm_i_b,a,b,full1T,full2T,tsq,jump,dim_big,dim_small,ft1,ft2,ncph
    logical :: second_trafo_step
    real(realk),pointer :: dumm(:)
    integer :: mv((nv*nv)/2),st,pos1

    scaleitby=0.5E0_realk
    second_trafo_step=.false.
    lsa=fa+la-1
    lsg=fg+lg-1
    !building the sigma matrices is a complicated matter, due to non-uniform
    !batchsizes. A distiction between the different possible cases has to be
    !made and the elements copied respectively
    case_sel=0
    if(fa>=fg)then
      if(lsa>=lsg)then
        !case 0 is a triagular submatrix:
        !                                full2=tlen
        ! ---|                          |----|
        ! \  |                          | \  |  full1=tlen
        !  \ |   build full from that   |  \ |
        !   \|                          |---\|
        case_sel=1
        full1=tlen
        full2=tlen
        full1T=tlen
        full2T=tlen
        l1=0
        l2=0
      elseif(lsa<lsg)then
        !case 1 is a  submatrix of the type:
        !                                   full2=lg-goffs
        ! ------|                          |------|
        ! \     |                          | \    |
        !  \    |   build full from that   |  \   |full1=la
        !   \---|                          |___\__|
        !                                  |   |0 |
        !            full1T=lg-goffs-tlen  |___|__|
        case_sel=2                         !full2T=la 
        full1=la
        full2=lg-goffs
        full1T=lg-goffs-tlen
        full2T=la
        l1=la
        l2=0
        second_trafo_step=.true.
      endif
    else
      if(lsa>=lsg)then
        !case 2 is a  submatrix of the type:
        !                                 full2=lg
        ! |---|                         |------|
        ! |   |                         |0 |   |
        ! |   |  build full from that   |__|   | full1=aoffs+tlen
        !  \  |                         |  |\  |
        !   \ |                         |  | \ |
        !    \|     full1T=lg           |__|__\|
        !                               full2T=aoffs
        case_sel=3
        full1=aoffs+tlen
        full2=tlen
        full1T=tlen
        full2T=aoffs
        l1=aoffs
        l2=aoffs
        second_trafo_step=.true.
      elseif(lsa<lsg)then
        if(lsa>=fg)then
          !case 3 is a  submatrix of the type: 
          !                                       full2=lg
          ! |-----|                         |--------|
          ! |     |                         |0 |     |
          ! |     |  build full from that   |__|     |
          !  \    |                         |  |\    |  full1 = la
          !   \   |                         |  | \   |
          !    \--|                         |  |__\__|
          !                     full1T=lg   |     | 0|
          !                                 |_____|__|
          !                                    full2T = la
          case_sel=4
          full1=la
          full2=lg
          full1T=lg
          full2T=la
          l1=aoffs
          l2=aoffs
          second_trafo_step=.true.
        else
          !case 4 is a  submatrix of dimensions la lg, and can thus be transformed in
          !full anyways
          case_sel=5
          full1=la
          full2=lg
          full1T=lg
          full2T=la
          l1=la
          l2=0
          second_trafo_step=.true.
        endif
      endif
    endif
    if(case_sel==0)call lsquit("ERROR(combine_and_transform_sigma):case not known",DECinfo%output)

    !print *,"-------------------------------------------------------------------------------------"
    !print *,"       case        dim1        dim2            la          lg        tlen        tred"
    !print *,case_sel,full1,full2,la,lg,tlen,tred
    !print *,"-------------------------------------------------------------------------------------"



    !Zero the elements to update for testing, not needed in a performance
    !implementation
    pos=nor*full1*full2
    if(second_trafo_step)pos=pos+nor*full1T*full2T
    w0(1:pos)=0.0E0_realk

    !set required variables
    ttri=tlen*(tlen+1)/2
    tsq = tlen * tlen
    pos=1
    occ=1
    dim_big=full1*full2
    dim_small=full1T*full2T

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP& SHARED(w0,w3,case_sel,nor,goffs,lg,la,full1,full1T,ttri,tred,&
    !$OMP& full2,full2T,tlen,l1,second_trafo_step,aoffs,dim_big,dim_small,l2)&
    !$OMP& PRIVATE(occ,gamm,gamm_i_b,pos,nel2cp,pos2,jump,ft1,ft2,ncph,pos21)
    !$OMP DO
    do occ=1,nor
      do gamm=1,lg-goffs
        gamm_i_b=gamm+goffs
        !SYMMETRIC COMBINATION OF THE SIGMAS
        if(case_sel==3.or.case_sel==4)then
          pos=1+(gamm-1)*full1+(occ-1)*full1*full2
        else
          pos=1+(gamm+aoffs-1)*full1+(occ-1)*full1*full2
        endif
        if(gamm>tlen)then
          !get the elements from the rectangular part of the batch
          !print *,"getel from rect"
          nel2cp=l1
          pos2=1+ttri+(gamm-tlen-1)*(la-aoffs)+(occ-1)*((la-aoffs)*(lg-goffs-tlen)+ttri)
          if(case_sel==4)then
            nel2cp=nel2cp+tlen
            pos2 = pos2 + tlen * aoffs + (gamm-tlen-1) * (la-tlen)+(occ-1)*(aoffs*lg)
          endif
          if(second_trafo_step)then
            jump = full1T
            ft1 = full1T
            ft2 = full2T
          else
            jump=full1
            ft1=full1
            ft2=full2
          endif

        else
          !get the elements from the triangular part of the batch
          pos2=1+(gamm*(gamm-1)/2)+(gamm-1)*aoffs+(occ-1)*tred
          nel2cp=l2+gamm
          jump = full1
          ft1=full1
          ft2=full2
        endif

        call dcopy(nel2cp,w3(pos2),1,w0(pos),1)
        !get corresponding position in sigma- and add to output
        pos21=pos2+tred*nor
        call daxpy(nel2cp,1.0E0_realk,w3(pos21),1,w0(pos),1)
    
        !ANTI-SYMMETRIC COMBINATION OF THE SIGMAS
        pos=gamm+aoffs+(occ-1)*ft1*ft2
        if(second_trafo_step.and.gamm>tlen) pos=pos+full1*full2*nor-tlen
        if(case_sel==3.or.case_sel==4)then
          !fill diagonal part
          if(gamm>tlen)then
            ncph=0
          else
            ncph=gamm
          endif
          call daxpy(ncph,1.0E0_realk,w3(pos2+aoffs),1,w0(aoffs+gamm+(occ-1)*full1*full2),full1)
          call daxpy(ncph,-1.0E0_realk,w3(pos21+aoffs),1,w0(aoffs+gamm+(occ-1)*dim_big),full1)
          !fill small matrix
          if(gamm>tlen)then
            ncph=nel2cp
          else
            ncph=nel2cp-gamm
          endif
          call daxpy(ncph,1.0E0_realk,w3(pos2),1,w0(gamm+(occ-1)*full1T*full2T+dim_big*nor),full1T)
          call daxpy(ncph,-1.0E0_realk,w3(pos21),1,w0(gamm+(occ-1)*full1T*full2T+dim_big*nor),full1T)
        else
          call daxpy(nel2cp,1.0E0_realk,w3(pos2),1,w0(pos),jump)
          call daxpy(nel2cp,-1.0E0_realk,w3(pos21),1,w0(pos),jump)
        endif
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call lsmpi_poke()


    !Print the individual contributions
    !if(.false.)then
    !  do occ=1,nor
    !    print *,"symmetric"
    !    print *,w3(1+(occ-1)*tred:tred+(occ-1)*tred)
    !    !print *,"antisymmetric"
    !    !print *,w3(nor*tred+1+(occ-1)*tred:nor*tred+tred+(occ-1)*tred)
    !    print *,"a g -> 1. and l.",1+(occ-1)*full1*full2,full1*full2+(occ-1)*full1*full2
    !    do alpha=1,full1
    !      print *,(w0(alpha+(gamm-1)*full1+(occ-1)*full1*full2),gamm=1,full2)
    !    enddo
    !    print *,""
    !    if(second_trafo_step)then
    !      print *,"g a -> 1. and l.",dim_big*nor+1+(occ-1)*full1T*full2T,dim_big*nor+full1T*full2T+(occ-1)*full1T*full2T
    !      do gamm=1,full1T
    !        print *,(w0(dim_big*nor+gamm+(alpha-1)*full1T+(occ-1)*full1T*full2T),alpha=1,full2T)
    !      enddo
    !      print *,""
    !    endif
    !    print *,""
    !  enddo
    !endif

    !add up the contributions for the sigma [ alpha gamma ] contributions
    ! get the order sigma[ gamma i j alpha ]
    call mat_transpose_pl(full1,full2*nor,1.0E0_realk,w0,0.0E0_realk,w2)
    call lsmpi_poke()
    !transform gamma -> b
    call dgemm('t','n',nv,nor*full1,full2,1.0E0_realk,xvirt(fg+goffs),nb,w2,full2,0.0E0_realk,w3,nv)
    call lsmpi_poke()
    !transform alpha -> a , order is now sigma [ a b i j]
    call dgemm('t','t',nv,nv*nor,full1,1.0E0_realk,xvirt(fa),nb,w3,nor*nv,0.0E0_realk,w2,nv)
    call lsmpi_poke()


    ! add up contributions in the residual with keeping track of i<j
    !OMP PARALLEL DEFAULT(NONE) SHARED(no,w2,omega,nv,scaleitby)&
    !OMP& PRIVATE(i,j,pos)
    !OMP DO 
   ! do i=1,no
   !   do j=1,no
   !     if(i<=j)then
   !       pos=1+((i+j*(j-1)/2)-1)*nv*nv
   !       call mat_transpose_pl(nv,nv,scaleitby,w2(pos:pos+nv*nv-1),&
   !            &1.0E0_realk,omega%elm1(1+(i-1)*nv*nv+(j-1)*no*nv*nv:))
   !     endif
   !     if(i>j)then
   !       pos=1+((j+i*(i-1)/2)-1)*nv*nv
   !       call daxpy(nv*nv,scaleitby,w2(pos),1,&
   !                 &omega%elm1(1+(i-1)*nv*nv+(j-1)*no*nv*nv),1)
   !     endif
   !   enddo
   ! enddo
    !dumm = 0.0E0_realk
    do j=no,1,-1
      do i=j,1,-1
        pos1=1+((i+j*(j-1)/2)-1)*nv*nv
        pos2=1+(i-1)*nv*nv+(j-1)*no*nv*nv
        if(j/=1) call dcopy(nv*nv,w2(pos1),1,w2(pos2),1)
      enddo
    enddo
    do j=no,1,-1
      do i=j,1,-1
          pos1=1+(i-1)*nv*nv+(j-1)*no*nv*nv
          pos2=1+(j-1)*nv*nv+(i-1)*no*nv*nv
          if(i/=j) call dcopy(nv*nv,w2(pos1),1,w2(pos2),1)
          call alg513(w2(pos1:nv*nv+pos1-1),nv,nv,nv*nv,mv,(nv*nv)/2,st)
      enddo
    enddo
    if(s==4.or.s==3.or.s==0)then
      call daxpy(int(no*no*nv*nv),scaleitby,w2,1,omega%elm1,1)
    elseif(s==2.or.s==1)then
      call array_add(omega,scaleitby,w2)
    endif
    !OMP END DO
    !OMP END PARALLEL
    call lsmpi_poke()

    !If the contributions are split in terms of the sigma matrix this adds the
    !sigma [gamma alpha] contributions
    if(second_trafo_step)then
      !get the order sigma [gamma i j alpha]
      pos2 = 1+full1*full2*nor
      if(case_sel==3.or.case_sel==4)then
        l1=fa
        l2=fg
      else
        l1=fa
        l2=fg+goffs+tlen
      endif
      call mat_transpose_pl(full1T,full2T*nor,1.0E0_realk,w0(pos2:full1T*full2T*nor+pos2-1),0.0E0_realk,w2)
      call lsmpi_poke()
      !transform gamma -> a
      call dgemm('t','n',nv,nor*full1T,full2T,1.0E0_realk,xvirt(l1),nb,w2,full2T,0.0E0_realk,w3,nv)
      call lsmpi_poke()
      !transform alpha -> , order is now sigma[b a i j]
      call dgemm('t','t',nv,nv*nor,full1T,1.0E0_realk,xvirt(l2),nb,w3,nor*nv,0.0E0_realk,w2,nv)
      call lsmpi_poke()
      ! add up contributions in the residual with a and b interchanged compared
      ! to the previous contribution
      !OMP PARALLEL DEFAULT(NONE) SHARED(no,w2,omega,nv,scaleitby)&
      !OMP& PRIVATE(i,j,pos)
      !OMP DO 
      !do i=1,no
      !  do j=1,no
      !    if(i<=j)then
      !      pos=1+((i+j*(j-1)/2)-1)*nv*nv
      !      call mat_transpose_pl(nv,nv,scaleitby,w2(pos:pos+nv*nv-1),&
      !            &1.0E0_realk,omega%elm1(1+(i-1)*nv*nv+(j-1)*no*nv*nv:))
      !    endif
      !    if(i>j)then
      !      pos=1+((j+i*(i-1)/2)-1)*nv*nv
      !      call daxpy(nv*nv,scaleitby,w2(pos),1,&
      !                &omega%elm1(1+(i-1)*nv*nv+(j-1)*no*nv*nv),1)
      !    endif
      !  enddo
      !enddo
      do j=no,1,-1
        do i=j,1,-1
          pos1=1+((i+j*(j-1)/2)-1)*nv*nv
          pos2=1+(i-1)*nv*nv+(j-1)*no*nv*nv
          if(j/=1) call dcopy(nv*nv,w2(pos1),1,w2(pos2),1)
        enddo
      enddo
      do j=no,1,-1
        do i=j,1,-1
            pos1=1+(i-1)*nv*nv+(j-1)*no*nv*nv
            pos2=1+(j-1)*nv*nv+(i-1)*no*nv*nv
            if(i/=j) call dcopy(nv*nv,w2(pos1),1,w2(pos2),1)
            call alg513(w2(pos1:nv*nv+pos1-1),nv,nv,nv*nv,mv,(nv*nv)/2,st)
        enddo
      enddo
      if(s==4.or.s==3.or.s==0)then
        call daxpy(no*no*nv*nv,scaleitby,w2,1,omega%elm1,1)
      elseif(s==2.or.s==1)then
        call array_add(omega,scaleitby,w2)
      endif
      !OMP END DO
      !OMP END PARALLEL
      call lsmpi_poke()
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Construct the B2 term from the intermediates in w0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !add up the contributions for the sigma [ alpha gamma ] contributions
    ! get the order sigma[ gamma i j alpha ]
    call mat_transpose_pl(full1,full2*nor,1.0E0_realk,w0,0.0E0_realk,w2)
    call lsmpi_poke()
    !transform gamma -> l
    call dgemm('t','n',no,nor*full1,full2,1.0E0_realk,xocc(fg+goffs),nb,w2,full2,0.0E0_realk,w3,no)
    call lsmpi_poke()
    !transform alpha -> a , order is now sigma [ k l i j]
    call dgemm('t','t',no,no*nor,full1,1.0E0_realk,xocc(fa),nb,w3,nor*no,1.0E0_realk,sio4,no)
    call lsmpi_poke()



    !If the contributions are split in terms of the sigma matrix this adds the
    !sigma [gamma alpha i j] contributions
    if(second_trafo_step)then
      !get the order sigma [gamma i j alpha]
      pos2 = 1+full1*full2*nor
      if(case_sel==3.or.case_sel==4)then
        l1=fa
        l2=fg
      else
        l1=fa
        l2=fg+goffs+tlen
      endif
      call mat_transpose_pl(full1T,full2T*nor,1.0E0_realk,w0(pos2:full1T*full2T*nor+pos2-1),0.0E0_realk,w2)
      call lsmpi_poke()
      !transform gamma -> l
      call dgemm('t','n',no,nor*full1T,full2T,1.0E0_realk,xocc(l1),nb,w2,full2T,0.0E0_realk,w3,no)
      call lsmpi_poke()
      !transform alpha -> k, order is now sigma[k l i j]
      call dgemm('t','t',no,no*nor,full1T,1.0E0_realk,xocc(l2),nb,w3,nor*no,1.0E0_realk,sio4,no)
      call lsmpi_poke()
    endif


  end subroutine combine_and_transform_sigma


  !> \brief Construct integral matrix which uses symmetries in the virtual
  !> indices
  !> \author Patrick Ettenhuber
  !> \date October 2012
  subroutine get_I_cged(Int_out,Int_in,m,nv)
    implicit none
    !> integral output with indices reduced to c>=d
    real(realk),intent(inout)::Int_out(:)
    !>full integral input m*c*d
    real(realk),intent(in) :: Int_in(:)
    !> leading dimension m and virtual dimension
    integer,intent(in)::m,nv
    integer ::d,pos,pos2,a,b,c,cged
    !OMP PARALLEL DEFAULT(NONE) SHARED(int_in,int_out,m,nv)&
    !OMP& PRIVATE(pos,pos2,d)
    pos =1
    !OMP DO
    do d=1,nv
      pos2=1+(d-1)*m+(d-1)*nv*m
      !print *,pos2,pos
      call dcopy(m*(nv-d+1),Int_in(pos2),1,Int_out(pos),1)
      pos=pos+m*(nv-d+1)
    enddo
    !OMP END DO
    !OMP END PARALLEL
  end subroutine get_I_cged


  !> \brief Construct symmetric and antisymmentric combinations of an itegral
  !matrix 
  !> \author Patrick Ettenhuber
  !> \date October 2012
  subroutine get_I_plusminus_le(w0,w1,w2,op,fa,fg,la,lg,nb,tlen,tred,goffs)
    implicit none
    !> blank workspace
    real(realk),intent(inout) :: w0(:),w2(:)
    !> workspace containing the integrals
    real(realk),intent(in) :: w1(:)
    !> integer specifying the first element in alpha and gamma batch
    integer,intent(in) :: fa,fg
    !> integer specifying the length of alpha and gamma batches
    integer,intent(in) :: la,lg
    !> number of upper triangular elements in the current batch
    integer,intent(in) :: tred
    !> number of elements in triangular part
    integer,intent(in) :: tlen
    !> number of ao indices
    integer,intent(in) :: nb
    !> offset in the gamma batch
    integer,intent(in) :: goffs
    !> character specifying which operation to use
    character, intent(in) :: op
    integer :: i,alpha,beta,gamm,delta,cagc,cagi,bs,bctr
    integer :: alpha_b,beta_b,gamm_b,delta_b,elements,aleg
    integer :: globg, globa, loca,eldiag,elsqre,nbnb,nrnb
    real(realk) ::chk,chk2,el
    real(realk),pointer :: trick(:,:,:)
    logical :: modb
    bs=int(sqrt(((8.0E6_realk)/1.6E1_realk)))
    !bs=5
    !print *,"block size",bs,(bs*bs*8)/1024.0E0_realk
    nbnb=(nb/bs)*bs
    modb=(mod(nb,bs)>0)
    bctr = bs-1
    cagi=tred

    call ass_D1to3(w2,trick,[nb,nb,cagi])
    if(op=='p')then
      call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
      aleg=0
      do gamm=0,lg-1
        do alpha=0,la-1
          if(fa+alpha<=fg+gamm)then
            !aleg = (alpha+(gamm*(gamm+1))/2) 
            eldiag = aleg*nb*nb
            elsqre = alpha*nb*nb+gamm*nb*nb*la
            !print *,alpha,gamm,1+eldiag,nb*nb+eldiag,1+elsqre,nb*nb+elsqre,aleg,cagi,nb*nb
            if(fa+alpha==fg+gamm)   call dscal(nb*nb,0.5E0_realk,w2(1+elsqre),1)
            call dcopy(nb*nb,w2(1+elsqre),1,w2(1+eldiag),1)
            !$OMP PARALLEL PRIVATE(el,delta_b,beta_b,beta,delta)&
            !$OMP& SHARED(bs,bctr,trick,nb,aleg,nbnb,modb)&
            !$OMP& DEFAULT(NONE)
            if(nbnb>0)then
              !$OMP DO
              do delta_b=1,nbnb,bs
                do beta_b=delta_b+bs,nbnb,bs
                  do delta=0,bctr
                    do beta=0,bctr
                      el=trick(beta+beta_b,delta_b+delta,aleg+1)
                         trick(beta+beta_b,delta_b+delta,aleg+1)=&
                        &trick(beta+beta_b,delta_b+delta,aleg+1) + &
                      &trick(delta_b+delta,beta+beta_b,aleg+1)
                       trick(delta_b+delta,beta+beta_b,aleg+1)=&
                      &trick(delta_b+delta,beta+beta_b,aleg+1) + &
                      &el
                    enddo
                  enddo
                enddo
              enddo
              !$OMP END DO NOWAIT
            endif
            if(nbnb>0.and.modb)then
              !$OMP DO
              do delta_b=1,nbnb,bs
                do delta=0,bctr
                  do beta=nbnb+1,nb
                    el=trick(beta,delta_b+delta,aleg+1)
                       trick(beta,delta_b+delta,aleg+1)=&
                      &trick(beta,delta_b+delta,aleg+1) + &
                    &trick(delta_b+delta,beta,aleg+1)
                     trick(delta_b+delta,beta,aleg+1)=&
                    &trick(delta_b+delta,beta,aleg+1) + &
                    &el
                  enddo
                enddo
              enddo
              !$OMP END DO NOWAIT
            endif
            if(nbnb>0)then
              !$OMP DO
              do delta_b=1,nbnb,bs
                do delta=0,bctr
                  do beta=delta+1,bctr
                    el=trick(beta+delta_b,delta_b+delta,aleg+1)
                       trick(beta+delta_b,delta_b+delta,aleg+1)=&
                      &trick(beta+delta_b,delta_b+delta,aleg+1) + &
                    &trick(delta_b+delta,beta+delta_b,aleg+1)
                     trick(delta_b+delta,beta+delta_b,aleg+1)=&
                    &trick(delta_b+delta,beta+delta_b,aleg+1) + &
                    &el
                  enddo
                   trick(delta+delta_b,delta_b+delta,aleg+1)=&
                  &trick(delta+delta_b,delta_b+delta,aleg+1) + &
                  &trick(delta_b+delta,delta+delta_b,aleg+1)
                enddo
              enddo
              !$OMP END DO NOWAIT
            endif
            !$OMP END PARALLEL 
            if(modb)then
              do delta=nbnb+1,nb
                do beta=delta+1,nb
                  el=trick(beta,delta,aleg+1)
                     trick(beta,delta,aleg+1)=&
                    &trick(beta,delta,aleg+1) + &
                  &trick(delta,beta,aleg+1)
                   trick(delta,beta,aleg+1)=&
                  &trick(delta,beta,aleg+1) + &
                  &el
                enddo
                trick(delta,delta,aleg+1)=&
                &trick(delta,delta,aleg+1) + &
                &trick(delta,delta,aleg+1)
              enddo
            endif
            aleg=aleg+1
          endif
        enddo
      enddo
      call array_reorder_3d(1.0E0_realk,w2,nb,nb,cagi,[2,3,1],0.0E0_realk,w0)
    endif


    !Calculate the antisymmetric combination in the same way
    if(op=='m')then
      call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
      aleg=0
      do gamm=0,lg-1
        do alpha=0,la-1
          if(fa+alpha<=fg+gamm)then
            !aleg = (alpha+(gamm*(gamm+1))/2) 
            eldiag = aleg*nb*nb
            elsqre = alpha*nb*nb+gamm*nb*nb*la
            call dcopy(nb*nb,w2(1+elsqre),1,w2(1+eldiag),1)
            !$OMP PARALLEL PRIVATE(el,delta_b,beta_b,beta,delta)&
            !$OMP& SHARED(bctr,bs,trick,nb,aleg,nbnb,modb)&
            !$OMP& DEFAULT(NONE)
            if(nbnb>0)then
              !$OMP DO
              do delta_b=1,nbnb,bs
                do beta_b=delta_b+bs,nbnb,bs
                  do delta=0,bctr
                    do beta=0,bctr
                      el=trick(beta+beta_b,delta_b+delta,aleg+1)
                         trick(beta+beta_b,delta_b+delta,aleg+1)=&
                        &trick(delta_b+delta,beta+beta_b,aleg+1)- &
                        &trick(beta+beta_b,delta_b+delta,aleg+1) 
                       trick(delta_b+delta,beta+beta_b,aleg+1)=&
                      &el - trick(delta_b+delta,beta+beta_b,aleg+1) 
                    enddo
                  enddo
                enddo
              enddo
              !$OMP END DO NOWAIT
            endif
            if(nbnb>0.and.modb)then
              !$OMP DO
              do delta_b=1,nbnb,bs
                do delta=0,bctr
                  do beta=nbnb+1,nb
                    el=trick(beta,delta_b+delta,aleg+1)
                       trick(beta,delta_b+delta,aleg+1)=&
                      &trick(delta_b+delta,beta,aleg+1)- &
                      &trick(beta,delta_b+delta,aleg+1) 
                     trick(delta_b+delta,beta,aleg+1)=&
                    &el -trick(delta_b+delta,beta,aleg+1)
                  enddo
                enddo
              enddo
              !$OMP END DO NOWAIT
            endif
            if(nbnb>0)then
              !$OMP DO
              do delta_b=1,nbnb,bs
                do delta=0,bctr
                  do beta=delta+1,bctr
                    el=trick(beta+delta_b,delta_b+delta,aleg+1)
                       trick(beta+delta_b,delta_b+delta,aleg+1)=&
                    &trick(delta_b+delta,beta+delta_b,aleg+1) - &
                      &trick(beta+delta_b,delta_b+delta,aleg+1)
                     trick(delta_b+delta,beta+delta_b,aleg+1)=&
                    &el - trick(delta_b+delta,beta+delta_b,aleg+1)
                  enddo
                   trick(delta_b+delta,delta+delta_b,aleg+1)=&
                  &trick(delta_b+delta,delta+delta_b,aleg+1) - &
                  &trick(delta+delta_b,delta_b+delta,aleg+1)
                enddo
              enddo
              !$OMP END DO NOWAIT
            endif
            !$OMP END PARALLEL 
            if(modb)then
              do delta=nbnb+1,nb
                do beta=delta+1,nb
                  el=trick(beta,delta,aleg+1)
                     trick(beta,delta,aleg+1)=&
                    &trick(delta,beta,aleg+1) - &
                  &trick(beta,delta,aleg+1)
                   trick(delta,beta,aleg+1)=&
                  &el- trick(delta,beta,aleg+1)
                enddo
                 trick(delta,delta,aleg+1)=&
                &trick(delta,delta,aleg+1) - &
                &trick(delta,delta,aleg+1)
              enddo
            endif
            aleg=aleg+1
          endif
        enddo
      enddo
      call array_reorder_3d(1.0E0_realk,w2,nb,nb,cagi,[2,3,1],0.0E0_realk,w0)
    endif
    nullify(trick)
  end subroutine get_I_plusminus_le
    

  !> \brief Reorder t to use symmetry in both occupied and virtual indices,
  !> thereby restricting the first virtual to be less equal to the second and
  !> make symmetric and antisymmetric combinations of these 
  !> \author Patrick Ettenhuber
  !> \date October 2012
  subroutine get_tpl_and_tmi(t2,nv,no,tpl,tmi)
    implicit none
    !> number of virtual and occupied orbitals
    integer, intent(in) :: nv, no
    !>input the full doubles ampitudes
    real(realk),intent(in) :: t2(nv,nv,no,no)
    !> output amplitudes with indices restricted c<d,i<j
    real(realk),pointer :: tpl(:)
    !> output amplitudes with indices restricted c>d,i<j
    real(realk),pointer :: tmi(:)
    integer :: dims(4)
    integer ::i,j,d, ilej,cged,crvd
    integer :: counter,counter2
    dims(1)=nv;dims(2)=nv;dims(3)=no;dims(4)=no
    !amplitudes are assumed to be ordered cdij, now construction of c<=d,i>=j
    if(dims(1)/=dims(2) .or. dims(3)/=dims(4))then
      call lsquit("ERROR(get_tpl):amplitudes are in the wrong order",DECinfo%output)
    endif
    !get the combined reduced virtual dimension
    crvd = dims(1) * (dims(1)+1) / 2
    !OMP PARALLEL PRIVATE(i,j,d,cged,ilej) SHARED(crvd,t2,tmi,tpl) DEFAULT(NONE)
    !OMP DO
    do j=1,dims(3)
      do i=1,j
        cged=1
        ilej = i + (j-1) * j / 2
        do d=1,dims(2)
          !copy the elements c>d into tpl and tmi
          call dcopy(dims(2)-d+1,t2(d,d,i,j),1,tpl(cged+(ilej-1)*crvd),1)
          call dcopy(dims(2)-d+1,t2(d,d,i,j),1,tmi(cged+(ilej-1)*crvd),1)

          !add and subtract the counterparts for tpl and tmi respectively
          call daxpy(dims(1)-d+1,-1.0E0_realk,t2(d,d,i,j),dims(1),tmi(cged+(ilej-1)*crvd),1)
          call daxpy(dims(1)-d+1,1.0E0_realk,t2(d,d,i,j),dims(1),tpl(cged+(ilej-1)*crvd),1)

          !take care of plus combination diagonal elements to be only half
          tpl(cged+(ilej-1)*crvd)=0.5*tpl(cged+(ilej-1)*crvd)
         
          !count in the triangular part of the matrix
          cged = cged+dims(2)-d+1
        enddo
      enddo
    enddo
    !OMP END DO
    !OMP END PARALLEL
  end subroutine get_tpl_and_tmi



  !> \brief calculate batch sizes automatically-->dirty but better than nothing
  !> \author Patrick Ettenhuber
  !> \date January 2012
  recursive subroutine get_max_batch_sizes(scheme,nb,nv,no,nba,nbg,&
  &minbsize,manual,iter,MemFree,first)
    implicit none
    integer, intent(inout) :: scheme
    integer, intent(in)    :: nb,nv,no
    integer :: iter
    integer, intent(inout) :: nba,nbg,minbsize
    real(realK),intent(in) :: MemFree
    real(realk)            :: mem_used,frac_of_total_mem,m
    logical,intent(in)     :: manual,first
    integer :: nnod,magic

    frac_of_total_mem=0.80E0_realk
    nba=minbsize
    nbg=minbsize
    nnod=1
    magic = DECinfo%MPIsplit/5
#ifdef VAR_LSMPI
    nnod=infpar%lg_nodtot
#endif
    !test for scheme with highest reqirements --> fastest
    mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,4,.false.)
    if(first)then
      if (mem_used>frac_of_total_mem*MemFree)then
#ifdef VAR_LSMPI
        !test for scheme with medium requirements
        mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,3,.false.)
        if (mem_used>frac_of_total_mem*MemFree)then
          !test for scheme with low requirements
          mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,2,.false.)
          if (mem_used>frac_of_total_mem*MemFree)then
            !test for scheme with low requirements
            mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,1,.false.)
#else
            !test for scheme with lowest reqirements
            mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,0,.false.)
#endif
            if (mem_used>frac_of_total_mem*MemFree)then
              write(DECinfo%output,*) "MINIMUM MEMORY REQUIREMENT IS NOT AVAILABLE"
              write(DECinfo%output,'("Fraction of free mem to be used:          ",f8.3," GB")')&
              &frac_of_total_mem*MemFree
              write(DECinfo%output,'("Memory required in memory saving scheme:  ",f8.3," GB")')mem_used
              mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,1,.false.)
              write(DECinfo%output,'("Memory required in intermediate scheme: ",f8.3," GB")')mem_used
              mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,2,.false.)
              write(DECinfo%output,'("Memory required in memory wasting scheme: ",f8.3," GB")')mem_used
              call lsquit("ERROR(CCSD): there is just not enough memory&
              &available",DECinfo%output)
#ifndef VAR_LSMPI
            else
              scheme=0
            endif
#else
            else
              scheme=1
            endif
          else
            scheme=2
          endif
        else
          scheme=3
        endif
#endif
      else
        scheme=4
      endif
    endif

    if(DECinfo%force_scheme)then
      scheme=DECinfo%en_mem
      print *,"!!FORCING CCSD!!"
      if(.not.DECinfo%solver_par)then
        if(scheme==3.or.scheme==2.or.scheme==1)then
          print *,"CHOSEN SCHEME DOES NOT WORK WITHOUT PARALLEL SOLVER, USE&
          & .CCSDsolver_par IN DALTON.INP"
          call lsquit("ERROR(get_ccsd_residual_integral_driven):invalid scheme",-1)
        endif
      endif
      if(scheme==4)then
        print *,"NON PDM-SCHEME WITH HIGH MEMORY REQUIREMENTS"
      elseif(scheme==3)then
        print *,"SCHEME WITH MEDIUM MEMORY REQUIREMENTS (PDM)"
      elseif(scheme==2)then
        print *,"SCHEME WITH LOW MEMORY REQUIREMENTS (PDM)"
      elseif(scheme==1)then
        print *,"SCHEME WITH LOW MEMORY REQUIREMENTS (PDM), HIGH SCALING LOW COMM"
      elseif(scheme==0)then
        print *,"NON PDM-SCHEME WITH LOW MEMORY REQUIREMENTS"
      endif
    endif

#ifndef VAR_LSMPI
      !scheme3 and  2 and 1 are pure mpi-scheme, this means that a redefinition in the case
      !of non mpi-builds is necessary
      if(scheme==3.or.scheme==2.or.scheme==1) scheme=0
#endif
    


    ! Attention this manual block should be only used for debugging, also you
    ! will have to ajust its properties to your system. The block is called
    ! with the keywordk ".manual_batchsizes" in DALTON.INP the next line
    ! then contains the alpha and gamma batch sizes separated by space
    if (manual) then
      ! KK and PE hacks -> only for debugging
      ! extended to mimic the behaviour of the mem estimation routine when memory is filled up
      if((DECinfo%ccsdGbatch==0).and.(DECinfo%ccsdAbatch==0)) then
        call get_max_batch_sizes(scheme,nb,nv,no,nba,nbg,minbsize,.false.,iter,MemFree,.false.)
      else
        nba=DECinfo%ccsdAbatch - iter *0
        nbg=DECinfo%ccsdGbatch - iter *0
      endif
      ! Use value given in input --> the zero can be adjusted to vary batch sizes during the iterations
      m = nbg-0*iter
      if (minbsize <  m) nbg = m
      if (minbsize >= m) nbg = minbsize
      m = nba-0*iter
      if (minbsize <  m) nba = m
      if (minbsize >= m) nba = minbsize

      if (nbg>=nb) nbg = nb
      if (nba>=nb) nba = nb

      mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,scheme,.false.)

      if (frac_of_total_mem*MemFree<mem_used) then
        print *, "ATTENTION your chosen batch sizes might be too large!!!"
      endif

    ! This block is routinely used in a calculation
    else
      !make batches larger until they do not fit anymore
      !determine gamma batch first
      do while ((frac_of_total_mem*MemFree>mem_used) .and. (nb>=nbg))
        nbg=nbg+1
        mem_used=get_min_mem_req(no,nv,nb,nba,nbg,3,scheme,.false.)
      enddo
      if (nbg>=nb)then
        nbg = nb
      else if (nbg<=minbsize)then
        nbg = minbsize
      else
        nbg=nbg-1
      endif
      !determine
      do while ((frac_of_total_mem*MemFree>mem_used) .and. (nb>=nba))
        nba=nba+1
        mem_used=get_min_mem_req(no,nv,nb,nba,nbg,3,scheme,.false.)
      enddo
      if (nba>=nb)then
         nba = nb
      else if (nba<=minbsize)then
         nba = minbsize
      else
         nba=nba-1
      endif
    endif
    !mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,first)
    mem_used=get_min_mem_req(no,nv,nb,nba,nbg,4,scheme,.false.)
    !if(first)write(DECinfo%output,*)"request batches with: N(gamma):",nbg," and N(alpha)",nba
    !if(first)write(DECinfo%output,*)"will use:",mem_used,"GB of ",MemFree, "GB available"
    !if(first)write(DECinfo%output,*)"------>",mem_used/MemFree * 100, "%"
    !divide the work if more nodes than jobs are available


    !if much more slaves than jobs are available, split the jobs to get at least
    !one for all the slaves
    if((nb/nba)*(nb/nbg)<magic*nnod.and.(nba>minbsize).and.nnod>1)then
      nba=(nb/(magic*nnod))
      if(nba<minbsize)nba=minbsize
    endif
    if((nb/nba)*(nb/nbg)<magic*nnod.and.(nba==minbsize).and.nnod>1)then
      do while((nb/nba)*(nb/nbg)<magic*nnod)
        nbg=nbg-1
        if(nbg<0)exit
      enddo
      if(nbg<minbsize)nbg=minbsize
    endif
  end subroutine get_max_batch_sizes



!> \brief calculate the memory requirement for the matrices in the ccsd routine
!> \author Patrick Ettenhuber
!> \date January 2012
  function get_min_mem_req(no,nv,nb,nba,nbg,choice,memintensive,print_stuff) result (memrq)
    implicit none
    integer, intent(in) :: no,nv,nb
    integer, intent(in) :: nba,nbg
    integer, intent(in) :: choice
    real(realk) :: memrq, memin, memout
    logical, intent(in) :: print_stuff
    integer, intent(in) :: memintensive
    integer :: nor, nvr, fe, ne , nnod, me, i
    integer :: d1(4),ntpm(4),tdim(4),mode,splt,ntiles,tsze
    integer(kind=ls_mpik) :: master
    integer :: l1,ml1,fai1,tl1
    integer :: l2,ml2,fai2,tl2
    integer :: l3,ml3,fai3,tl3
    integer :: l4,ml4,fai4,tl4
    integer :: nloctiles
    integer :: cd , e2
    nor = no*(no+1)/2
    nvr = nv*(nv+1)/2
    master = 0
    mode=4
    d1(1)=nv;d1(2)=nv;d1(3)=no;d1(4)=no
    nnod = 1
    me = 0
#ifdef VAR_LSMPI
    nnod = infpar%lg_nodtot
#endif
    l1   = (nv*no) / nnod
    l2   = (nv*nv) / nnod
    l3   = (nv*no*no) / nnod
    l4   = (nv*nv*no) / nnod
    ml1  = mod(nv*no,nnod)
    ml2  = mod(nv*nv,nnod)
    ml3  = mod(nv*no*no,nnod)
    ml4  = mod(nv*nv*no,nnod)
    fai1 = me * l1 + 1
    fai2 = me * l2 + 1
    fai3 = me * l3 + 1
    fai4 = me * l4 + 1
    tl1  = l1
    tl2  = l2
    tl3  = l3
    tl4  = l4

    if(ml1>0)then
      if(me<ml1)then
        fai1 = fai1 + me
        tl1  = l1 + 1
      else
        fai1 = fai1 + ml1
        tl1  = l1
      endif
    endif
    tl1 = tl1 * no * nv
    if(ml2>0)then
      if(me<ml2)then
        fai2 = fai2 + me
        tl2  = l2 + 1
      else
        fai2 = fai2 + ml2
        tl2  = l2
      endif
    endif
    tl2 = tl2 * no * no
    if(ml3>0)then
      if(me<ml3)then
        fai3 = fai3 + me
        tl3  = l3 + 1
      else
        fai3 = fai3 + ml3
        tl3  = l3
      endif
    endif
    tl3 = tl3 * nv
    if(ml4>0)then
      if(me<ml4)then
        fai4 = fai4 + me
        tl4  = l4 + 1
      else
        fai4 = fai4 + ml4
        tl4  = l4
      endif
    endif
    tl4 = tl4 * no
    !calculate minimum memory requirement
    ! u+3*integrals
    if(memintensive==4)then
      !govov + u 2 + Omega 2 +  H +G  
      memrq = 1.0E0_realk*(3*no*no*nv*nv+ nb*nv+nb*no)
      !gvoov gvvoo
      memrq=memrq+ 2.0E0_realk*nv*nv*no*no
      !allocation of matries ONLY used inside the loop
      !uigcj sio4
      memin  = 1.0E0_realk*(no*no*nv*nbg+no*no*nor)
      !tpl tmi
      memin  = memin + nor*nvr*2
      !w0
      memin = memin +&
      &nb*nb*nba*nbg
      !w1
      memin = memin +&
      &max(max(max(nb*nb*nba*nbg,nv*nv*no*nba),no*no*nv*nbg),no*no*nv*nba)
      !w2
      memin = memin +&
      &max(max(nb*nb*nba*nbg,nv*nv*no*no),nor*no*no)
      !w3
      memin = memin +&
      &max(max(max(max(max(max(max(nv*no*nba*nbg,no*no*nba*nbg),no*no*nv*nba),&
      &2*nor*nba*nbg),nor*nv*nba),nor*nv*nbg),no*nor*nba),no*nor*nbg)
      ! allocation of matrices ONLY used outside loop
      ! w1 + FO + w2 + w3
      memout = 1.0E0_realk*(max(nv*nv*no*no,nb*nb)+nb*nb+2*no*no*nv*nv)
      !memrq=memrq+max(memin,memout)
    elseif(memintensive==3)then
      !govov stays in pdm and is dense in second part
      ! u 2 + omega2 + H +G  + keep space for one update batch
      call get_int_dist_info(int(nv*nv*no*no,kind=8),fe,ne,master)
      memrq = 1.0E0_realk*(2*no*no*nv*nv+ nb*nv+nb*no) + ne
      !gvoov gvvoo (govov, allocd outside)
      memrq=memrq+ 2.0E0_realk*ne 
      !call array_default_batches([no,nv,no,nv],4,tdim,splt)
      !call array_get_ntpm([no,nv,no,nv],tdim,4,ntpm,ntiles)
      !allocation of matries ONLY used inside the loop
      !uigcj sio4
      memin  = 1.0E0_realk*(no*no*nv*nbg+no*no*nor)
      !tpl tmi
      memin  = memin + nor*nvr*2
      !w0
      memin = memin +&
      &nb*nb*nba*nbg
      !w1
      memin = memin +&
      &max(max(max(nb*nb*nba*nbg,nv*nv*no*nba),no*no*nv*nbg),no*no*nv*nba)
      !w2
      memin = memin +&
      &max(max(nb*nb*nba*nbg,nv*nv*no*no),nor*no*no)
      !w3
      memin = memin +&
      &max(max(max(max(max(max(max(nv*no*nba*nbg,no*no*nba*nbg),no*no*nv*nba),&
      &2*nor*nba*nbg),nor*nv*nba),nor*nv*nbg),no*nor*nba),no*nor*nbg)
      ! allocation of matrices ONLY used outside loop
      ! w1 + FO + w2 + w3 + govov
      memout = 1.0E0_realk*(max(nv*nv*no*no,nb*nb)+max(nb*nb,max(2*tl1,tl2)))
      !memrq=memrq+max(memin,memout)
    elseif(memintensive==2)then
      call array_default_batches(d1,mode,tdim,splt)
      call array_get_ntpm(d1,tdim,mode,ntpm,ntiles)
      nloctiles=ntiles/nnod
      if(mod(ntiles,nnod)>0)nloctiles=nloctiles+1
      tsze = 1
      do i = 1, mode
        tsze = tsze * tdim(i)
      enddo
      !govov stays in pdm and is dense in second part
      ! u 2 + H +G + space for one update tile 
      memrq = 1.0E0_realk*(tsze*nloctiles+ nb*nv+nb*no) + tsze
      !gvoov gvvoo
      memrq=memrq+ 2.0E0_realk*tsze*nloctiles
      !allocation of matries ONLY used inside the loop
      !uigcj sio4
      memin  = 1.0E0_realk*(no*no*nv*nbg+no*no*nor)
      !tpl tmi
      memin  = memin + nor*nvr*2
      !w0
      memin = memin +&
      &nb*nb*nba*nbg
      !w1
      memin = memin +&
      &max(max(max(nb*nb*nba*nbg,nv*nv*no*nba),no*no*nv*nbg),no*no*nv*nba)
      !w2
      memin = memin +&
      &max(max(nb*nb*nba*nbg,nv*nv*no*no),nor*no*no)
      !w3
      memin = memin +&
      &max(max(max(max(max(max(max(nv*no*nba*nbg,no*no*nba*nbg),no*no*nv*nba),&
      &2*nor*nba*nbg),nor*nv*nba),nor*nv*nbg),no*nor*nba),no*nor*nbg)
      ! allocation of matrices ONLY used outside loop
      ! w1 + FO + w2 + w3
      !in cd terms w2 and w3 have tl1, in b2 w2 has tl2
      cd = max(2*tl1,tl2)
      ! in e2 term w2 has max(tl2,tl3) and w3 has max(no2,nv2)
      e2 = max(tl3,tl4) + max(no*no,nv*nv)
      memout = 1.0E0_realk*(max(nv*nv*no*no,nb*nb)+max(nb*nb,max(cd,e2)))
      !memrq=memrq+max(memin,memout)
    elseif(memintensive==1)then
      print *,"ATTENTION, mem estimation not yet correct"
      !govov stays in pdm and is dense in second part
      ! u 2 + H +G  
      memrq = 1.0E0_realk*(no*no*nv*nv+ nb*nv+nb*no)
      !gvoov gvvoo
      call get_int_dist_info(int(nv*nv*no*no,kind=long),fe,ne,master)
      memrq=memrq+ 2.0E0_realk*ne
      !allocation of matries ONLY used inside the loop
      !uigcj sio4
      memin  = 1.0E0_realk*(no*no*nv*nbg+no*no*nor)
      !tpl tmi
      memin  = memin + nor*nvr*2
      !w0
      memin = memin +&
      &nb*nb*nba*nbg
      !w1
      memin = memin +&
      &max(max(max(nb*nb*nba*nbg,nv*nv*no*nba),no*no*nv*nbg),no*no*nv*nba)
      !w2
      memin = memin +&
      &max(max(nb*nb*nba*nbg,nv*nv*no*no),nor*no*no)
      !w3
      memin = memin +&
      &max(max(max(max(max(max(max(nv*no*nba*nbg,no*no*nba*nbg),no*no*nv*nba),&
      &2*nor*nba*nbg),nor*nv*nba),nor*nv*nbg),no*nor*nba),no*nor*nbg)
      ! allocation of matrices ONLY used outside loop
      ! w1 + FO + w2 + w3 + govov
      memout = 1.0E0_realk*(max(nv*nv*no*no,nb*nb)+nb*nb+3*no*no*nv*nv)
      !memrq=memrq+max(memin,memout)
    elseif(memintensive==0)then
      !govov + u 2 + H +G  
      memrq = 1.0E0_realk*(2*no*no*nv*nv+ nb*nv+nb*no)
      !allocation of matries ONLY used inside the loop
      !tgdij uigcj sio4 
      memin  = 2.0E0_realk*no*no*nv*nbg + no*no*nor
      !tpl tmi
      memin  = memin + nor*nvr*2
      !w0
      memin = memin +&
      &nb*nb*nba*nbg
      !w1
      memin = memin +&
      &max(max(nb*nb*nba*nbg,nv*nv*no*nba),no*no*nv*nbg)
      !w2
      memin = memin +&
      &max(max(max(nb*nb*nba*nbg,nv*nv*no*no),nor*no*no),no*no*nv*nba)
      !w3
      memin = memin +&
      &max(max(max(max(max(max(max(nv*no*nba*nbg,no*no*nba*nbg),no*no*nv*nba),&
      &2*nor*nba*nbg),nor*nv*nba),nor*nv*nbg),no*nor*nba),no*nor*nbg)
      ! allocation of matrices ONLY used outside loop
      ! w1 + FO 
      memout = 1.0E0_realk*max(nv*nv*no*no,nb*nb)+nb*nb
      !memrq=memrq+max(memin,memout)
    else
      call lsquit("ERROR(get_min_mem_req):requested memory scheme not known",-1)
    endif

    if(print_stuff) then
      write(DECinfo%output,*) "Memory requirements:"
      write(DECinfo%output,*) "Basic  :",(memrq *8.0E0_realk)/(1.024E3_realk**3)
      write(DECinfo%output,*) "Part B :",(memin *8.0E0_realk)/(1.024E3_realk**3)
      write(DECinfo%output,*) "Part C :",(memout*8.0E0_realk)/(1.024E3_realk**3)
    endif
    select case(choice)
      case(1)
        memrq = memrq
      case(2)
        memrq = memrq + memout
      case(3)
        memrq = memrq + memin
      case(4)
        memrq = memrq + max(memin,memout)
    end select
    memrq =((memrq*8.0E0_realk)/(1.024E3_realk**3))
  end function get_min_mem_req



  !> \brief Singles residual for CCSD model
  subroutine getSinglesResidualCCSD(omega1,u,gao,pqfock,qpfock, & 
              xocc,xvirt,yocc,yvirt,nocc,nvirt)

    implicit none
    type(array2), intent(inout) :: omega1
    type(array4), intent(inout) :: u,gao
    type(array2), intent(inout) :: pqfock,qpfock
    type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt
    type(array2) :: tmp
    type(array4) :: qqpq, pppq
    integer, intent(in) :: nocc,nvirt
    integer :: a,i,k,l,c,d
    real(realk) :: starttime,endtime,aStart,aEnd,bStart,bEnd,cStart,cEnd, &
              dStart,dEnd

   
    aStart=0.0E0_realk; aEnd=0.0E0_realk
    bStart=0.0E0_realk; bEnd=0.0E0_realk
    cStart=0.0E0_realk; cEnd=0.0E0_realk
    dStart=0.0E0_realk; dEnd=0.0E0_realk
    starttime=0.0E0_realk; endtime=0.0E0_realk

    call cpu_time(starttime)

#ifdef SINGLES_EXTRA_SIMPLE    

    qqpq = get_gmo_simple(gao,xvirt,yvirt,xocc,yvirt)

    do i=1,nocc
      do a=1,nvirt

        do c=1,nvirt
        do k=1,nocc
        do d=1,nvirt
        omega1%val(a,i) = omega1%val(a,i) + & 
                u%val(c,k,d,i)*qqpq%val(a,d,k,c)
        end do
        end do
        end do

      end do
    end do
    print *,'A1 ',omega1*omega1

    call array4_free(qqpq)

    pppq = get_gmo_simple(gao,xocc,yocc,xocc,yvirt)

    do i=1,nocc
      do a=1,nvirt

        do c=1,nvirt
        do k=1,nocc
        do l=1,nocc

        omega1%val(a,i) = omega1%val(a,i) - &
            u%val(a,k,c,l) * pppq%val(k,i,l,c)

        end do
        end do
        end do

      end do
    end do
    print *,'B1 ',omega1*omega1

    call array4_free(pppq)

    do i=1,nocc
      do a=1,nvirt
      
        do c=1,nvirt
          do k=1,nocc
            omega1%val(a,i) = omega1%val(a,i) + &
                u%val(a,i,c,k) * pqfock%val(k,c)
          end do
        end do  

      end do
    end do
    print *,'C1 ',omega1*omega1

    omega1%val = omega1%val + qpfock%val
    print *,'D1 ',omega1*omega1

#else

    tmp = array2_init([nvirt,nocc])

    ! A1
    call cpu_time(aStart)
    qqpq = get_gmo_simple(gao,xocc,yvirt,yvirt,xvirt)
    call array4_reorder(qqpq,[2,1,3,4]) ! qqpq[ic,ab] -> qqpq[ci,ab]
    call array4_contract3(qqpq,u,tmp)
    call array4_free(qqpq)
    call array2_add_to(omega1,1.0E0_realk,tmp)
    call cpu_time(aEnd)
    if(DECinfo%cc_driver_debug) then
      print *,'A1 done, norm : ',omega1*omega1
    end if

    ! B1
    call cpu_time(bStart)
    pppq = get_gmo_simple(gao,xocc,yocc,xocc,yvirt)
    call array4_reorder(u,[4,3,2,1])  ! u[ai,bj] -> u[jb,ia]
    call array4_reorder(pppq,[3,4,1,2])
    call array2_zero(tmp)
    call array4_contract3(u,pppq,tmp)
    call array2_add_to(omega1,-1.0E0_realk,tmp)
    call array4_free(pppq)
    call cpu_time(bEnd)
    if(DECinfo%cc_driver_debug) then
      print *,'B1 done, norm: ',omega1*omega1
    end if

    ! C1
    call cpu_time(cStart)
    call array2_zero(tmp)
    call array4_reorder(u,[1,2,4,3]) ! u[jb,ia] -> u[jb,ai] 
    call array4_contract_array2(u,pqfock,tmp)
    call array4_reorder(u,[3,4,2,1])
    call array2_add_to(omega1,1.0E0_realk,tmp)
    call cpu_time(cEnd)
    if(DECinfo%cc_driver_debug) then
      print *,'C1 done, norm: ',omega1*omega1
    end if

    ! D1
    call cpu_time(dStart)
    call array2_add_to(omega1,1.0E0_realk,qpfock)
    call cpu_time(dEnd)
    if(DECinfo%cc_driver_debug) then
      print *,'D1 done, norm: ',omega1*omega1
    end if

    call array2_free(tmp)
    call cpu_time(endtime)
    if(DECinfo%PL>1) then
      write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD A1 : ',aEnd-aStart,' s'
      write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD B1 : ',bEnd-bStart,' s'
      write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD C1 : ',cEnd-cStart,' s'
      write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD D1 : ',dEnd-dStart,' s'
      write(DECinfo%output,'(a,f16.3,a)')  &
         ' time :: CCSD singles : ',endtime-starttime,' s'
    end if

#endif

    return
  end subroutine getSinglesResidualCCSD
  
  !> \brief Precondition singles 
  function precondition_singles_newarr(omega1,ppfock,qqfock) result(prec)

    implicit none
    type(array), intent(in) :: omega1,ppfock,qqfock
    type(array) :: prec
    integer, dimension(2) :: dims
    integer :: a,i
    if(omega1%mode/=2.or.ppfock%mode/=2.or.qqfock%mode/=2)then
      call lsquit("ERROR(precondition_singles_newarr):wrong number of modes&
      & for this operation",DECinfo%output)
    endif

    dims = omega1%dims
    prec = array_init(dims,2)

    do a=1,dims(1)
      do i=1,dims(2)
      
        prec%elm2(a,i) = omega1%elm2(a,i)/( ppfock%elm2(i,i) - qqfock%elm2(a,a) ) 

      end do
    end do
  end function precondition_singles_newarr
  function precondition_singles_oldarr(omega1,ppfock,qqfock) result(prec)
    implicit none
    type(array2), intent(in) :: omega1,ppfock,qqfock
    type(array2) :: prec
    integer, dimension(2) :: dims
    integer :: a,i

    dims = omega1%dims
    prec = array2_init(dims)

    do a=1,dims(1)
      do i=1,dims(2)
      
        prec%val(a,i) = omega1%val(a,i)/( ppfock%val(i,i) - qqfock%val(a,a) ) 

      end do
    end do
  end function precondition_singles_oldarr


  !> \brief T1 transformations
  subroutine getT1transformation(t1,xocc,xvirt,yocc,yvirt, &
       & ypo,ypv,yho,yhv)

    implicit none
    type(array2), intent(inout) :: t1
    type(array2), intent(inout) :: ypo,ypv,yho,yhv
    type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt

    ! Occupied X and Y matrices
    call getT1transformation_occ(t1,xocc,yocc,ypo,yho,yhv)

    ! Virtual X and Y matrices
    call getT1transformation_virt(t1,xvirt,yvirt,ypo,ypv,yhv)


  end subroutine getT1transformation


  !> \brief T1 transformations for occupied X and Y matrices
  subroutine getT1transformation_occ(t1,xocc,yocc,ypo,yho,yhv)

    implicit none
    type(array2), intent(inout) :: t1
    type(array2), intent(inout) :: ypo,yho,yhv
    type(array2), intent(inout) :: xocc,yocc

    ! xocc
    call array2_zero(xocc)
    call array2_copy(xocc,ypo)

    ! yocc
    call array2_zero(yocc)
    call array2_copy(yocc,yho)
    call array2_matmul(yhv,t1,yocc,'n','n',1.0E0_realk,1.0E0_realk)


  end subroutine getT1transformation_occ



  !> \brief T1 transformations for virtual X and Y matrices
  subroutine getT1transformation_virt(t1,xvirt,yvirt,ypo,ypv,yhv)

    implicit none
    type(array2), intent(inout) :: t1
    type(array2), intent(inout) :: ypo,ypv,yhv
    type(array2), intent(inout) :: xvirt,yvirt

    ! xvirt
    call array2_zero(xvirt)
    call array2_copy(xvirt,ypv)
    call array2_matmul(ypo,t1,xvirt,'n','t',-1.0E0_realk,1.0E0_realk)

    ! yvirt
    call array2_zero(yvirt)
    call array2_copy(yvirt,yhv)

  end subroutine getT1transformation_virt



  !> \brief Get inactive fock matrix in ao basis
  function getInactiveFock(h1,gao,xocc,yocc,nocc,nbasis) result(ifock)

    implicit none
    type(array4) :: tmp1,tmp2
    type(array2) :: ifock
    type(array2), intent(in) :: h1
    type(array2), intent(inout) :: xocc,yocc
    type(array4), intent(inout) :: gao
    integer, intent(in) :: nocc,nbasis
    integer :: mu,nu,i

    ifock = array2_init([nbasis,nbasis])

    ! transform to (oo,mu nu)
    tmp1 = array4_init([nocc,nbasis,nbasis,nbasis])
    call array4_read(gao)
    call array4_contract1(gao,yocc,tmp1,.true.)
    call array4_dealloc(gao)
    call array4_reorder(tmp1,[2,1,3,4]) ! (sig o,mu nu)

    tmp2 = array4_init([nocc,nocc,nbasis,nbasis])
    call array4_contract1(tmp1,xocc,tmp2,.true.) ! (oo,mu nu)

    ! sum over 'o'
    do nu=1,nbasis
       do mu=1,nbasis
          do i=1,nocc
             ifock%val(mu,nu) = ifock%val(mu,nu) + 2.0E0_realk*tmp2%val(i,i,mu,nu)
          end do
       end do
    end do

    call array4_free(tmp2)
    print *,'ifock coul ',ifock*ifock

    ! transformation to (mu o,o nu)
    tmp2 = array4_init([nocc,nocc,nbasis,nbasis])
    call array4_reorder(tmp1,[3,2,1,4])
    call array4_contract1(tmp1,xocc,tmp2,.true.)
    call array4_reorder(tmp2,[3,2,1,4])
    call array4_free(tmp1)

    do nu=1,nbasis
       do mu=1,nbasis
          do i=1,nocc
             ifock%val(mu,nu) = ifock%val(mu,nu) + tmp2%val(mu,i,i,nu)
          end do
       end do
    end do

    call array2_add_to(ifock,1.0E0_realk,h1)

    return
  end function getInactiveFock

  !> Get inactive Fock matrix (simple version)
  function getInactiveFock_simple(h1,gao,xocc,yocc,nocc,nbasis) result(ifock)
    implicit none
    type(array4), intent(inout) :: gao
    type(array2) :: ifock
    type(array2), intent(inout) :: xocc,yocc,h1
    integer, intent(in) :: nbasis,nocc
    type(array2) :: unit
    type(array4) :: pqoo, pooq
    integer :: i,p,q
    real(realk) :: tmp
    ifock = array2_init([nbasis,nbasis])

    ! coul
    pqoo = get_oopq(gao,xocc,yocc)
    do q=1,nbasis
       do p=1,nbasis
          tmp=0E0_realk
          do i=1,nocc
             tmp = tmp + pqoo%val(i,i,p,q)
          end do
          ifock%val(p,q) = ifock%val(p,q) + 2.0E0_realk * tmp
       end do
    end do
    call array4_free(pqoo)

    ! ex
    pooq = get_exchange_as_oopq(gao,yocc,xocc)
    do q=1,nbasis
       do p=1,nbasis
          tmp = 0E0_realk
          do i=1,nocc
             tmp = tmp + pooq%val(i,i,p,q)
          end do
          ifock%val(p,q) = ifock%val(p,q) - tmp
       end do
    end do
    call array4_free(pooq)

    call array2_add_to(ifock,1.0E0_realk,h1)

    return
  end function getInactiveFock_simple

  !> \brief Get fock correction
  function getFockCorrection(fock_full,ifock) result(fock_correction)

    implicit none
    type(array2) :: fock_correction
    type(array2), intent(in) :: fock_full,ifock

    fock_correction = array2_add(1.0E0_realk,fock_full,-1.0E0_realk,ifock)

    return
  end function getFockCorrection

  !> \brief Simple Fock from RI integrals
  function getInactiveFockFromRI(l_ao,xocc,yocc,h1) result(this)

    implicit none
    type(array2) :: this
    type(ri), intent(in) :: l_ao
    type(array2), intent(in) :: h1,xocc,yocc
    type(ri) :: IJ,alphaI,Ibeta
    integer :: nocc,naux,nbas,l,i
    real(realk) :: trace
    real(realk), pointer :: tmpfock(:,:)

    nbas = xocc%dims(1)
    nocc = xocc%dims(2)
    naux = l_ao%dims(3)

    IJ = ri_init([nocc,nocc,naux])
    alphaI = ri_init([nbas,nocc,naux])
    Ibeta = ri_init([nocc,nbas,naux])

    ! transform
    do l=1,naux
       IJ%val(:,:,l) = matmul(matmul(transpose(xocc%val),l_ao%val(:,:,l)), &
            yocc%val)
       alphaI%val(:,:,l) = matmul(l_ao%val(:,:,l),yocc%val)
       Ibeta%val(:,:,l) = matmul(transpose(xocc%val),l_ao%val(:,:,l))
    end do

    call mem_alloc(tmpfock,nbas,nbas)
    tmpfock = 0.0E0_realk

    do l=1,naux

       ! 2g_mu_nu_i_i
       trace=0E0_realk
       do i=1,nocc
          trace=trace+IJ%val(i,i,l)
       end do
       tmpfock=tmpfock+2E0_realk*trace*l_ao%val(:,:,l)

       ! -g_mu_i_i_nu
       tmpfock=tmpfock-matmul(alphaI%val(:,:,l),Ibeta%val(:,:,l))

    end do
    tmpfock=tmpfock+h1%val

    this = array2_init([nbas,nbas],tmpfock)
    call mem_dealloc(tmpfock)

    ! free
    call ri_free(IJ)
    call ri_free(alphaI)
    call ri_free(Ibeta)

    return
  end function getInactiveFockFromRI

  !> \brief Get T1 transformed Fock matrices
  subroutine getFockMatrices(ifock,xocc,xvirt,yocc,yvirt, &
       ppfock,qqfock,pqfock,qpfock,nbasis,nocc,nvirt)
    implicit none
    type(array2), intent(inout) :: ifock,xocc,xvirt,yocc,yvirt, &
         ppfock,qqfock,pqfock,qpfock
    integer, intent(in) :: nbasis,nocc,nvirt
    type(array2) :: p_tmp, q_tmp

    p_tmp = array2_init([nocc,nbasis])
    q_tmp = array2_init([nvirt,nbasis])

    call array2_transpose(xocc)
    call array2_transpose(xvirt)

    p_tmp%val = matmul(xocc%val,ifock%val)
    q_tmp%val = matmul(xvirt%val,ifock%val)

    ppfock%val = matmul(p_tmp%val,yocc%val)
    pqfock%val = matmul(p_tmp%val,yvirt%val)

    qpfock%val = matmul(q_tmp%val,yocc%val)
    qqfock%val = matmul(q_tmp%val,yvirt%val)

    call array2_free(p_tmp)
    call array2_free(q_tmp)

    call array2_transpose(xocc)
    call array2_transpose(xvirt)

    return
  end subroutine getFockMatrices


  !> \brief Get T1 transformed Fock matrix in the AO basis using
  !> the full molecular T1 amplitudes and the fullmolecule structure.
  subroutine fullmolecular_get_AOt1Fock(mylsitem,MyMolecule,t1,fockt1)
    implicit none
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Singles amplitudes for full molecular system
    !> Although this is intent(inout), it is unchanged at output
    type(array2),intent(inout) :: t1
    !> T1-transformed Fock matrix in the AO basis (also initialized here)
    type(array2),intent(inout) :: fockt1
    type(array2) :: ypo,yho,ypv
    integer :: nbasis,nocc,nvirt
    integer,dimension(2) :: bo,bv

    ! Dimensions
    nbasis = MyMolecule%nbasis
    nocc = MyMolecule%numocc
    nvirt = MyMolecule%numvirt
    bo(1) = nbasis
    bo(2) = nocc
    bv(1) = nbasis
    bv(2) = nvirt

    ! Initialize stuff in array2 format
    ypo = array2_init(bo,MyMolecule%ypo)
    yho = array2_init(bo,MyMolecule%ypo)
    ypv = array2_init(bv,MyMolecule%ypv)

    ! Get T1-transformed Fock matrix in AO basis
    call Get_AOt1Fock(mylsitem,t1,fockt1,nocc,nvirt,nbasis,ypo,yho,ypv)

    ! Free stuff
    call array2_free(ypo)
    call array2_free(yho)
    call array2_free(ypv)


  end subroutine fullmolecular_get_AOt1Fock



  !> \brief Get T1 transformed Fock matrix in the AO basis using
  !> the fragment T1 amplitudes (dimension: occ AOS, virt AOS) and fragment MOs.
  subroutine fragment_get_AOt1Fock(MyFragment,fockt1)
    implicit none
    !> Fragment info
    type(ccatom), intent(inout) :: MyFragment
    !> T1-transformed Fock matrix in the AO basis (also initialized here)
    type(array2),intent(inout) :: fockt1
    type(array2) :: t1
    type(array2) :: ypo,yho,ypv
    integer :: nbasis,nocc,nvirt
    integer,dimension(2) :: bo,bv,vo

    ! Dimensions
    ! **********
    nbasis = MyFragment%number_basis
    nocc = MyFragment%noccAOS
    nvirt = MyFragment%nunoccAOS
    bo(1) = nbasis
    bo(2) = nocc
    bv(1) = nbasis
    bv(2) = nvirt
    vo(1) = nvirt
    vo(2) = nocc

    ! Initialize stuff in array2 format
    ! *********************************
    ypo = array2_init(bo,MyFragment%ypo)
    yho = array2_init(bo,MyFragment%ypo)   ! particle and hole coefficients are identical
    ypv = array2_init(bv,MyFragment%ypv)

    ! Set t1 equal to the t1 associated with the fragment
    ! ***************************************************

    ! Sanity check
    if(.not. MyFragment%t1_stored) then
       call lsquit('fragment_get_AOt1Fock: &
            & t1 amplitudes for fragment are not stored!',DECinfo%output)
    end if
    t1 = array2_init(vo,MyFragment%t1)


    ! Get T1-transformed Fock matrix in AO basis
    fockt1=array2_init([nbasis,nbasis])
    call Get_AOt1Fock(MyFragment%mylsitem,t1,fockt1,nocc,nvirt,nbasis,ypo,yho,ypv)

    ! Free stuff
    call array2_free(ypo)
    call array2_free(yho)
    call array2_free(ypv)


  end subroutine fragment_get_AOt1Fock


  subroutine Get_AOt1Fock_arraywrapper(mylsitem,t1,fockt1,nocc,nvirt,nbasis,ypo,yho,ypv)
    implicit none
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Singles amplitudes for fragment or full molecule
    !> Although this is intent(inout), it is unchanged at output
    type(array2),intent(inout) :: t1
    !> T1-transformed Fock matrix in the AO basis (also initialized here)
    type(array),intent(inout) :: fockt1
    !> Number of occupied orbitals (fragment AOS or full molecule)
    integer,intent(in) :: nocc
    !> Number of virtual orbitals (fragment AOS or full molecule)
    integer,intent(in) :: nvirt
    !> Number of basis functions (atomic fragment extent or full molecule)
    integer,intent(in) :: nbasis
    !> Occupied MOs (particle)
    type(array),intent(inout) :: ypo
    !> Occupied MOs (hole - currently particle=hole always)
    type(array),intent(inout) :: yho
    !> Virtual MOs (hole)
    type(array),intent(inout) :: ypv

    type(array2) :: fockt1_a2
    type(array2) :: ypo_a2
    type(array2) :: yho_a2
    type(array2) :: ypv_a2

    fockt1_a2%dims=fockt1%dims
    ypo_a2%dims=ypo%dims
    yho_a2%dims=yho%dims
    ypv_a2%dims=ypv%dims

    call ass_D1to2(fockt1%elm1,fockt1_a2%val,fockt1%dims)
    call ass_D1to2(ypo%elm1,ypo_a2%val,ypo%dims)
    call ass_D1to2(yho%elm1,yho_a2%val,yho%dims)
    call ass_D1to2(ypv%elm1,ypv_a2%val,ypv%dims)

    call Get_AOt1Fock(mylsitem,t1,fockt1_a2,nocc,nvirt,nbasis,ypo_a2,yho_a2,ypv_a2)


    nullify(fockt1_a2%val)
    nullify(ypo_a2%val)
    nullify(yho_a2%val)
    nullify(ypv_a2%val)

    !call print_norm(fockt1)
  end subroutine Get_AOt1Fock_arraywrapper

  !> \brief Get T1 transformed Fock matrix in the AO basis using
  !> either fragment or full molecular T1 amplitudes.
  subroutine Get_AOt1Fock_oa(mylsitem,t1,fockt1,nocc,nvirt,nbasis,ypo,yho,ypv)
    implicit none
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Singles amplitudes for fragment or full molecule
    !> Although this is intent(inout), it is unchanged at output
    type(array2),intent(inout) :: t1
    !> T1-transformed Fock matrix in the AO basis (also initialized here)
    type(array2),intent(inout) :: fockt1
    !> Number of occupied orbitals (fragment AOS or full molecule)
    integer,intent(in) :: nocc
    !> Number of virtual orbitals (fragment AOS or full molecule)
    integer,intent(in) :: nvirt
    !> Number of basis functions (atomic fragment extent or full molecule)
    integer,intent(in) :: nbasis
    !> Occupied MOs (particle)
    type(array2),intent(inout) :: ypo
    !> Occupied MOs (hole - currently particle=hole always)
    type(array2),intent(inout) :: yho
    !> Virtual MOs (hole)
    type(array2),intent(inout) :: ypv
    type(array2) :: xocc,yocc
    type(matrix) :: fockmat,dens,h1
    integer :: i,idx1,idx2
    integer,dimension(2) :: bo,bb

    ! Dimensions
    bo(1) = nbasis
    bo(2) = nocc
    bb(1) = nbasis
    bb(2) = nbasis


    ! Get T1-modified occupied orbitals xocc and yocc
    xocc=array2_init(bo)
    yocc=array2_init(bo)
    call getT1transformation_occ(t1,xocc,yocc,ypo,yho,ypv)

    ! Get effective T1 density: Yocc Xocc^T
    call mat_init(dens,nbasis,nbasis)
    call dgemm('n','t',nbasis,nbasis,nocc,1.0E0_realk,yocc%val,&
         & nbasis,xocc%val,nbasis,0.0E0_realk,dens%elms,nbasis)
    call array2_free(xocc)
    call array2_free(yocc)

    ! Get two-electron part of T1-transformed Fock matrix
    call mat_init(fockmat,nbasis,nbasis)
    call mat_zero(fockmat)
    call dec_fock_transformation(fockmat,dens,MyLsitem,.false.)
    call mat_free(dens)

    ! Get one-electron contribution to Fock matrix
    call mat_init(h1,nbasis,nbasis)
    call mat_zero(h1)
    call ii_get_h1(DECinfo%output,DECinfo%output,mylsitem%setting,h1)

    ! Add one- and two-electron contributions of Fock matrix
    call mat_daxpy(1.0E0_realk,h1,fockmat)
    call mat_free(h1)

    ! Convert to output array2 format
    !fockt1 = array2_init(bb)
    do i=1,nbasis
       idx1=(i-1)*nbasis+1
       idx2=i*nbasis
       fockt1%val(1:nbasis,i) = fockmat%elms(idx1:idx2)
    end do
    call mat_free(fockmat)


    !call print_norm(fockt1)
  end subroutine Get_AOt1Fock_oa


  subroutine get_fock_matrix_for_dec_arraywrapper(nbasis,dens,mylsitem,fock_array,add_one_electron)
    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Density matrix for fragment (or full molecule), in any case it equals Cocc Cocc^T
    real(realk),intent(in) :: dens(nbasis,nbasis)
    !> Ls item information for fragment (or full molecule)
    type(lsitem), intent(inout) :: mylsitem
    !> Fock matrix in array2 form
    type(array), intent(inout) :: fock_array
    !> Add one electron matrix to Fock matrix
    logical,intent(in) :: add_one_electron
    integer :: nocc,i,idx1,idx2
    type(matrix) :: D,fock,h1

    ! Density matrix in matrix form
    call mat_init(D,nbasis,nbasis)
    call mat_set_from_full(dens,1.0E0_realk,D)

    ! Get Fock matrix for the fragment (or full molecule) density matrix
    call mat_init(fock,nbasis,nbasis)
    call mat_zero(fock)
    call dec_fock_transformation(fock,D,MyLsitem,.true.)
    call mat_free(D)


    if(add_one_electron) then
       ! Get one-electron contribution to Fock matrix
       call mat_init(h1,nbasis,nbasis)
       call mat_zero(h1)
       call ii_get_h1(DECinfo%output,DECinfo%output,mylsitem%setting,h1)
       
       ! Add one- and two-electron contributions of Fock matrix
       call mat_daxpy(1.0E0_realk,h1,fock)
       call mat_free(h1)
    end if

    ! Convert to output array2 format
    call dcopy(nbasis*nbasis,fock%elms,1,fock_array%elm1,1)
    call mat_free(fock)
  end subroutine get_fock_matrix_for_dec_arraywrapper

  !> \brief Get T1 transformed Fock matrix in AO basis for fragment, i.e.
  !> using a density matrix constructed using ONLY the occupied AOS
  !> orbitals in the fragment and the singles amplitudes for the fragment.
  !> (Also works for full molecule but not in which case the standard HF fock matrix is returned)
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine get_fock_matrix_for_dec_oa(nbasis,dens,mylsitem,fock_array2,add_one_electron)
    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Density matrix for fragment (or full molecule), in any case it equals Cocc Cocc^T
    real(realk),intent(in) :: dens(nbasis,nbasis)
    !> Ls item information for fragment (or full molecule)
    type(lsitem), intent(inout) :: mylsitem
    !> Fock matrix in array2 form
    type(array2), intent(inout) :: fock_array2
    !> Add one electron matrix to Fock matrix
    logical,intent(in) :: add_one_electron
    integer :: nocc,i,idx1,idx2
    type(matrix) :: D,fock,h1


    ! Density matrix in matrix form
    call mat_init(D,nbasis,nbasis)
    call mat_set_from_full(dens,1.0E0_realk,D)

    ! Get Fock matrix for the fragment (or full molecule) density matrix
    call mat_init(fock,nbasis,nbasis)
    call mat_zero(fock)
    call dec_fock_transformation(fock,D,MyLsitem,.true.)
    call mat_free(D)


    if(add_one_electron) then
       ! Get one-electron contribution to Fock matrix
       call mat_init(h1,nbasis,nbasis)
       call mat_zero(h1)
       call ii_get_h1(DECinfo%output,DECinfo%output,mylsitem%setting,h1)
       
       ! Add one- and two-electron contributions of Fock matrix
       call mat_daxpy(1.0E0_realk,h1,fock)
       call mat_free(h1)
    end if

    ! Convert to output array2 format
    do i=1,nbasis
       idx1=(i-1)*nbasis+1
       idx2=i*nbasis
       fock_array2%val(1:nbasis,i) = fock%elms(idx1:idx2)
    end do
    call mat_free(fock)

  end subroutine get_fock_matrix_for_dec_oa



end module ccsd_module


!> \brief slave function for data preparation
!> \author Patrick Ettenhuber
!> \date March 2012
#ifdef VAR_LSMPI
subroutine ccsd_data_preparation()
    use precision
    !use tensor_def_module
    use dec_typedef_module
    use typedeftype,only:lsitem,array
    use infpar_module
    use lsmpi_type, only:ls_mpibcast,ls_mpibcast_chunks
    use daltoninfo, only:ls_free
    use memory_handling, only: mem_alloc, mem_dealloc
    use tensor_interface_module, only: array_init,array_free,&
        &memory_allocate_array_dense,memory_deallocate_array_dense,memory_deallocate_window
    ! DEC DEPENDENCIES (within deccc directory) 
    ! *****************************************
    use decmpi_module, only: mpi_communicate_ccsd_calcdata
    use array2_simple_operations,only: array2_free,array2_init
    use array4_simple_operations,only: array4_free,array4_init
    use ccsd_module, only:get_ccsd_residual_integral_driven
    implicit none
    !type(mp2_batch_construction) :: bat
    type(array2) :: deltafock,ppfock,qqfock,pqfock,qpfock,omega1,fock
    type(array2) :: xocc,xvirt,yocc,yvirt
    type(array)  :: om2,t2,govov
    type(lsitem) :: MyLsItem
    logical :: solver_par
    integer :: nbas,nocc,nvirt,scheme
    integer(kind=long) :: nelms
    integer      :: iter,k,n4,i
    real(realk),pointer  :: xodata(:),xvdata(:),yodata(:),yvdata(:),&
                           & df(:),f(:),ppf(:),qqf(:),pqf(:),qpf(:),om1(:),t2data(:),&
                           & t2d(:,:,:,:),xod(:,:),xvd(:,:),yod(:,:),yvd(:,:)



    !note that for the slave all allocatable arguments are just dummy indices
    !the allocation and broadcasting happens in here
    call mpi_communicate_ccsd_calcdata(om2,t2,govov,xodata,xvdata,yodata,yvdata,&
    &MyLsItem,nbas,nvirt,nocc,iter,scheme,solver_par)
   
    DECinfo%solver_par=solver_par
    if(solver_par)then
      call memory_allocate_array_dense(t2)
      if(scheme==0..or.scheme==4)then
        call memory_allocate_array_dense(govov)
      endif
    else
      !for not solver_par there is only the option 0 or 4
      t2   =array_init([nvirt,nvirt,nocc,nocc],4)
      govov=array_init([nocc,nvirt,nocc,nvirt],4)
      om2  =array_init([nvirt,nvirt,nocc,nocc],4)
    endif
    
    ! Quantities, that need to be defined and setset
    ! ********************************************************
     

    !split messages in 2GB parts, compare to counterpart in
    !mpi_communicate_ccsd_calcdate
    k=250000000

    nelms = nbas*nocc
    call mem_alloc(yodata,nelms)
    call mem_alloc(xodata,nelms)
    call ls_mpibcast_chunks(xodata,nelms,infpar%master,infpar%lg_comm,k)
    call ls_mpibcast_chunks(yodata,nelms,infpar%master,infpar%lg_comm,k)

    nelms = nbas*nvirt
    call mem_alloc(xvdata,nelms)
    call mem_alloc(yvdata,nelms)
    call ls_mpibcast_chunks(xvdata,nelms,infpar%master,infpar%lg_comm,k)
    call ls_mpibcast_chunks(yvdata,nelms,infpar%master,infpar%lg_comm,k)
    

    nelms = int(nvirt*nvirt*nocc*nocc,kind=8)
    call ls_mpibcast_chunks(t2%elm1,nelms,infpar%master,infpar%lg_comm,k)
    if(iter/=1.and.(scheme==0.or.scheme==4))then
      call ls_mpibcast_chunks(govov%elm1,nelms,infpar%master,infpar%lg_comm,k)
    endif


    ! Quantities, that need to be defined but not set
    ! ********************************************************
    nelms=nbas*nbas
    call mem_alloc(df,nelms)
    call mem_alloc(f,nelms)
    
    nelms=nocc*nocc
    call mem_alloc(ppf,nelms)

    nelms=nvirt*nocc
    call mem_alloc(pqf,nelms)
    call mem_alloc(qpf,nelms)
    call mem_alloc(om1,nelms)
  
    nelms=nvirt*nvirt
    call mem_alloc(qqf,nelms)

    ! Calculate contribution to integrals/amplitudes for slave
    ! ********************************************************
    call get_ccsd_residual_integral_driven(df,om2,t2,f,govov,nocc,nvirt,&
       ppf,qqf,pqf,qpf,xodata,xvdata,yodata,yvdata,nbas,MyLsItem,om1,iter)

    ! FREE EVERYTHING
    ! ***************
    if(solver_par)then
      call memory_deallocate_array_dense(om2)
      call memory_deallocate_array_dense(t2)
      if(scheme==0..or.scheme==4)then
        call memory_deallocate_array_dense(govov)
      endif
    else
      call array_free(om2)
      call array_free(t2)
      call array_free(govov)
    endif
    call mem_dealloc(df)
    call mem_dealloc(f)
    call mem_dealloc(ppf)
    call mem_dealloc(pqf)
    call mem_dealloc(qpf)
    call mem_dealloc(qqf)
    call mem_dealloc(om1)
    call mem_dealloc(xodata)
    call mem_dealloc(yodata)
    call mem_dealloc(yvdata)
    call mem_dealloc(xvdata)
    call ls_free(MyLsItem)
end subroutine ccsd_data_preparation

subroutine calculate_E2_and_permute_slave()
    use precision
    use dec_typedef_module
    use typedeftype,only:lsitem,array
    use infpar_module
    use lsmpi_type, only:ls_mpibcast
    use daltoninfo, only:ls_free
    use memory_handling, only: mem_alloc, mem_dealloc
    ! DEC DEPENDENCIES (within deccc directory) 
    ! *****************************************
    use decmpi_module, only: share_E2_with_slaves
    use ccsd_module, only:calculate_E2_and_permute
    implicit none
    real(realk),pointer :: ppf(:),qqf(:),xo(:),yv(:),Gbi(:),Had(:)
    integer :: no,nv,nb,s
    type(array) :: t2,omega2
    real(realk),pointer :: w1(:)
    
    call share_E2_with_slaves(ppf,qqf,t2,xo,yv,Gbi,Had,no,nv,nb,omega2,s)
    call mem_alloc(ppf,no*no)
    call mem_alloc(qqf,nv*nv)
    call ls_mpibcast(ppf,no*no,infpar%master,infpar%lg_comm)
    call ls_mpibcast(qqf,nv*nv,infpar%master,infpar%lg_comm)
    call mem_alloc(w1,int(no*no*nv*nv,kind=8))
    call calculate_E2_and_permute(ppf,qqf,w1,t2,xo,yv,Gbi,Had,no,nv,nb,omega2,s)

    call mem_dealloc(ppf)
    call mem_dealloc(qqf)
    call mem_dealloc(xo)
    call mem_dealloc(yv)
    call mem_dealloc(Gbi)
    call mem_dealloc(Had)
    call mem_dealloc(w1)
end subroutine calculate_E2_and_permute_slave
#endif
