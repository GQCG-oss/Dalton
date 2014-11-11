!> @file
!> Subroutines related with construction of CC doubles residual
!> \author Patrick Ettenhuber, Marcin Ziolkowski, and Kasper Kristensen
module ccsd_module
  
  use,intrinsic :: iso_c_binding,only:c_ptr,c_f_pointer,c_associated,c_null_ptr

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
  use II_XC_interfaceModule
#ifdef VAR_ICHOR
   use IchorErimoduleHost
#endif
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif

  use lsparameters!, only: AORdefault
  use tensor_interface_module

    ! DEC DEPENDENCIES (within deccc directory)   
    ! *****************************************
  use cc_tools_module
  use dec_workarounds_module
#ifdef VAR_MPI
  use decmpi_module!, only: mpi_communicate_ccsd_calcdata,distribute_mpi_jobs
#endif

    use dec_fragment_utils
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
!         & tensor_change_atype_to_d,print_norm
    use ccintegrals!, only: get_gmo_simple,getL,dec_fock_transformation
    use pno_ccsd_module


    public :: getDoublesResidualMP2_simple,&
         & ccsd_residual_wrapper, get_ccsd_residual_integral_driven, &
         & getFockCorrection, getInactiveFockFromRI,getInactiveFock_simple, &
         & precondition_singles, precondition_doubles,get_aot1fock, get_fock_matrix_for_dec, &
         & gett1transformation, fullmolecular_get_aot1fock,calculate_E2_and_permute, &
         & get_max_batch_sizes, ccsd_energy_full_occ,print_fragment_energies_full, &
         & mo_work_dist, check_job, get_mo_ccsd_residual, &
         & wrapper_get_ccsd_batch_sizes, yet_another_ccsd_residual,&
         & RN_RESIDUAL_INT_DRIVEN, RN_YET_ANOTHER_RES
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
    
    INTEGER, PARAMETER :: RN_RESIDUAL_INT_DRIVEN = 1
    INTEGER, PARAMETER :: RN_YET_ANOTHER_RES     = 2

contains



function precondition_doubles_newarr(omega2,ppfock,qqfock,loc) result(prec)

   implicit none
   type(tensor), intent(in) :: omega2
   type(tensor), intent(inout) :: ppfock, qqfock
   logical, intent(in) :: loc
   type(tensor) :: prec
   integer, dimension(4) :: dims
   integer :: a,i,b,j
   real(realk) :: time_prec

   call time_start_phase(PHASE_WORK, twall = time_prec)

   if(omega2%mode/=4.or.ppfock%mode/=2.or.qqfock%mode/=2)then
      call lsquit("ERROR(precondition_doubles_newarr):wrong number of modes&
         & for this operation",DECinfo%output)
   endif
   dims = omega2%dims

   !make sure all data is local
   if(loc)then
      call tensor_init(prec,dims,4)

      !$OMP PARALLEL DEFAULT(NONE) SHARED(prec,dims,omega2,ppfock,qqfock) &
      !$OMP PRIVATE(i,j,a,b)

      !$OMP DO COLLAPSE(4)
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
      !$OMP END DO
      !$OMP END PARALLEL

   else
      !make sure all data is in the correct for this routine, that is omega2 is
      !TT_TILED_DIST and ppfock%addr_p_arr and qqfock%addr_p_arr are associated

      call tensor_minit(prec,dims,4,atype="TDAR",tdims=omega2%tdim,fo=omega2%offset)
      call tensor_change_atype_to_rep(ppfock,loc)
      call tensor_change_atype_to_rep(qqfock,loc)
      call precondition_doubles_parallel(omega2,ppfock,qqfock,prec)
      call tensor_change_atype_to_d(ppfock)
      call tensor_change_atype_to_d(qqfock)

   endif

   if(DECinfo%PL>2)then
      call time_start_phase(PHASE_WORK, ttot = time_prec, &
         &labelttot = 'PREC: PRECONDITION DOUBLES:', output = DECinfo%output)
   else
      call time_start_phase(PHASE_WORK)
   endif
end function precondition_doubles_newarr


!> \brief Preconditioning of doubles residual interface
!> \author Kasper Kristensen
!> \date October 2010
function precondition_doubles_oldarr(omega2,ppfock,qqfock) result(prec)

   implicit none
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

   implicit none
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

    implicit none
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





  !subroutine get_ccsd_residual_integral_driven_oldtensor_wrapper(deltafock,omega2,t2,fock,govov,nocc,nvirt,&
  !     & ppfock,qqfock,pqfock,qpfock,xocc,xvirt,yocc,yvirt,nbas,MyLsItem, omega1,iter)
  !     implicit none
  !     type(array2),intent(in) :: deltafock
  !     type(array4), intent(inout) :: omega2,t2
  !     type(array2), intent(inout) :: ppfock, qqfock,pqfock,qpfock,omega1,fock
  !     type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt
  !     type(lsitem), intent(inout) :: MyLsItem
  !     integer,intent(in) :: nbas
  !     integer,intent(in) :: nocc
  !     integer,intent(in) :: nvirt
  !     integer,intent(in) :: iter
  !     real(realk),target,intent(inout) :: govov(nvirt*nocc*nvirt*nocc)
  !     real(realk),pointer :: t2_p(:),xo_p(:),xv_p(:),yo_p(:),yv_p(:)
  !     type(tensor) :: t2a,ga,o2

  !    ! get t2 and omega2 into the correct order
  !    call array4_reorder(t2,[1,3,2,4])
  !    call ass_D4to1(t2%val,t2_p,t2%dims)
  !    t2a=tensor_init(t2%dims,4)
  !    call memory_deallocate_tensor_dense(t2a)
  !    call ass_D4to1(t2%val,t2a%elm1,t2%dims)
  !    call assoc_ptr_arr(t2a)

  !    call array4_reorder(omega2,[1,3,2,4])
  !    o2=tensor_init(omega2%dims,4)
  !    call memory_deallocate_tensor_dense(o2)
  !    call ass_D4to1(omega2%val,o2%elm1,omega2%dims)
  !    call assoc_ptr_arr(o2)
  !
  !    ga=tensor_init([nocc,nvirt,nocc,nvirt],4)
  !    call memory_deallocate_tensor_dense(ga)
  !    ga%elm1 => govov
  !    call assoc_ptr_arr(ga)

  !    call ass_D2to1(xocc%val,xo_p,xocc%dims)
  !    call ass_D2to1(xvirt%val,xv_p,xvirt%dims)
  !    call ass_D2to1(yocc%val,yo_p,yocc%dims)
  !    call ass_D2to1(yvirt%val,yv_p,yvirt%dims)
  !    !call get_ccsd_residual_integral_driven(deltafock%val,omega2%val,t2_p,fock%val,govov,nocc,nvirt,&
  !    ! & ppfock%val,qqfock%val,pqfock%val,qpfock%val,xo_p,xv_p,yo_p,yv_p,nbas,MyLsItem,&
  !    ! & omega1%val,iter)
  !    call get_ccsd_residual_integral_driven(deltafock%val,o2,t2a,fock%val,ga,nocc,nvirt,&
  !     & ppfock%val,qqfock%val,pqfock%val,qpfock%val,xo_p,xv_p,yo_p,yv_p,nbas,MyLsItem,&
  !     & omega1%val,iter,.true.)

  !    nullify(t2a%elm1)
  !    nullify(t2a%elm4)
  !    call tensor_free(t2a)
  !    nullify(o2%elm1)
  !    nullify(o2%elm4)
  !    call tensor_free(o2)
  !    nullify(ga%elm1)
  !    nullify(ga%elm4)
  !    call tensor_free(ga)

  !    !reorder ampitudes and residuals for use within the solver
  !    call array4_reorder(t2,[1,3,2,4])
  !    call array4_reorder(omega2,[1,3,2,4])
  !    nullify(xo_p)
  !    nullify(yo_p)
  !    nullify(xv_p)
  !    nullify(yv_p)
  !end subroutine get_ccsd_residual_integral_driven_oldtensor_wrapper
  subroutine ccsd_residual_wrapper(ccmodel,delta_fock,omega2,t2,&
             & fock,iajb,no,nv,ppfock,qqfock,pqfock,qpfock,xo,&
             & xv,yo,yv,nb,MyLsItem,omega1,t1,pgmo_diag,pgmo_up,&
             & MOinfo,mo_ccsd,pno_cv,pno_s,nspaces,iter,local,use_pnos,rest,frag)
    implicit none
    !> CC model
    integer,intent(in)    :: ccmodel
    type(tensor)              :: delta_fock
    type(tensor)              :: omega2
    type(tensor)              :: t2
    type(tensor)              :: fock
    type(tensor)              :: iajb
    integer,intent(in)       :: no,nv
    type(tensor)              :: ppfock
    type(tensor)              :: qqfock
    type(tensor)              :: pqfock
    type(tensor)              :: qpfock
    type(tensor)              :: xo
    type(tensor)              :: xv
    type(tensor)              :: yo
    type(tensor)              :: yv
    integer,intent(in)       :: nb
    integer,intent(in)       :: nspaces
    type(lsitem)             :: MyLsItem
    type(tensor)              :: omega1
    type(tensor)              :: t1
    type(tensor)              :: pgmo_diag
    type(tensor)              :: pgmo_up
    type(MObatchInfo)        :: MOinfo
    logical,intent(in)       :: mo_ccsd
    type(PNOSpaceInfo)       :: pno_cv(:), pno_S(:)
    integer,intent(in)       :: iter
    logical,intent(in)       :: local
    logical,intent(in)       :: use_pnos
    logical,intent(inout)    :: rest
    type(decfrag),intent(in),optional :: frag
    !internal variables
    logical :: parent
    integer :: lg_me,lg_nnod
    type(tensor) :: tloc,oloc,gloc
    integer   :: tdi(4), odi(4)
    integer   :: ttd(4), otd(4)
    character(4) :: ats, aos

#ifdef MOD_UNRELEASED
    if(use_pnos)then

       !TODO: remove these sortings
       call tensor_init( tloc, [nv,no,nv,no], 4 )
       call tensor_cp_data( t2,   tloc, order = [1,3,2,4] )
       tdi = t2%dims
       if(.not.local)ttd = t2%tdim
       ats = t2%atype
       call tensor_free( t2 )

       call tensor_init( oloc, [nv,no,nv,no], 4 )
       call tensor_zero( oloc )
       odi = omega2%dims
       if(.not.local)otd = omega2%tdim
       aos = omega2%atype
       call tensor_free( omega2 )


       call tensor_init(gloc,[nv,no,nv,no],4)
       call tensor_cp_data( iajb, gloc )

       call get_ccsd_residual_pno_style(t1%elm1,tloc%elm1,omega1%elm1,&
          &oloc%elm1,gloc%elm1,no,nv,nb,xo%elm1,xv%elm1,yo%elm1,yv%elm1,mylsitem,&
          &present(frag),pno_cv,pno_S,nspaces,ppfock%elm1,&
          &qqfock%elm1,delta_fock%elm1,iter,f=frag)

       !TODO: remove these sortings
       call tensor_minit( t2, tdi, 4, local = local, atype = ats, tdims=ttd )
       call tensor_cp_data( tloc, t2,     order = [1,3,2,4] )
       call tensor_free(tloc)

       call tensor_minit( omega2, odi, 4, local = local, atype = aos, tdims=otd )
       call tensor_cp_data( oloc, omega2, order = [1,3,2,4] )
       call tensor_free(oloc)

       call tensor_free(gloc)

    else

       if (mo_ccsd) then 

          call get_mo_ccsd_residual(ccmodel,pgmo_diag,pgmo_up,t1,omega1,t2,omega2,&
             & iajb,nb,no,nv,iter,MOinfo,mylsitem,xo%elm2,xv%elm2,yo%elm2,yv%elm2,&
             & delta_fock,fock,ppfock, pqfock,qpfock,qqfock,local)

       else 

          if (DECinfo%force_scheme.and. DECinfo%en_mem == 0) then

             call yet_another_ccsd_residual(ccmodel,delta_fock%elm1,omega2,t2,&
                & fock%elm1,iajb,no,nv,ppfock%elm1,qqfock%elm1,pqfock%elm1,qpfock%elm1,&
                & xo%elm1,xv%elm1,yo%elm1,yv%elm1,nb,MyLsItem,omega1%elm1,iter,local,&
                & rest=rest)

          else

             call get_ccsd_residual_integral_driven(ccmodel,delta_fock%elm1,omega2,t2,&
                & fock%elm1,iajb,no,nv,ppfock%elm1,qqfock%elm1,pqfock%elm1,qpfock%elm1,&
                & xo%elm1,xv%elm1,yo%elm1,yv%elm1,nb,MyLsItem,omega1%elm1,iter,local,&
                & rest=rest)

          endif

       end if
    endif
#else
    call get_ccsd_residual_integral_driven(ccmodel,delta_fock%elm1,omega2,t2,&
       & fock%elm1,iajb,no,nv,ppfock%elm1,qqfock%elm1,pqfock%elm1,qpfock%elm1,xo%elm1,&
       & xv%elm1,yo%elm1,yv%elm1,nb,MyLsItem,omega1%elm1,iter,local,rest=rest)
#endif

 end subroutine ccsd_residual_wrapper



  !> \brief Get CCSD residual in an integral direct fashion.
  !> \author Patrick Ettenhuber
  !> \date December 2012
  !
  !deltafock = on input deltafock is the difference fock matrix of one fragment with
  !respect to the fock matrix of the full molecule, this only has relevance for
  !DEC calculations, else it is 0
  !
  !omega2 = is the residual defined as type array. if local = .true. the array
  !is assumed to be in local memory stored in the elms1 variable, else it is
  !assumed to be in parallel distributed memory. the order is [a,b,i,j]
  !
  !t2 = are the amplitudes as array type, the "local" variable gives the assumed
  !data distribution as for omega2, the order is [a,b,i,j]
  !
  !fock = is the ao Fock matrix which is only needed in case of a CC2
  !calculation
  !
  !govov = is an mo-integral matrix in the distribution dictated by "local" as
  !t2 and omega2
  !
  !no = number of occupied orbitals in the system
  !
  !nv = number of virtual orbitals in the system
  !
  !nb = number of basis functions in the system
  !
  !ppfock = is the t1 transformed occ occ fock matrix on output
  !
  !qqfock = is the t1 transformed virt virt fock matrix on output
  !
  !pqfock = is the t1 transformed occ virt fock matrix on output
  !
  !qpfock = is the t1 transformed virt occ fock matrix on output
  !
  !xo,xv,yo,yv = are the lambda particle and hole matrices for the occupied and
  !virtual parts
  !
  !MyLSITEM = integral information
  !
  !omega1 = is the singles residual as [a,i]
  !
  !iter = the iteration number of the current cc iteration
  !
  !local = tells the routine (how) to handle the memory distriution
  !
  !rest = tells the routine wheter the calculation has been restarted from
  !amplitude files
  subroutine get_ccsd_residual_integral_driven(ccmodel,deltafock,omega2,t2,fock,govov,no,nv,&
        ppfock,qqfock,pqfock,qpfock,xo,xv,yo,yv,nb,MyLsItem, omega1,iter,local,rest)
     implicit none

     !> CC model
     integer,intent(in) :: ccmodel
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
     type(tensor),intent(inout) :: omega2
     !> the current guess amplitudes
     !real(realk),pointer,intent(in) ::t2(:)
     type(tensor),intent(inout) :: t2
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
     type(tensor),intent(inout) :: govov
     logical, intent(in) :: local
     !logical that specifies whether the amplitudes were read
     logical, optional, intent(inout) :: rest

     ! elementary types needed for the calculation
     type(mpi_realk)      :: gvvoo,gvoov,tpl,tmi,w0,w1,w2,w3,uigcj
     real(realk), pointer :: Had(:), t2_d(:,:,:,:), Gbi(:)
     type(c_ptr) :: Hadc,t2_dc, Gbic
     integer(kind=ls_mpik) :: Hadw,t2_dw,Gbiw,gvvoow,gvoovw

     integer(kind=8) :: w0size,w1size,w2size,w3size,neloc

     ! Variables for mpi
     logical :: master,lg_master,parent
     integer :: fintel,nintel,fe,ne
     integer(kind=ls_mpik) :: nnod
     real(realk) :: startt, stopp
     integer(kind=ls_mpik) :: ierr

     integer :: sio4_mode, sio4_dims(4),sio4_tdim(4) 
     type(tensor) :: u2, sio4
     type(tensor) :: gvoova,gvvooa
     !special arrays for scheme=1
     type(tensor) :: t2jabi,u2kcjb
     integer,pointer       :: tasks(:)
     type(c_ptr)           :: tasksc
     integer(kind=ls_mpik) :: tasksw,taskslw
     integer(kind=ls_mpik) :: lg_me,lg_nnod
     integer(kind=8)       :: len81,len82
     integer(kind=4)       :: len41,len42
     integer               :: lenI1,lenI2
     integer,parameter :: inflen=5
     real(realk)       :: inf(inflen)
#ifdef VAR_MPI
     ! stuff for direct communication
     integer(kind=ls_mpik) :: gvvoo_w, gvoov_w
     integer(kind=ls_mpik) :: hstatus, nctr,mode
     integer :: rcnt(infpar%lg_nodtot),dsp(infpar%lg_nodtot)
     character*(MPI_MAX_PROCESSOR_NAME) :: hname
     real(realk),pointer :: buf1(:), buf2(:), buf3(:)
     !integer(kind=ls_mpik),pointer :: tasksw(:)
     integer(kind=8) :: maxts,nbuff
#endif
     logical :: lock_outside

     type(matrix) :: Dens, iFock
     ! CHECKING and MEASURING variables
     integer(kind=long) :: maxsize64,dummy64
     integer :: myload,nelms,n4
     real(realk) :: tcpu, twall,tcpu1,twall1,tcpu2,twall2, deb1,deb2,ActuallyUsed
     real(realk) :: MemFree,MemFree2,MemFree3,MemFree4
     real(realk) :: tcpu_end,twall_end,time_a, time_c, time_d,time_singles
     real(realk) :: time_doubles,timewall_start,wait_time,max_wait_time,min_wait_time,ave_wait_time
     integer     :: scheme
     integer(kind=8) :: els2add
     logical :: memfound

     ! variables used for BATCH construction and INTEGRAL calculation
     integer :: alphaB,gammaB,dimAlpha,dimGamma
     integer :: dim1,dim2,dim3,K,MinAObatch
     integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
     integer :: iorb,nthreads,idx,residual_nr
#ifdef VAR_ICHOR
     type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
     type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
     integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
     logical :: MoTrans, NoSymmetry,SameMol
#else
     Character(80)        :: FilenameCS,FilenamePS
     Character(80),pointer:: BatchfilenamesCS(:,:)
     Character(80),pointer:: BatchfilenamesPS(:,:)
     logical :: FoundInMem,doscreen
     integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
     integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:)
     TYPE(DECscreenITEM)    :: DecScreen
#endif
     integer, pointer :: batchdimAlpha(:), batchdimGamma(:)     
     type(batchtoorb), pointer :: batch2orbAlpha(:)
     type(batchtoorb), pointer :: batch2orbGamma(:)
     Character            :: INTSPEC(5)
     logical :: fullRHS
     integer :: MaxAllowedDimAlpha,MaxActualDimAlpha,nbatchesAlpha
     integer :: MaxAllowedDimGamma,MaxActualDimGamma,nbatchesGamma

     integer :: a,b,i,j,l,m,n,c,d,fa,fg,la,lg,worksize
     integer :: nb2,nb3,nv2,no2,b2v,o2v,v2o,no3,vs,os
     integer(kind=8) :: nb4,o2v2,no4,buf_size
     integer :: tlen,tred,nor,nvr,goffs,aoffs
     integer :: prev_alphaB,mpi_buf,ccmodel_copy
     logical :: jobtodo,first_round,dynamic_load,restart,print_debug
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !TEST AND DEVELOPMENT VARIABLES!!!!!
     real(realk) :: op_start,op_stop, dt_last_phase
     real(realk) :: time_init,    time_init_work,    time_init_comm, time_cndonly
     real(realk) :: time_intloop, time_intloop_work, time_intloop_comm, time_intloop_idle, time_intloop_int
     real(realk) :: time_intloop_B1work, time_intloop_B1comm, time_intloop_stop
     real(realk) :: time_fock_mat,commtime,time_reduction2,time_Bonly
     real(realk) :: time_Bcnd,time_Bcnd_work,time_Bcnd_comm, time_Bcnd_idle
     real(realk) :: time_cnd,time_cnd_work,time_cnd_comm,time_get_ao_fock, time_get_mo_fock
     real(realk) :: time_Esing,time_Esing_work,time_Esing_comm, time_Esing_idle
     real(realk) :: unlock_time, waiting_time, flushing_time
     real(realk) :: phase_counters_int_dir(nphases)
     integer :: testmode(4)
     integer(kind=long) :: xyz,zyx1,zyx2
     logical :: debug
     character(tensor_MSG_LEN) :: msg
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef VAR_OMP
     integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
     character(4) :: def_atype
     integer, parameter :: bs = 1

     !init timing variables
     call time_start_phase(PHASE_WORK, twall = twall)

     time_init_work       = 0.0E0_realk
     time_init_comm       = 0.0E0_realk
     time_intloop_work    = 0.0E0_realk
     time_intloop_comm    = 0.0E0_realk
     time_intloop_int     = 0.0E0_realk
     time_intloop_B1work  = 0.0E0_realk
     time_intloop_B1comm  = 0.0E0_realk
     time_intloop_idle    = 0.0E0_realk
     time_cnd_work        = 0.0E0_realk
     time_cnd_comm        = 0.0E0_realk
     time_Bcnd_work       = 0.0E0_realk
     time_Bcnd_comm       = 0.0E0_realk
     time_Esing_work      = 0.0E0_realk
     time_Esing_comm      = 0.0E0_realk
     commtime             = 0.0E0_realk
#ifdef VAR_MPI
     unlock_time          = time_lsmpi_win_unlock 
     waiting_time         = time_lsmpi_wait
     flushing_time        = time_lsmpi_win_flush
#else
     unlock_time          = 0.0E0_realk 
     waiting_time         = 0.0E0_realk
     flushing_time        = 0.0E0_realk
#endif


     ! Set default values for the path throug the routine
     ! **************************************************
     restart                  = .false.
     if(present(rest))restart = rest
     scheme                   = 0
     dynamic_load             = DECinfo%dyn_load
     startt                   = 0.0E0_realk
     stopp                    = 0.0E0_realk
     print_debug              = (DECinfo%PL>3)
     debug                    = .false.

     ! Set some shorthand notations
     ! ****************************
     nb2                      = nb*nb
     nb3                      = nb2*nb
     nb4                      = int((i8*nb3)*nb,kind=long)
     nv2                      = nv*nv
     no2                      = no*no
     no3                      = no2*no
     no4                      = int((i8*no3)*no,kind=long)
     b2v                      = nb2*nv
     o2v                      = no2*nv
     v2o                      = nv2*no
     o2v2                     = int((i8*nv2)*no2,kind=long)
     nor                      = no*(no+1)/2
     nvr                      = nv*(nv+1)/2
     vs                       = t2%tdim(1)
     os                       = t2%tdim(3)
     residual_nr              = RN_RESIDUAL_INT_DRIVEN

     
     ! Memory info
     ! ***********
     call get_currently_available_memory(MemFree)

     ! Set integral info
     ! *****************
     INTSPEC(1)               = 'R' !R = Regular Basis set on the 1th center 
     INTSPEC(2)               = 'R' !R = Regular Basis set on the 2th center 
     INTSPEC(3)               = 'R' !R = Regular Basis set on the 3th center 
     INTSPEC(4)               = 'R' !R = Regular Basis set on the 4th center 
     INTSPEC(5)               = 'C' !C = Coulomb operator
#ifdef VAR_ICHOR
     iprint = 0           !print level for Ichor Integral code
     MoTrans = .FALSE.    !Do not transform to MO basis! 
     NoSymmetry = .FALSE. !Use Permutational Symmetry! 
     SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
     !Determine the full number of AO batches - not to be confused with the batches of AOs
     !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
     iAO = 1
     call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
#else
     doscreen = MyLsItem%setting%scheme%cs_screen.OR.MyLsItem%setting%scheme%ps_screen
#endif

     ! Set MPI related info
     ! ********************
     master                   = .true.
     lg_master                = .true.
     parent                   = .true.
     lg_me                    = int(0,kind=ls_mpik)
     lg_nnod                  = 1
#ifdef VAR_MPI
     lg_me                    = infpar%lg_mynum
     lg_nnod                  = infpar%lg_nodtot
     lock_outside             = DECinfo%CCSD_NO_DEBUG_COMM
     mode                     = MPI_MODE_NOCHECK

     parent                   = (infpar%parent_comm==MPI_COMM_NULL)

     lg_master                = (lg_me == 0)
     master                   = lg_master.and.parent

     call get_int_dist_info(o2v2,fintel,nintel)

     StartUpSlaves: if(master .and. lg_nnod>1) then
        call time_start_phase(PHASE_COMM)
        call ls_mpibcast(CCSDDATA,infpar%master,infpar%lg_comm)
        ccmodel_copy = ccmodel
        call mpi_communicate_ccsd_calcdata(ccmodel_copy,omega2,t2,govov,xo,xv,yo,yv,MyLsItem,nb,nv,no,iter,local,residual_nr)
        call time_start_phase(PHASE_WORK, at = time_init_comm)
     endif StartUpSlaves

     if(.not.local)then
        t2%access_type     = AT_ALL_ACCESS
        call tensor_lock_wins( t2, 's', mode , all_nodes = alloc_in_dummy )
        call memory_allocate_tensor_dense( t2 )
        buf_size = min(int((MemFree*0.8*1024.0_realk**3)/(8.0*t2%tsize)),5)*t2%tsize
        call mem_alloc(buf1,buf_size)
        call tensor_gather(1.0E0_realk,t2,0.0E0_realk,t2%elm1,o2v2,wrk=buf1,iwrk=buf_size)
        call mem_dealloc(buf1)
        call tensor_unlock_wins( t2, all_nodes = alloc_in_dummy, check =.not.alloc_in_dummy )
     endif
     call lsmpi_reduce_realk_min(MemFree,infpar%master,infpar%lg_comm)
#endif


     !print*,"HACK:random t"
     !if(master)call random_number(t2%elm1)


     if(master.and.print_debug)then
        call print_norm(xo,int(nb*no,kind=8)," NORM(xo)    :")
        call print_norm(xv,int(nb*nv,kind=8)," NORM(xv)    :")
        call print_norm(yo,int(nb*no,kind=8)," NORM(yo)    :")
        call print_norm(yv,int(nb*nv,kind=8)," NORM(yv)    :")
     endif

     ! Initialize stuff
#ifdef VAR_ICHOR
     nullify(AOGammabatchinfo)
     nullify(AOalphabatchinfo)    
#else
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
#endif
     nullify(Had)
     nullify(Gbi)
#ifdef VAR_MPI
     nullify(tasks)
#endif

     if(master) then
        !==================================================
        !                  Batch construction             !
        !==================================================

        ! Get free memory and determine maximum batch sizes
        ! -------------------------------------------------
#ifdef VAR_ICHOR
        !Determine the minimum allowed AObatch size MinAObatch
        !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
        !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch = 6 (the 2*(Px,Py,Pz))
        !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
        !'R'  !Specifies that it is the Regular AO basis that should be batched
        iAO = 4 !the center that the batching should occur on (they are all the same in this case)  
        call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
#else
        call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
#endif
        call get_max_batch_sizes(scheme,nb,bs,nv,vs,no,os,MaxAllowedDimAlpha,MaxAllowedDimGamma,&
        &MinAObatch,DECinfo%manual_batchsizes,iter,MemFree,.true.,els2add,local,.false.,mylsitem%setting,intspec)

        !SOME WORDS ABOUT THE CHOSEN SCHEME:
        ! Depending on the availability of memory on the nodes a certain scheme
        ! for the calculations is chosen, hereby the schemes 4-1 are the MPI
        ! schemes with decreasing memory reqirements. Hereby scheme 4 and 0 can be
        ! used without MPI. Scheme 0 should never be chosen with MPI and
        ! schemes 2,3 should never be chosen without MPI

        ! scheme 4: everything is treated locally, only the main integral driven
        !           loop is MPI-parallel
        ! scheme 3: treat govov, gvoov and gvvoo in the main part in PDM
        ! scheme 2: additionally to 3 also the amplitudes, u, the residual are
        !           treated in PDM, the strategy is to only use one V^2O^2 in 
        !           local mem
        ! scheme 1: All 4 dimensional quantities are stored in PDM

#ifndef VAR_MPI
        if(scheme==3.or.scheme==2.or.scheme==1) call lsquit("ERROR(ccsd_residual_integral_driven):wrong choice of scheme",-1)
#endif
     endif


     !all communication for MPI prior to the loop
#ifdef VAR_MPI
     call time_start_phase(PHASE_COMM, at = time_init_work )
        

     call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
     call ls_mpi_buffer(scheme,infpar%master)
     call ls_mpi_buffer(print_debug,infpar%master)
     call ls_mpi_buffer(dynamic_load,infpar%master)
     call ls_mpi_buffer(restart,infpar%master)
     call ls_mpi_buffer(MaxAllowedDimAlpha,infpar%master)
     call ls_mpi_buffer(MaxAllowedDimGamma,infpar%master)
     call ls_mpi_buffer(els2add,infpar%master)
     call ls_mpi_buffer(lock_outside,infpar%master)
     call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

     call time_start_phase(PHASE_WORK, at = time_init_comm)

     hstatus = 80
     CALL MPI_GET_PROCESSOR_NAME(hname,hstatus,ierr)


     govov%access_type  = AT_ALL_ACCESS
     omega2%access_type = AT_ALL_ACCESS

     if( alloc_in_dummy .and. scheme == 2)then
        call tensor_lock_wins(omega2,'s',all_nodes = .true.)
     endif
#endif

     !if the residual is handeled as dense, allocate and zero it, adjust the
     !access parameters to the data
     if(omega2%itype/=TT_DENSE.and.(scheme==3.or.scheme==4))then
        call memory_allocate_tensor_dense_pc(omega2)
        omega2%itype=TT_DENSE
     endif
     call tensor_zero(omega2)


     ! ************************************************
     ! * Determine batch information for Gamma batch  *
     ! ************************************************
#ifdef VAR_ICHOR
    iAO = 4 !Gamma is the 4. Center of the 4 center two electron coulomb integral
    !Determine how many batches of AOS based on the MaxAllowedDimGamma, the requested
    !size of the AO batches. iAO is the center that the batching should occur on. 
    !'R'  !Specifies that it is the Regular AO basis that should be batched 
    call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
         & nbatchesGamma,DECinfo%output)
    call mem_alloc(AOGammabatchinfo,nbatchesGamma)
    !Construct the batches of AOS based on the MaxAllowedDimGamma, the requested
    !size of the AO batches - MaxAllowedDimGamma must be unchanged since the call 
    !to determine_Ichor_nbatchesofAOS
    !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
    !So MaxActualDimGamma must be less og equal to MaxAllowedDimGamma
    call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
         & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)
#else
     ! Orbital to batch information
     ! ----------------------------
     call mem_alloc(orb2batchGamma,nb)
     call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
        & nb,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,&
        &nbatchesGamma,orb2BatchGamma,'R')
#endif
     if(master.and.DECinfo%PL>1)write(DECinfo%output,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma,&
        & 'with maximum size',MaxActualDimGamma

#ifndef VAR_ICHOR
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
#endif

     ! ************************************************
     ! * Determine batch information for Alpha batch  *
     ! ************************************************

#ifdef VAR_ICHOR
    iAO = 3 !Alpha is the 3. Center of the 4 center two electron coulomb integral
    !Determine how many batches of AOS based on the MaxAllowedDimAlpha, the requested
    !size of the AO batches. iAO is the center that the batching should occur on. 
    !'R'  !Specifies that it is the Regular AO basis that should be batched 
    call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
         & nbatchesAlpha,DECinfo%output)
    call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
    !Construct the batches of AOS based on the MaxAllowedDimAlpha, the requested
    !size of the AO batches - MaxAllowedDimAlpha must be unchanged since the call 
    !to determine_Ichor_nbatchesofAOS
    !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
    !So MaxActualDimAlpha must be less og equal to MaxAllowedDimAlpha
    call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
         & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
#else
     ! Orbital to batch information
     ! ----------------------------
     call mem_alloc(orb2batchAlpha,nb)
     call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
        & nb,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,'R')
#endif

     if(master.and.DECinfo%PL>1)write(DECinfo%output,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha&
        &, 'with maximum size',MaxActualDimAlpha

#ifndef VAR_ICHOR
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
#endif

     ! ************************************************
     ! *  Allocate matrices used in the batched loop  *
     ! ************************************************



     ! PRINT some information about the calculation
     ! --------------------------------------------
     if(master.and.DECinfo%PL>1) then
        if(scheme==4) write(DECinfo%output,'("Using memory intensive scheme (NON-PDM)")')
        if(scheme==3) write(DECinfo%output,'("Using memory intensive scheme with direct updates")')
        if(scheme==2) write(DECinfo%output,'("Using memory intensive scheme only 1x V^2O^2")')
        if(scheme==1) write(DECinfo%output,'("Using Dmitry s scheme")')
        ActuallyUsed=get_min_mem_req(no,os,nv,vs,nb,bs,MaxActualDimAlpha,MaxActualDimGamma,iter,3,scheme,.false.,&
           &mylsitem%setting,intspec)
        write(DECinfo%output,'("Using",1f8.4,"% of available Memory in part B on master")')ActuallyUsed/MemFree*100
        ActuallyUsed=get_min_mem_req(no,os,nv,vs,nb,bs,MaxActualDimAlpha,MaxActualDimGamma,iter,2,scheme,.false.,&
           &mylsitem%setting,intspec)
        write(DECinfo%output,'("Using",1f8.4,"% of available Memory in part C on master")')ActuallyUsed/MemFree*100
        ActuallyUsed=get_min_mem_req(no,os,nv,vs,nb,bs,MaxActualDimAlpha,MaxActualDimGamma,iter,4,scheme,.true.,&
           &mylsitem%setting,intspec)
     endif

     ! Use the dense amplitudes
     ! ------------------------

     !get the t+ and t- for the Kobayshi-like B2 term
     call mem_alloc( tpl, int(i8*nor*nvr,kind=long), simple = .true. )
     call mem_alloc( tmi, int(i8*nor*nvr,kind=long), simple = .true. )

!#ifdef VAR_MPI
!     call tensor_unlock_wins(t2,.true.)
!#endif
     if(scheme == 1) then
        print *, "Dmitry, I stop the program here, if your scheme is called.&
        & Here begins your work with distributing the tpl and tmi, also u2, &
        & needs to be constructed correctly"
        stop 0
     endif

     !if I am the working process, then
     call get_tpl_and_tmi_fort(t2%elm1,nv,no,tpl%d,tmi%d)

     if(master.and.print_debug)then
        call print_norm(tpl%d,int(nor*nvr,kind=8)," NORM(tpl)   :")
        call print_norm(tmi%d,int(nor*nvr,kind=8)," NORM(tmi)    :")
     endif


     !get u2 in pdm or local
     if(scheme==2)then
#ifdef VAR_MPI
        call memory_deallocate_tensor_dense(t2)

        call time_start_phase(PHASE_COMM, at = time_init_work )

        call tensor_ainit( u2, [nv,nv,no,no], 4, local=local, atype='TDAR', tdims=[vs,vs,os,os] )
        call tensor_add( u2,  2.0E0_realk, t2, a = 0.0E0_realk, order=[2,1,3,4] )
        call tensor_add( u2, -1.0E0_realk, t2, order=[2,1,4,3] )
        
        if( alloc_in_dummy )call tensor_lock_wins(u2,'s',all_nodes = .true.)

        call time_start_phase(PHASE_WORK, at = time_init_comm )
#endif

     else
        call tensor_ainit(u2, [nv,nv,no,no], 4, local=local, atype='LDAR' )
        !calculate u matrix: t[c d i j] -> t[d c i j], 2t[d c i j] - t[d c j i] = u [d c i j]
        call array_reorder_4d(  2.0E0_realk, t2%elm1,nv,nv,no,no,[2,1,3,4],0.0E0_realk,u2%elm1)
        call array_reorder_4d( -1.0E0_realk, t2%elm1,nv,nv,no,no,[2,1,4,3],1.0E0_realk,u2%elm1)
     endif

     if(print_debug) call print_norm(u2," NORM(u2)    :",print_on_rank=0)

     call mem_alloc(Had,nv*nb)
     call mem_alloc(Gbi,nb*no)


     if( CCmodel > MODEL_CC2 )then

        if(scheme==4)then
           write(def_atype,'(A4)')'LDAR'
        else if(scheme==2.or.scheme==3)then
           write(def_atype,'(A4)')'TDAR'
        endif
        call tensor_ainit(gvvooa, [nv,no,no,nv],4, local=local, atype=def_atype, tdims=[vs,os,os,vs])
        call tensor_ainit(gvoova, [nv,no,nv,no],4, local=local, atype=def_atype, tdims=[vs,os,vs,os])
        call tensor_zero(gvvooa)
        call tensor_zero(gvoova)
!        call mem_alloc(sio4,int(i8*nor*no2,kind=long))
!#ifdef VAR_MPI
!        call lsmpi_win_create(sio4%d,sio4w,int(i8*nor*no2,kind=long),infpar%lg_comm)
!#endif
        if(scheme == 4 .or. scheme == 3)then

           sio4_mode = 3
           sio4_dims(1:sio4_mode) = [no,no,nor]
           sio4_tdim(1:sio4_mode) = [os,os,nor]
           write(def_atype,'(A4)')'LDAR'

        else if(scheme == 2)then

           sio4_mode = 4
           sio4_dims(1:sio4_mode) = [no,no,no,no]
           sio4_tdim(1:sio4_mode) = [os,os,os,os]
           write(def_atype,'(A4)')'TDAR'
        endif

        call tensor_ainit(sio4,sio4_dims(1:sio4_mode),sio4_mode,local=local,atype=def_atype,tdims = sio4_tdim(1:sio4_mode))
        call tensor_zero(sio4)


#ifdef VAR_MPI
        if( alloc_in_dummy )then
           if(scheme==2.or.scheme==3)then
              call tensor_lock_wins(gvvooa,'s',all_nodes = .true.)
              call tensor_lock_wins(gvoova,'s',all_nodes = .true.)
           endif
           if(scheme==2)then
              call tensor_lock_wins(sio4,  's',all_nodes = .true.)
           endif
        endif
#endif

     endif


     !zero the matrix
     !$OMP WORKSHARE
     Had = 0.0E0_realk
     Gbi = 0.0E0_realk
     !$OMP END WORKSHARE

     !call get_currently_available_memory(MemFree2)
     !call get_available_memory(DECinfo%output,MemFree4,memfound,suppress_print=.true.)

     ! allocate working arrays depending on the batch sizes
     w0size = get_wsize_for_ccsd_int_direct(0,no,os,nv,vs,nb,0,&
        &MaxActualDimAlpha,MaxActualDimGamma,scheme,mylsitem%setting,intspec)
     call mem_alloc( w0, w0size , simple = .false. )

     w1size = get_wsize_for_ccsd_int_direct(1,no,os,nv,vs,nb,0,&
        &MaxActualDimAlpha,MaxActualDimGamma,scheme,mylsitem%setting,intspec)
     call mem_alloc( w1, w1size , simple = .false.)

     w2size = get_wsize_for_ccsd_int_direct(2,no,os,nv,vs,nb,0,&
        &MaxActualDimAlpha,MaxActualDimGamma,scheme,mylsitem%setting,intspec)
     call mem_alloc( w2, w2size , simple = .false. )

     w3size = get_wsize_for_ccsd_int_direct(3,no,os,nv,vs,nb,0,&
        &MaxActualDimAlpha,MaxActualDimGamma,scheme,mylsitem%setting,intspec)
     call mem_alloc( w3, w3size , simple = .false. )

     !call get_currently_available_memory(MemFree3)

#ifdef VAR_MPI
     !print *,infpar%lg_mynum,"have",MemFree,MemFree2,MemFree4,MemFree3
     !call lsmpi_barrier(infpar%lg_comm)

     !print *,infpar%lg_mynum,"first touching w0",w0%n,(w0%n*8.0E0_realk)/1024.0**3
     !w0%d=0.0E0_realk
     !call lsmpi_barrier(infpar%lg_comm)
     !print *,infpar%lg_mynum,"first touching w1",w1%n,(w1%n*8.0E0_realk)/1024.0**3
     !w1%d=0.0E0_realk
     !call lsmpi_barrier(infpar%lg_comm)
     !print *,infpar%lg_mynum,"first touching w2",w2%n,(w2%n*8.0E0_realk)/1024.0**3
     !w2%d=0.0E0_realk
     !call lsmpi_barrier(infpar%lg_comm)
     !print *,infpar%lg_mynum,"first touching w3",w3%n,(w3%n*8.0E0_realk)/1024.0**3
     !w3%d=0.0E0_realk
     !call lsmpi_barrier(infpar%lg_comm)
#endif

     !allocate semi-permanent storage arrays for loop
     !print *,"allocing help things:",o2v*MaxActualDimGamma*2,&
     !      &(8.0E0_realk*o2v*MaxActualDimGamma*2)/(1024.0E0_realk*1024.0E0_realk*1024.0E0_realk)
     call mem_alloc( uigcj, int((i8*o2v)*MaxActualDimGamma,kind=8))


#ifdef VAR_ICHOR
     !Calculate Screening integrals 
     SameMOL = .TRUE. !Specifies same molecule on all centers 
     call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
#else
     ! This subroutine builds the full screening matrix.
     call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
          & nbatchesAlpha,nbatchesGamma,INTSPEC)
     IF(mylsitem%setting%scheme%cs_screen .OR. mylsitem%setting%scheme%ps_screen)THEN
        call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
             & nb,nbatchesAlpha,nbatchesGamma,&
             & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
             & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
        call II_getBatchOrbitalScreenK(DecScreen,mylsitem%setting,&
             & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
             & batchindexAlpha,batchindexGamma,&
             & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
     ENDIF
     !setup LHS screening - the full AO basis is used so we can use the
     !                      full matrices:        FilenameCS and FilenamePS
     !Note that it is faster to calculate the integrals in the form
     !(dimAlpha,dimGamma,nbasis,nbasis) so the full AO basis is used on the RHS
     !but the integrals is stored and returned in (nbasis,nbasis,dimAlpha,dimGamma)
#endif

#ifdef VAR_OMP
     nthreads=OMP_GET_MAX_THREADS()
     if(master.and.DECinfo%PL>2)write(DECinfo%output,*)&
        & 'Starting CCSD residuals - OMP. Number of threads: ', OMP_GET_MAX_THREADS()
#else
     nthreads=1
     if(master.and.DECinfo%PL>2)write(DECinfo%output,*) &
        &'Starting CCSD integral/amplitudes - NO OMP!'
#endif


#ifdef VAR_MPI
     if(.not.dynamic_load)then

        ! Calculate the batches for a good load balance
        lenI2 = nbatchesAlpha*nbatchesGamma
        !call mem_alloc( tasks, tasksc, lenI2 )
        call mem_alloc( tasks, lenI2 )

        myload = 0
        tasks  = 0

#ifdef VAR_ICHOR
        call mem_alloc(batchdimAlpha,nbatchesAlpha)
        do idx=1,nbatchesAlpha
           batchdimAlpha(idx) = AOAlphabatchinfo(idx)%dim 
        enddo
        call mem_alloc(batch2orbAlpha,nbatchesAlpha)
        do idx=1,nbatchesAlpha
           call mem_alloc(batch2orbAlpha(idx)%orbindex,1)
           batch2orbAlpha(idx)%orbindex(1) = AOAlphabatchinfo(idx)%orbstart
           batch2orbAlpha(idx)%norbindex = 1
        end do
        call mem_alloc(batchdimGamma,nbatchesGamma)
        do idx=1,nbatchesGamma
           batchdimGamma(idx) = AOGammabatchinfo(idx)%dim 
        enddo
        call mem_alloc(batch2orbGamma,nbatchesGamma)
        do idx=1,nbatchesGamma
           call mem_alloc(batch2orbGamma(idx)%orbindex,1)
           batch2orbGamma(idx)%orbindex(1) = AOGammabatchinfo(idx)%orbstart
           batch2orbGamma(idx)%norbindex = 1
        end do
        call distribute_mpi_jobs(tasks,nbatchesAlpha,nbatchesGamma,&
             & batchdimAlpha,batchdimGamma,myload,lg_nnod,lg_me,scheme,&
             & no,nv,nb,batch2orbAlpha,batch2orbGamma)
        call mem_dealloc(batchdimAlpha)
        call mem_dealloc(batchdimGamma)
        do idx=1,nbatchesAlpha
           call mem_dealloc(batch2orbAlpha(idx)%orbindex)
        end do
        call mem_dealloc(batch2orbAlpha)
        do idx=1,nbatchesGamma
           call mem_dealloc(batch2orbGamma(idx)%orbindex)
        end do
        call mem_dealloc(batch2orbGamma)
#else
        call distribute_mpi_jobs(tasks,nbatchesAlpha,nbatchesGamma,batchdimAlpha,&
             &batchdimGamma,myload,lg_nnod,lg_me,scheme,no,nv,nb,batch2orbAlpha,&
             &batch2orbGamma)
#endif

     else

        lenI2 = nbatchesGamma

        call mem_alloc( tasks, tasksc, lenI2 ) 

        tasks = 0
        if(lg_me == 0) tasks(1) = lg_nnod

        call lsmpi_win_create(tasks,tasksw,nbatchesGamma,infpar%lg_comm)
#ifdef VAR_HAVE_MPI3
        call lsmpi_win_lock_all(tasksw,ass=mode)
#endif

     endif

     
     if(master.and.DECinfo%PL>2)then
        write(*,'("CCSD time in lsmpi_win_unlock phase A",g10.3)') time_lsmpi_win_unlock - unlock_time
        write(*,'("CCSD time in lsmpi_wait       phase A",g10.3)') time_lsmpi_wait       - waiting_time
        write(*,'("CCSD time in lsmpi_win_flush  phase A",g10.3)') time_lsmpi_win_flush  - flushing_time
     endif
     unlock_time   = time_lsmpi_win_unlock 
     waiting_time  = time_lsmpi_wait
     flushing_time = time_lsmpi_win_flush
#endif

     myload = 0

     !TIMING
     call time_start_phase(PHASE_WORK, at = time_init_work, twall = time_intloop )
     time_init = time_init_work + time_init_comm
     if(master.and.DECinfo%PL>2)then
        write( *,'("-------------------------------------------------------------------")')
        write( *,'("CCSD residual init work:",g10.3,"s, comm:",g10.3,"s")')time_init_work,time_init_comm
     endif
     call time_phases_get_current(current_wt=phase_counters_int_dir)

     fullRHS=(nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)

     !**********************************
     ! Begin the loop over gamma batches
     !**********************************

     first_round=.false.
     if(dynamic_load)first_round=.true.

     BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
#ifdef VAR_ICHOR
        dimGamma     = AOGammabatchinfo(gammaB)%dim      ! Dimension of gamma batch
        GammaStart   = AOGammabatchinfo(gammaB)%orbstart ! First orbital index in gamma batch
        GammaEnd     = AOGammabatchinfo(gammaB)%orbEnd   ! Last orbital index in gamma batch
        AOGammaStart = AOGammabatchinfo(gammaB)%AOstart  ! First AO batch index in gamma batch
        AOGammaEnd   = AOGammabatchinfo(gammaB)%AOEnd    ! Last AO batch index in gamma batch
#else
        dimGamma     = batchdimGamma(gammaB)                         ! Dimension of gamma batch
        GammaStart   = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
        GammaEnd     = batch2orbGamma(gammaB)%orbindex(dimGamma)     ! Last index in gamma batch
#endif
        !short hand notation
        fg         = GammaStart
        lg         = dimGamma

        !Lambda^h [gamma d] u[d c i j] = u [gamma c i j]
        if(scheme==2)then
#ifdef VAR_MPI
           call time_start_phase(PHASE_COMM, at = time_intloop_work )
           !AS LONG AS THE INTEGRALS ARE WRITTEN IN W1 we might unlock here
           if(gammaB /= 1)then
              if( alloc_in_dummy )then
                 call lsmpi_win_flush(omega2%wi(1),local=.false.)
              else
                 if( lock_outside )call tensor_unlock_wins(omega2,.true.)
              endif
           endif

           call tensor_convert(u2,w2%d,wrk=w1%d,iwrk=w1%n)

           if( alloc_in_dummy )call lsmpi_win_flush(u2%wi(1),local=.true.)

           call time_start_phase(PHASE_WORK, at = time_intloop_comm)

           call dgemm('n','n',lg,o2v,nv,1.0E0_realk,yv(fg),nb,w2%d,nv,0.0E0_realk,w1%d,lg)
#endif 
        else
           call dgemm('n','n',lg,o2v,nv,1.0E0_realk,yv(fg),nb,u2%elm1,nv,0.0E0_realk,w1%d,lg)
        endif
        !u [gamma c i j ] -> u [i gamma c j]
        call array_reorder_4d(1.0E0_realk,w1%d,lg,nv,no,no,[3,1,2,4],0.0E0_realk,uigcj%d)

        alphaB=0
           

        !**********************************
        ! Begin the loop over alpha batches
        !**********************************


        BatchAlpha: do while(alphaB<=nbatchesAlpha) ! AO batches

           !check if the current job is to be done by current node
           call time_start_phase(PHASE_COMM, at = time_intloop_work)
           call check_job(scheme,first_round,dynamic_load,alphaB,gammaB,nbatchesAlpha,&
              &nbatchesGamma,tasks,tasksw,print_debug)
           call time_start_phase(PHASE_WORK, at = time_intloop_comm)

           !break the loop if alpha become too large, necessary to account for all
           !of the mpi and non mpi schemes, this is accounted for, because static,
           !and dynamic load balancing are enabled
           if(alphaB>nbatchesAlpha) exit

#ifdef VAR_ICHOR
           dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
           AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
           AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
           AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
           AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
#else
           dimAlpha   = batchdimAlpha(alphaB)                              ! Dimension of alpha batch
           AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
           AlphaEnd   = batch2orbAlpha(alphaB)%orbindex(dimAlpha)          ! Last index in alpha batch
#endif

           !short hand notation
           fa         = AlphaStart
           la         = dimAlpha
           myload     = myload + la * lg


           !u[k gamma  c j] * Lambda^p [alpha j] ^T = u [k gamma c alpha]
           call dgemm('n','t',no*nv*lg,la,no,1.0E0_realk,uigcj%d,no*nv*lg,xo(fa),nb,0.0E0_realk,w1%d,nv*no*lg)
           call lsmpi_poke()
           !Transpose u[k gamma c alpha]^T -> u[c alpha k gamma]
           call array_reorder_4d(1.0E0_realk,w1%d,no,lg, nv,la,[3,4,1,2],0.0E0_realk,w3%d)
           call lsmpi_poke()

           !print*,"GAMMA:",fg,nbatchesGamma,"ALPHA:",fa,nbatchesAlpha
           !print*,"--------------------------------------------------"

           !setup RHS screening - here we only have a set of AO basisfunctions
           !                      so we use the batchscreening matrices.
           !                      like BatchfilenamesCS(alphaB,gammaB)
           !Note that it is faster to calculate the integrals in the form
           !(dimAlpha,dimGamma,nbasis,nbasis) so the subset of the AO basis is used on the LHS
           !but the integrals is stored and returned in (nbasis,nbasis,dimAlpha,dimGamma)
           call time_start_phase(PHASE_WORK, at = time_intloop_work)
#ifdef VAR_ICHOR
           dim1 = nb*nb*dimAlpha*dimGamma   ! dimension for integral array
           call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,nb,nb,dimAlpha,dimGamma,&
                & w1%d,INTSPEC,FULLRHS,1,nAObatches,1,nAObatches,AOAlphaStart,AOAlphaEnd,&
                & AOGammaStart,AOGammaEnd,MoTrans,nb,nb,dimAlpha,dimGamma,NoSymmetry)
#else
           IF(doscreen) Mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
           IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGab(alphaB,gammaB)%p
           ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
           ! ************************************************************************************
           dim1 = nb*nb*dimAlpha*dimGamma   ! dimension for integral array
           ! Store integral in tmp1(1:dim1) array in (beta,delta,alphaB,gammaB) order
           call LSTIMER('START',tcpu1,twall1,DECinfo%output)
           !Mylsitem%setting%scheme%intprint=6
           call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, Mylsitem%setting, w1%d,batchindexAlpha(alphaB),&
              &batchindexGamma(gammaB),&
              &batchsizeAlpha(alphaB),batchsizeGamma(gammaB),nb,nb,dimAlpha,dimGamma,fullRHS,INTSPEC)
           !Mylsitem%setting%scheme%intprint=0
#endif
           call LSTIMER('START',tcpu2,twall2,DECinfo%output)

           call time_start_phase(PHASE_COMM, at = time_intloop_int)
#ifdef VAR_MPI
           !AS LONG AS THE INTEGRALS ARE WRITTEN IN W1 we might unlock here
           if( alloc_in_dummy )then
              if( scheme==2 )call lsmpi_win_flush(omega2%wi(1),local=.true.)
           else
              if( lock_outside .and. scheme==2 )call tensor_unlock_wins(omega2,.true.)
           endif
#endif
           call time_start_phase(PHASE_WORK, at = time_intloop_comm)

           call array_reorder_4d(1.0E0_realk,w1%d,nb,nb,la,lg,[4,2,3,1],0.0E0_realk,w0%d)
           call lsmpi_poke()

           ! I [gamma delta alpha beta] * Lambda^p [beta l] = I[gamma delta alpha l]
           call dgemm('n','n',lg*la*nb,no,nb,1.0E0_realk,w0%d,lg*nb*la,xo,nb,0.0E0_realk,w2%d,lg*nb*la)
           call lsmpi_poke()
           !Transpose I [gamma delta alpha l]^T -> I [alpha l gamma delta]
           call array_reorder_4d(1.0E0_realk,w2%d,lg,nb,la,no,[3,4,1,2],0.0E0_realk,w1%d)
           call lsmpi_poke()

           !u [b alpha k gamma] * I [alpha k gamma delta] =+ Had [a delta]
           call dgemm('n','n',nv,nb,lg*la*no,1.0E0_realk,w3%d,nv,w1%d,lg*la*no,1.0E0_realk,Had,nv)
           call lsmpi_poke()

           !VVOO
           if ( Ccmodel > MODEL_CC2 ) then
              !I [alpha  i gamma delta] * Lambda^h [delta j]          = I [alpha i gamma j]
              call dgemm('n','n',la*no*lg,no,nb,1.0E0_realk,w1%d,la*no*lg,yo,nb,0.0E0_realk,w3%d,la*no*lg)
              call lsmpi_poke()
              ! gvvoo = (vv|oo) constructed from w2%d                 = I [alpha i j  gamma]
              call array_reorder_4d(1.0E0_realk,w3%d,la,no,lg,no,[1,2,4,3],0.0E0_realk,w2%d)
              call lsmpi_poke()
              !I [alpha  i j gamma] * Lambda^h [gamma b]            = I [alpha i j b]
              call dgemm('n','n',la*no2,nv,lg,1.0E0_realk,w2%d,la*no2,yv(fg),nb,0.0E0_realk,w3%d,la*no2)
              call lsmpi_poke()
              !Lambda^p [alpha a]^T * I [alpha i j b]             =+ gvvoo [a i j b]
              if(scheme==4)then
                 call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w3%d,la,1.0E0_realk,gvvooa%elm1,nv)
              else if(scheme==3.or.scheme==2)then
#if VAR_MPI
                 if(.not. alloc_in_dummy .and. lock_outside) call tensor_lock_wins(gvvooa,'s',mode)
                 call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w3%d,la,0.0E0_realk,w2%d,nv)
                 call time_start_phase(PHASE_COMM, at = time_intloop_work)
                 call tensor_add(gvvooa,1.0E0_realk,w2%d,wrk=w0%d,iwrk=w0%n)
                 call time_start_phase(PHASE_WORK, at = time_intloop_comm)
#endif
              endif
              call lsmpi_poke()
           endif

           ! I [alpha l gamma delta] * Lambda^h [delta c] = I[alpha l gamma c]
           call dgemm('n','n',lg*la*no,nv,nb,1.0E0_realk,w1%d,la*no*lg,yv,nb,0.0E0_realk,w3%d,la*no*lg)
           call lsmpi_poke()
           !I [alpha l gamma c] * u [l gamma c j]  =+ Gbi [alpha j]
           call dgemm('n','n',la,no,nv*no*lg,1.0E0_realk,w3%d,la,uigcj%d,nv*no*lg,1.0E0_realk,Gbi(fa),nb)
           call lsmpi_poke()

           if ( Ccmodel > MODEL_CC2 ) then

              !Reorder I [alpha j gamma b]                      -> I [alpha j b gamma]
              call array_reorder_4d(1.0E0_realk,w3%d,la,no,lg,nv,[1,2,4,3],0.0E0_realk,w2%d)
              ! gvoov = (vo|ov) constructed from w2%d               = I [alpha j b  gamma]
              !I [alpha  j b gamma] * Lambda^h [gamma i]          = I [alpha j b i]
              call dgemm('n','n',la*no*nv,no,lg,1.0E0_realk,w2%d,la*no*nv,yo(fg),nb,0.0E0_realk,w1%d,la*no*nv)
              call lsmpi_poke()

              !Lambda^p [alpha a]^T * I [alpha j b i]             =+ gvoov [a j b i]
              if(scheme==4)then

                 call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w1%d,la,1.0E0_realk,gvoova%elm1,nv)

              else if(scheme==3.or.scheme==2)then
#ifdef VAR_MPI
                 call dgemm('t','n',nv,o2v,la,1.0E0_realk,xv(fa),nb,w1%d,la,0.0E0_realk,w2%d,nv)
                 call time_start_phase(PHASE_COMM, at = time_intloop_work)
                 if( .not. alloc_in_dummy.and.lock_outside )call tensor_lock_wins(gvoova,'s',mode)
                 call tensor_add(gvoova,1.0E0_realk,w2%d,wrk = w3%d, iwrk=w3%n)
                 call time_start_phase(PHASE_WORK, at = time_intloop_comm)
#endif
              endif
              call lsmpi_poke()
           endif

           call time_start_phase(PHASE_WORK, at = time_intloop_work)
#ifdef VAR_ICHOR
           !Build (batchA,full,batchC,full)
           call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,dimAlpha,nb,dimGamma,nb,&
                & w1%d,INTSPEC,FULLRHS,AOAlphaStart,AOAlphaEnd,1,nAObatches,AOGammaStart,AOGammaEnd,&
                & 1,nAObatches,MoTrans,dimAlpha,nb,dimGamma,nb,NoSymmetry)
#else
           IF(doscreen)Mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGabKLHS(alphaB)%p
           IF(doscreen)Mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGabKRHS(gammaB)%p

           call II_GET_DECPACKED4CENTER_K_ERI(DECinfo%output,DECinfo%output, &
              & Mylsitem%setting,w1%d,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
              & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),dimAlpha,nb,dimGamma,nb,INTSPEC,fullRHS)
#endif
           call lsmpi_poke()

           call time_start_phase(PHASE_COMM, at = time_intloop_int)
#ifdef VAR_MPI
           if( Ccmodel>MODEL_CC2 )then
              if( alloc_in_dummy )then
                 if(scheme==2.or.scheme==3)then
                    call lsmpi_win_flush(gvvooa%wi(1),local=.true.)
                    call lsmpi_win_flush(gvoova%wi(1),local=.true.)
                 endif
              else
                 if((scheme==2.or.scheme==3) .and. lock_outside)then
                    call tensor_unlock_wins(gvvooa,.true.)
                    call tensor_unlock_wins(gvoova,.true.)
                 endif
              endif
           endif
#endif
           call time_start_phase(PHASE_WORK, at = time_intloop_comm)

           if( Ccmodel > MODEL_CC2 )then
              if(fa<=fg+lg-1)then
                 !CHECK WHETHER THE TERM HAS TO BE DONE AT ALL, i.e. when the first
                 !element in the alpha batch has a smaller index as the last element in
                 !the gamma batch, chose the trafolength as minimum of alpha batch-length
                 !and the difference between first element of alpha batch and last element
                 !of gamma batch
                 call get_a22_and_prepb22_terms_ex(w0%d,w1%d,w2%d,w3%d,tpl%d,tmi%d,no,nv,nb,fa,fg,la,lg,&
                    &xo,yo,xv,yv,omega2,sio4,scheme,[w0%n,w1%n,w2%n,w3%n],lock_outside,&
                    &time_intloop_B1work, time_intloop_B1comm, scal=0.5E0_realk  )


                 !start a new timing phase after these terms
                 call time_start_phase(PHASE_WORK)

              endif
           endif


           !(w0%d):I[ delta gamma alpha beta] <- (w1%d):I[ alpha beta gamma delta ]
           call array_reorder_4d(1.0E0_realk,w1%d,la,nb,lg,nb,[2,3,1,4],0.0E0_realk,w0%d)
           call lsmpi_poke()
           ! (w3%d):I[i gamma alpha beta] = Lambda^h[delta i] I[delta gamma alpha beta]
           call dgemm('t','n',no,lg*la*nb,nb,1.0E0_realk,yo,nb,w0%d,nb,0.0E0_realk,w2%d,no)
           call lsmpi_poke()
           ! (w0%d):I[i gamma alpha j] = (w3%d):I[i gamma alpha beta] Lambda^h[beta j]
           call dgemm('n','n',no*lg*la,no,nb,1.0E0_realk,w2%d,no*lg*la,yo,nb,0.0E0_realk,w0%d,no*lg*la)
           call lsmpi_poke()
           if( Ccmodel > MODEL_CC2 )then
              select case(scheme)
              case(4,3)
                 ! (w3%d):I[alpha gamma i j] <- (w0%d):I[i gamma alpha j]
                 call add_int_to_sio4(w0%d,w2%d,w3%d,nor,no,nv,nb,fa,fg,la,lg,xo,sio4%elm1)
              case(2)
#ifdef VAR_MPI
                 ! (w3):I[ gamma i j alpha] <- (w0):I[i gamma alpha  j]
                 call array_reorder_4d(1.0E0_realk,w0%d,no,lg,la,no,[2,1,4,3],0.0E0_realk,w3%d)
                 ! (w2):I[ l i j alpha] <- (w3):Lambda^p [gamma l ]^T I[gamma i j alpha]
                 call dgemm('t','n',no,no*no*la,lg,1.0E0_realk,xo(fg),nb,w3%d,lg,0.0E0_realk,w2%d,no)
                 ! (w3):I[ k l i j] <- (w2):Lambda^p [alpha l ]^T I[ l i j , alpha]^T
                 call dgemm('t','t',no,no*no*no,la,1.0E0_realk,xo(fa),nb,w2%d,no*no*no,0.0E0_realk,w3%d,no)
                 if(.not.alloc_in_dummy)call tensor_lock_wins(sio4,'s',mode)
                 call tensor_add(sio4,1.0E0_realk,w3%d,order = [1,2,3,4],wrk=w2%d,iwrk=w2%n)
                 if( alloc_in_dummy )then
                    call lsmpi_win_flush(sio4%wi(1),local=.true.)
                 else
                    call tensor_unlock_wins(sio4,.true.)
                 endif
#endif
              end select
           endif
           call lsmpi_poke()


           ! (w2%d):I[gamma i j alpha] <- (w0%d):I[i gamma alpha j]
           call array_reorder_4d(1.0E0_realk,w0%d,no,lg,la,no,[2,1,4,3],0.0E0_realk,w2%d)
           call lsmpi_poke()
           ! (w3%d):I[b i j alpha] = Lamda^p[gamma b] (w2%d):I[gamma i j alpha]
           call dgemm('t','n',nv,no2*la,lg,1.0E0_realk,xv(fg),nb,w2%d,lg,0.0E0_realk,w3%d,nv)
           call lsmpi_poke()
           ! Omega += Lambda^p[alpha a]^T (w3%d):I[b i j alpha]^T
           if(scheme==2)then
#ifdef VAR_MPI
              call time_start_phase(PHASE_COMM, at = time_intloop_work )
              if( lock_outside .and. .not. alloc_in_dummy )call tensor_lock_wins(omega2,'s',mode)
              call time_start_phase(PHASE_WORK, at = time_intloop_comm )
              call dgemm('t','t',nv,o2v,la,0.5E0_realk,xv(fa),nb,w3%d,o2v,0.0E0_realk,w2%d,nv)
              call time_start_phase(PHASE_COMM, at = time_intloop_work )
              call tensor_add(omega2,1.0E0_realk,w2%d,wrk=w0%d,iwrk=w0%n)
              call time_start_phase(PHASE_WORK, at = time_intloop_comm )
#endif
           else
              call dgemm('t','t',nv,o2v,la,0.5E0_realk,xv(fa),nb,w3%d,o2v,1.0E0_realk,omega2%elm1,nv)
           endif
           call lsmpi_poke()

        end do BatchAlpha
     end do BatchGamma


     ! Free integral stuff
     ! *******************
#ifdef VAR_ICHOR
     call FREE_SCREEN_ICHORERI()
     call mem_dealloc(AOGammabatchinfo)
     call mem_dealloc(AOAlphabatchinfo)
#else
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
#endif
     ! free arrays only needed in the batched loops
#ifdef VAR_MPI
     call time_start_phase(PHASE_COMM, at = time_intloop_work )
     if(lock_outside.and.scheme==2)then
        call tensor_unlock_wins(omega2, all_nodes = alloc_in_dummy, check =.not.alloc_in_dummy)
     endif

     if(alloc_in_dummy)then
        if((scheme==3.or.scheme==2))then
           !SWITCH ONE SIDED EPOCH FROM ACCUMULATES TO GETS
           call tensor_unlock_wins(gvvooa, all_nodes = .true.)
           call tensor_unlock_wins(gvoova, all_nodes = .true.)
        endif
        if(scheme==2)then
           call tensor_unlock_wins(sio4, all_nodes = .true.)
        endif
     endif
     call time_start_phase(PHASE_WORK, at = time_intloop_comm )
#endif


     call mem_dealloc(uigcj)
     call mem_dealloc(tpl)
     call mem_dealloc(tmi)
     ! free working matrices and adapt to new requirements
     call mem_dealloc(w0)
     call mem_dealloc(w1)
     call mem_dealloc(w2)
     call mem_dealloc(w3)

#ifdef VAR_MPI
     maxts=0
     nbuff=0
     if((scheme==4.or.scheme==3).and.(ccmodel > MODEL_CC2).and.(.not.local))then
        maxts = max(govov%tsize,maxts)
        nbuff=nbuff+1
     endif
     if(scheme==3)then
        maxts = max(gvoova%tsize,maxts)
        maxts = max(gvvooa%tsize,maxts)
        nbuff=nbuff+1
     endif

     if(nbuff /= 0)then
        buf_size = max(1,min((int((MemFree * 1024.0E0_realk**3)/(8.0E0_realk)) - nbuff * o2v2) / maxts,5))*maxts
     else
        buf_size=0
     endif

     if((scheme==4.or.scheme==3).and.(ccmodel > MODEL_CC2).and.(.not.local))then

        call mem_alloc(buf1,buf_size)

        if(lock_outside)then
           call tensor_lock_wins( govov, 's', mode, all_nodes = alloc_in_dummy )
        endif

        call memory_allocate_tensor_dense( govov )

        call time_start_phase(PHASE_COMM, twall = time_Bcnd )
        call tensor_gather(1.0E0_realk,govov,0.0E0_realk,govov%elm1,o2v2,wrk=buf1,iwrk=buf_size)
        call time_start_phase(PHASE_WORK, twall = time_Bcnd )

        if(.not.lock_outside)then
           call mem_dealloc(buf1)
        endif

     endif

     call get_currently_available_memory(MemFree)

     ! Finish the MPI part of the Residual calculation
     call time_start_phase(PHASE_IDLE, at = time_intloop_work )

     !!!!!!!!!!!!!!!!!!!!!!!!!DO NOT TOUCH THIS BARRIER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call lsmpi_barrier(infpar%lg_comm)
     !!!!!!!!!!!!!!!!!!!!!!!!!DO NOT TOUCH THIS BARRIER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call time_start_phase(PHASE_COMM, at = time_intloop_idle, twall = commtime )

     if(alloc_in_dummy.and.scheme==2)then
        call tensor_lock_wins(omega2,'s',all_nodes=.true.)
     endif

     call time_phases_get_diff(current_wt=phase_counters_int_dir)
     call lsmpi_local_reduction(phase_counters_int_dir,nphases,infpar%master)

     max_wait_time = time_intloop_idle
     min_wait_time = time_intloop_idle


     if(scheme==3)then

        if(lock_outside)then
           call tensor_lock_wins( gvoova , 's', mode, all_nodes = alloc_in_dummy )
           call tensor_lock_wins( gvvooa , 's', mode, all_nodes = alloc_in_dummy )
        endif

        call mem_alloc( gvvoo, o2v2, simple=.true. )
        call mem_alloc( gvoov, o2v2, simple=.true. )

        call mem_alloc(buf2,buf_size)
        call mem_alloc(buf3,buf_size)


        if( Ccmodel>MODEL_CC2 )then
           call tensor_gather(1.0E0_realk,gvoova,0.0E0_realk,gvoov%d,o2v2,wrk=buf2,iwrk=buf_size)
           call tensor_gather(1.0E0_realk,gvvooa,0.0E0_realk,gvvoo%d,o2v2,wrk=buf3,iwrk=buf_size)
        endif

        if(.not.lock_outside)then
           call mem_dealloc(buf2)
           call mem_dealloc(buf3)
        endif

     endif

#endif


#ifdef VAR_MPI

     !GET TIMING INFORMATION
#ifdef VAR_LSDEBUG
     if(print_debug)write(*,'("--rank",I2,", load: ",I5,", w-time:",f15.4)') &
        &infpar%mynum,myload,time_intloop_idle
#endif

     call lsmpi_local_reduction(time_intloop_idle,infpar%master)
     call lsmpi_reduce_realk_max(max_wait_time,infpar%master,infpar%lg_comm)
     call lsmpi_reduce_realk_min(min_wait_time,infpar%master,infpar%lg_comm)

     ave_wait_time = time_intloop_idle/(infpar%nodtot*1.0E0_realk)
     if(master.and.print_debug)then
        write(*,'("----------------------------------------------------------")')
        write(*,'("sum: ",f15.4," 0: ",f15.4," Max: ",f15.4)') time_intloop_idle,&
           &ave_wait_time,max_wait_time
     endif

     startt=MPI_wtime()

     if(infpar%lg_nodtot>1.or.scheme==3) then

        ! The following block is structured like this due to performance reasons
        !***********************************************************************
        if(Ccmodel > MODEL_CC2)then

           if(scheme /= 2) call lsmpi_allreduce(sio4%elm1,int((i8*nor)*no2,kind=8),infpar%lg_comm)

           if(scheme==4)then

              call lsmpi_allreduce(gvvooa%elm1,o2v2,infpar%lg_comm)
              call lsmpi_allreduce(gvoova%elm1,o2v2,infpar%lg_comm)

           endif

        endif

     end if


     if(.not.dynamic_load)then
        call mem_dealloc(tasks)
     else
#ifdef VAR_HAVE_MPI3
        call lsmpi_win_unlock_all(tasksw)
#endif
        call lsmpi_win_free(tasksw)
        call mem_dealloc(tasks,tasksc)
     endif

     if(master.and.DECinfo%PL>2)then
        write(*,'("CCSD time in lsmpi_win_unlock phase B",g10.3)') time_lsmpi_win_unlock - unlock_time
        write(*,'("CCSD time in lsmpi_wait       phase B",g10.3)') time_lsmpi_wait       - waiting_time
        write(*,'("CCSD time in lsmpi_win_flush  phase B",g10.3)') time_lsmpi_win_flush  - flushing_time
     endif
     unlock_time   = time_lsmpi_win_unlock 
     waiting_time  = time_lsmpi_wait
     flushing_time = time_lsmpi_win_flush
#endif

     call time_start_phase(PHASE_WORK, at = time_intloop_comm , twall = time_intloop_stop)

     ! Print timings for the first part
     !*********************************

     commtime     = time_intloop_stop - commtime
     time_intloop = time_intloop_stop - time_intloop
     if(master.and.DECinfo%PL>2)then
        write( *,'("-------------------------------------------------------------------")')
        write( *,'("CCSD total int loop    :",g10.3,"s, comm:",g10.3,"s in collective")')time_intloop,commtime
        write( *,'("     total work time   :",g10.3,"s")') time_intloop_work+time_intloop_B1work
        write( *,'("     total comm time   :",g10.3,"s")') time_intloop_comm+time_intloop_B1comm
        write( *,'("     total ints time   :",g10.3,"s")') time_intloop_int
        write( *,'("     ave work time     :",g10.3,"s")') phase_counters_int_dir(PHASE_WORK_IDX)/float(lg_nnod)
        write( *,'("     ave comm time     :",g10.3,"s")') phase_counters_int_dir(PHASE_COMM_IDX)/float(lg_nnod)
        write( *,'("     ave idle time     :",g10.3,"s")') phase_counters_int_dir(PHASE_IDLE_IDX)/float(lg_nnod)
        write( *,'("     max/ave/min idle  :",g10.3,"s",g10.3,"s",g10.3,"s")') max_wait_time,ave_wait_time,min_wait_time
        write( *,'("     B1 work time      :",g10.3,"s")') time_intloop_B1work
        write( *,'("     B1 comm time      :",g10.3,"s")') time_intloop_B1comm
        if(DECinfo%PL>3)then
           write( *,'("     total time diff   :",g10.3,"s")') time_intloop -&
              & (time_intloop_work+time_intloop_B1work + &
              & time_intloop_comm+time_intloop_B1comm + time_intloop_int + &
              & time_intloop_idle)
        endif
     endif


     ! Reallocate 1 temporary array
     maxsize64 = max(int((i8*nv2)*no2,kind=8),int(nb2,kind=8))
     maxsize64 = max(maxsize64,int((i8*nv2)*nor,kind=8))
     call mem_alloc(w1,maxsize64,simple=.true.)


     call ccsd_debug_print(ccmodel,1,master,local,scheme,print_debug,o2v2,w1,&
        &omega2,govov,gvvooa,gvoova)


     if(Ccmodel>MODEL_CC2)then

        call time_start_phase(PHASE_WORK, twall = time_Bcnd )

        !get B2.2 contributions
        !**********************
        call get_B22_contrib_mo(sio4,t2,w1%d,w2%d,no,nv,omega2,scheme,lock_outside,&
           &time_Bcnd_work,time_Bcnd_comm)


        call tensor_free(sio4)

        call ccsd_debug_print(ccmodel,2,master,local,scheme,print_debug,o2v2,w1,&
           &omega2,govov,gvvooa,gvoova)

#ifdef VAR_MPI
        call time_start_phase(PHASE_COMM, at = time_Bcnd_work )
        if((scheme==4.or.scheme==3).and..not.local)then

           call tensor_unlock_wins(govov, all_nodes = alloc_in_dummy, check =.not.alloc_in_dummy )
           if(lock_outside)call mem_dealloc(buf1)

        endif


        if(scheme==3)then

           if(lock_outside)then
              call tensor_unlock_wins(gvoova, all_nodes = alloc_in_dummy, check =.not.alloc_in_dummy )
              call tensor_unlock_wins(gvvooa, all_nodes = alloc_in_dummy, check =.not.alloc_in_dummy)
              call mem_dealloc(buf2)
              call mem_dealloc(buf3)
           endif
           gvoova%elm1 => gvoov%d
           gvvooa%elm1 => gvvoo%d

        endif

        call time_start_phase(PHASE_WORK, at = time_Bcnd_comm )
#endif
        call time_start_phase(PHASE_WORK, twall = time_cndonly )
        time_Bonly = time_cndonly - time_Bcnd


        !Get the C2 and D2 terms
        !***********************
        if(scheme==4.or.scheme==3)then
           call get_cnd_terms_mo_3n4(w1%d,w2%d,w3%d,t2,u2,govov,gvoova,gvvooa,no,nv,omega2,&
              &scheme,lock_outside,els2add,time_cnd_work,time_cnd_comm)
        else if(scheme==2)then
           call get_cnd_terms_mo_2(w1%d,w2%d,w3%d,t2,u2,govov,gvoova,gvvooa,no,nv,omega2,&
              &scheme,lock_outside,local)
        endif

        call ccsd_debug_print(ccmodel,3,master,local,scheme,print_debug,o2v2,w1,&
           &omega2,govov,gvvooa,gvoova)

        !DEALLOCATE STUFF
        if(scheme==4)then
           call tensor_free(gvoova)
           call tensor_free(gvvooa)
#ifdef VAR_MPI
        else if(scheme==3)then
           gvvooa%elm1 => null()
           gvoova%elm1 => null()
           call tensor_free(gvoova)
           call tensor_free(gvvooa)
           call mem_dealloc(gvoov)
           call mem_dealloc(gvvoo)
        else if(scheme==2)then
           call tensor_free(gvoova)
           call tensor_free(gvvooa)
           if(alloc_in_dummy)then
              call tensor_unlock_wins(omega2,all_nodes=.true.)
              call tensor_unlock_wins(u2,all_nodes=.true.)
           endif
#endif
        endif

        call time_start_phase(PHASE_WORK, ttot = time_cndonly )
        time_Bcnd      = time_Bonly     + time_cndonly
        time_Bcnd_work = time_Bcnd_work + time_cnd_work
        time_Bcnd_comm = time_Bcnd_comm + time_cnd_comm


        !TIMING INFO
        if(master.and.DECinfo%PL>2)then
           write( *,'("-------------------------------------------------------------------")')
           write( *,'("CCSD total B+cnd       :",g10.3,"s")')time_Bcnd
           write( *,'("     work time         :",g10.3,"s")') time_Bcnd_work
           write( *,'("     comm time         :",g10.3,"s")') time_Bcnd_comm
           if(DECinfo%PL>3)then
              write( *,'("     total B           :",g10.3,"s")') time_Bonly
              write( *,'("     total cnd         :",g10.3,"s")') time_cndonly
           endif
        endif

     endif


     !IN CASE OF MPI (AND CORRECT SCHEME) REDUCE TO MASTER
     !*****************************************************
#ifdef VAR_MPI
     call time_start_phase(PHASE_COMM)

     if(infpar%lg_nodtot>1) then
        if(scheme==4.or.scheme==3)&
           &call lsmpi_local_reduction(omega2%elm1,o2v2,infpar%master)
        call lsmpi_local_reduction(Gbi,nb*no,infpar%master)
        call lsmpi_local_reduction(Had,nb*nv,infpar%master)
     endif

     call time_start_phase(PHASE_WORK, dt = time_reduction2)

     !convert stuff
     !set for correct access again, save as i a j b
     if(.not.local)then
        if((master.and..not.(scheme==2)).or.scheme==3) call tensor_deallocate_dense(govov)
        govov%itype      = TT_TILED_DIST
     endif

     govov%access_type  = AT_MASTER_ACCESS
     omega2%access_type = AT_MASTER_ACCESS
     t2%access_type     = AT_MASTER_ACCESS
     if(scheme==2)then
        u2%access_type     = AT_MASTER_ACCESS
     endif
#endif


     !call get_currently_available_memory(MemFree2)
     !call get_available_memory(6,MemFree3,memfound,.true.)
     !print *,infpar%lg_mynum,"slaves return",MemFree2,MemFree3
     !call lsmpi_barrier(infpar%lg_comm)

     ! slaves should exit the subroutine after the main work is done
     if(.not. master) then
        call mem_dealloc(w1)
        call mem_dealloc(Had)
        call mem_dealloc(Gbi)
        if(scheme==4.or.scheme==3)call tensor_free(u2)
        if(scheme==4)then
           call memory_deallocate_tensor_dense(govov)
        endif

        return

     endif

#ifdef VAR_LSDEBUG
     if(print_debug)then
        call print_norm(Gbi,int((i8*no)*nb,kind=8)," NORM(Gbi)       :")
        call print_norm(Had,int((i8*nv)*nb,kind=8)," NORM(Had)       :")
        call print_norm(omega2,                    " NORM(omega2 s-o):")
        call print_norm(govov,                     " NORM(govov s-o) :")
     endif
#endif

     call time_start_phase(PHASE_WORK, twall = time_get_ao_fock)

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
     !    IF(DECinfo%DFTreference)THEN
     !       call II_get_xc_fock_mat_full(DECinfo%output,DECinfo%output,&
     !           ....)
     !    ENDIF

     !use dens as temporay array 


     call ii_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsItem%setting,&
        & Dens%elms,nb,nb,AORdefault,AORdefault)
     ! Add one- and two-electron contributions to Fock matrix
     call daxpy(nb2,1.0E0_realk,Dens%elms,1,iFock%elms,1)
     !Free the density matrix
     call mat_free(Dens)

     ! KK: Add long-range Fock correction
     call daxpy(nb2,1.0E0_realk,deltafock,1,iFock%elms,1)

     call time_start_phase(PHASE_WORK, ttot = time_get_ao_fock, twall = time_get_mo_fock)

     if(print_debug)then
        call print_norm(deltafock,int((i8*nb)*nb,kind=8), " NORM(deltafock):")
        call print_norm(iFock%elms,int((i8*nb)*nb,kind=8)," NORM(iFock)    :")
     endif



     !Transform inactive Fock matrix into the different mo subspaces
     if (Ccmodel>MODEL_CC2) then
        ! -> Foo
        call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,iFock%elms,nb,0.0E0_realk,w1%d,no)
        call dgemm('n','n',no,no,nb,1.0E0_realk,w1%d,no,yo,nb,0.0E0_realk,ppfock,no)
        ! -> Fov
        call dgemm('n','n',no,nv,nb,1.0E0_realk,w1%d,no,yv,nb,0.0E0_realk,pqfock,no)
        ! -> Fvo
        call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,iFock%elms,nb,0.0E0_realk,w1%d,nv)
        call dgemm('n','n',nv,no,nb,1.0E0_realk,w1%d,nv,yo,nb,0.0E0_realk,qpfock,nv)
        ! -> Fvv
        call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1%d,nv,yv,nb,0.0E0_realk,qqfock,nv)
     else
        ! -> Foo
        call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,fock,nb,0.0E0_realk,w1%d,no)
        call dgemm('n','n',no,no,nb,1.0E0_realk,w1%d,no,yo,nb,0.0E0_realk,ppfock,no)
        ! -> Fov
        call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,iFock%elms,nb,0.0E0_realk,w1%d,no)
        call dgemm('n','n',no,nv,nb,1.0E0_realk,w1%d,no,yv,nb,0.0E0_realk,pqfock,no)
        ! -> Fvo
        call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,iFock%elms,nb,0.0E0_realk,w1%d,nv)
        call dgemm('n','n',nv,no,nb,1.0E0_realk,w1%d,nv,yo,nb,0.0E0_realk,qpfock,nv)
        ! -> Fvv
        call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,fock,nb,0.0E0_realk,w1%d,nv)
        call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1%d,nv,yv,nb,0.0E0_realk,qqfock,nv)
     endif



     if(print_debug)then
        call print_norm(ppfock,int((i8*no)*no,kind=8)," NORM(ppfock)   :")
        call print_norm(pqfock,int((i8*no)*nv,kind=8)," NORM(pqfock)   :")
        call print_norm(qpfock,int((i8*no)*nv,kind=8)," NORM(qpfock)   :")
        call print_norm(qqfock,int((i8*nv)*nv,kind=8)," NORM(qqfock)   :")
     endif

     !Free the AO fock matrix
     call mat_free(iFock)

     call time_start_phase(PHASE_WORK, ttot = time_get_mo_fock, twall = time_Esing)

     if(master.and.DECinfo%PL>2)then
        write( *,'("-------------------------------------------------------------------")')
        write( *,'("CCSD total fock        :",g10.3,"s")') time_get_ao_fock + time_get_mo_fock
        if(DECinfo%PL>3)then
           write( *,'("     total AO only     :",g10.3,"s")') time_get_ao_fock
           write( *,'("     total MO only     :",g10.3,"s")') time_get_mo_fock
        endif
     endif



     !CCD can be achieved by not using singles residual updates here
     if(.not. DECinfo%CCDhack)then

#ifdef VAR_MPI
        if(scheme==2)then
           call tensor_lock_wins(u2,'s', mode, all_nodes = alloc_in_dummy )
           call mem_alloc(w2,u2%tsize * 3_long)
        endif
#endif
        !calculate singles I term
        ! Reorder u [c a i k] -> u [a i c k]
        if(scheme==4.or.scheme==3)then
           call array_reorder_4d(1.0E0_realk,u2%elm1,nv,nv,no,no,[2,3,4,1],0.0E0_realk,w1%d)
        else if(scheme==2)then
           call tensor_convert(u2,w1%d,[2,3,4,1],wrk=w2%d,iwrk=w2%n)
        endif


        !GET SINGLES CONTRIBUTIONS
        !*************************

        !calculate singles J term
        ! F [a i] = Omega [a i]
        call dcopy(no*nv,qpfock,1,omega1,1)

#ifdef VAR_MPI
        if(scheme==2)then
           call tensor_unlock_wins(u2, check=.not.alloc_in_dummy, all_nodes = alloc_in_dummy )
           call mem_dealloc(w2)
        endif
#endif

        ! u [a i k c] * F[k c] =+ Omega [a i]
        call dgemv('n',nv*no,nv*no,1.0E0_realk,w1%d,nv*no,pqfock,1,1.0E0_realk,omega1,1)

        !calculate singles G term
        ! Lambda^p [alpha a]^T Gbi [alpha i] =+ Omega [a i]
        call dgemm('t','n',nv,no,nb,1.0E0_realk,xv,nb,Gbi,nb,1.0E0_realk,omega1,nv)
        !calculate singles H term
        ! (-1) Had [a delta] * Lambda^h [delta i] =+ Omega[a i]
        call dgemm('n','n',nv,no,nb,-1.0E0_realk,Had,nv,yo,nb,1.0E0_realk,omega1,nv)

     endif


     !GET DOUBLES E2 TERM - AND INTRODUCE PERMUTATIONAL SYMMMETRY
     !***********************************************************
     call calculate_E2_and_permute(ccmodel,ppfock,qqfock,w1%d,t2,xo,yv,Gbi,Had,no,nv,nb,&
        &omega2,o2v2,scheme,print_debug,lock_outside,time_Esing_work,time_Esing_comm)

     call mem_dealloc(Had)
     call mem_dealloc(Gbi)

     call time_start_phase(PHASE_WORK, ttot = time_Esing)

     if(master.and.DECinfo%PL>2)then
        write( *,'("-------------------------------------------------------------------")')
        write( *,'("CCSD total E2+singles  :",g10.3,"s")') time_Esing
        if(DECinfo%PL>3)then
           write( *,'("     comm              :",g10.3,"s")') time_Esing_comm
           write( *,'("     work              :",g10.3,"s")') time_Esing_work
        endif
     endif

#ifdef VAR_MPI
     if((.not.local).and.(scheme==4.or.scheme==3))then
        call tensor_mv_dense2tiled(omega2,.true.)
        call tensor_mv_dense2tiled(t2,.true.)
     endif
     if(master.and.DECinfo%PL>2)then
        write(*,'("CCSD time in lsmpi_win_unlock phase C",g10.3)') time_lsmpi_win_unlock - unlock_time
        write(*,'("CCSD time in lsmpi_wait       phase C",g10.3)') time_lsmpi_wait       - waiting_time
        write(*,'("CCSD time in lsmpi_win_flush  phase C",g10.3)') time_lsmpi_win_flush  - flushing_time
     endif
#endif

     call mem_dealloc(w1)
     call tensor_free(u2)

     call time_start_phase(PHASE_WORK, ttot = twall)
     if(master.and.DECinfo%PL>2)then
        write( *,'("-------------------------------------------------------------------")')
        write( *,'("CCSD residual total    ",g10.3,"s")') twall
     endif

     call LSTIMER('START',tcpu_end,twall_end,DECinfo%output)


     if(print_debug)then
        call print_norm(omega1,int((i8*no)*nv,kind=8)," NORM(omega1):")
        call print_norm(omega2,                       " NORM(omega2):")
     endif

  end subroutine get_ccsd_residual_integral_driven

  subroutine yet_another_ccsd_residual(ccmodel,deltafock,omega2,t2,fock,govov,no,nv,&
        ppfock,qqfock,pqfock,qpfock,xo_f,xv_f,yo_f,yv_f,nb,MyLsItem, omega1,iter,local,rest)
     implicit none

     !> CC model
     integer,intent(in) :: ccmodel
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
     type(tensor),intent(inout) :: omega2
     !> the current guess amplitudes
     !real(realk),pointer,intent(in) ::t2(:)
     type(tensor),intent(inout) :: t2
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
     real(realk),pointer :: xo_f(:),xv_f(:),yo_f(:),yv_f(:)

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
     type(tensor),intent(inout) :: govov
     logical, intent(in) :: local
     !logical that specifies whether the amplitudes were read
     logical, optional, intent(inout) :: rest

     ! elementary types needed for the calculation
     type(mpi_realk)      :: w0,w1,w2,w3
     real(realk), pointer :: t2_d(:,:,:,:)
     type(c_ptr) :: Hadc,t2_dc, Gbic
     integer(kind=ls_mpik) :: Hadw,t2_dw,Gbiw,gvvoow,gvoovw

     integer(kind=8) :: w0size,w1size,w2size,w3size,neloc
     integer :: nors, nvrs, order4(4), order2(2)

     ! Variables for mpi
     logical :: master,lg_master,parent
     integer :: fintel,nintel,fe,ne
     real(realk) :: startt, stopp
     integer(kind=ls_mpik) :: ierr

     integer :: sio4_mode, sio4_dims(4),sio4_tdim(4) 
     type(tensor) :: u2, sio4
     type(tensor) :: gvoov,gvvoo
     type(tensor) :: Cint, Eint
     type(tensor) :: int1, int2, int3, int4
     type(tensor) :: uigcj,tpl,tmi
     type(tensor) :: Had, Gbi, Gbi_local
     type(tensor) :: xo, xo_fa, xo_fg
     type(tensor) :: yo, yo_fa, yo_fg
     type(tensor) :: xv, xv_fa, xv_fg
     type(tensor) :: yv, yv_fa, yv_fg
     !special arrays for scheme=1
     type(tensor) :: t2jabi,u2kcjb
     integer,pointer       :: tasks(:)
     type(c_ptr)           :: tasksc
     integer(kind=ls_mpik) :: tasksw,taskslw
     integer(kind=ls_mpik) :: me,nnod
     integer(kind=8)       :: len81,len82
     integer(kind=4)       :: len41,len42
     integer               :: lenI1,lenI2
     integer,parameter :: inflen=5
     real(realk)       :: inf(inflen)
#ifdef VAR_MPI
     ! stuff for direct communication
     integer(kind=ls_mpik) :: gvvoo_w, gvoov_w
     integer(kind=ls_mpik) :: hstatus, nctr,mode
     integer :: rcnt(infpar%lg_nodtot),dsp(infpar%lg_nodtot)
     character*(MPI_MAX_PROCESSOR_NAME) :: hname
     real(realk),pointer :: buf1(:), buf2(:), buf3(:),p2(:,:)
     !integer(kind=ls_mpik),pointer :: tasksw(:)
     integer(kind=8) :: maxts,nbuff
#endif
     logical :: lock_outside

     type(tensor) :: iFock
     ! CHECKING and MEASURING variables
     integer(kind=long) :: maxsize64,dummy64
     integer :: myload,nelms,n4
     real(realk) :: tcpu, twall,tcpu1,twall1,tcpu2,twall2, deb1,deb2,ActuallyUsed
     real(realk) :: MemFree,MemFreeMin,MemUsed
     real(realk) :: tcpu_end,twall_end,time_a, time_c, time_d,time_singles
     real(realk) :: time_doubles,timewall_start,wait_time,max_wait_time,min_wait_time,ave_wait_time
     integer     :: scheme
     integer(kind=8) :: els2add
     logical :: memfound

     ! variables used for BATCH construction and INTEGRAL calculation
     logical :: save_cs_screen, save_ps_screen
     integer :: alphaB,gammaB,dimAlpha,dimGamma
     integer :: dim1,dim2,dim3,K,MinAObatch
     integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
     integer :: iorb,nthreads,idx
#ifdef VAR_ICHOR
     type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
     type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
     integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
     logical :: MoTrans, NoSymmetry,SameMol
#else
     Character(80)        :: FilenameCS,FilenamePS
     Character(80),pointer:: BatchfilenamesCS(:,:)
     Character(80),pointer:: BatchfilenamesPS(:,:)
     logical :: FoundInMem,doscreen
     integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
     integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:)
     TYPE(DECscreenITEM)    :: DecScreen
#endif
     integer, pointer :: batchdimAlpha(:), batchdimGamma(:)     
     type(batchtoorb), pointer :: batch2orbAlpha(:)
     type(batchtoorb), pointer :: batch2orbGamma(:)
     Character            :: INTSPEC(5)
     logical :: fullRHS
     integer :: MaxActualDimAlpha,nbatchesAlpha
     integer :: MaxActualDimGamma,nbatchesGamma
     integer :: ndimA,ndimB,ndimC,ndimD
     integer :: ndimAs,ndimBs,ndimCs,ndimDs
     integer :: startA,startB,startC,startD

     integer :: a,b,i,j,l,m,n,c,d,fa,fg,la,lg,worksize
     integer :: nb2,nb3,nv2,no2,b2v,o2v,v2o,no3,vs,os,bs,as,gs
     integer(kind=8) :: nb4,o2v2,no4,buf_size
     integer :: tlen,tred,nor,nvr,goffs,aoffs
     integer :: prev_alphaB,mpi_buf,ccmodel_copy
     logical :: jobtodo,first_round,dynamic_load,restart,print_debug
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !TEST AND DEVELOPMENT VARIABLES!!!!!
     real(realk) :: op_start,op_stop, dt_last_phase
     real(realk) :: time_init,    time_init_work,    time_init_comm, time_cndonly
     real(realk) :: time_intloop, time_intloop_work, time_intloop_comm, time_intloop_idle, time_intloop_int
     real(realk) :: unlock_time_max  ,unlock_time_min  , unlock_time 
     real(realk) :: waiting_time_max ,waiting_time_min , waiting_time
     real(realk) :: flushing_time_max,flushing_time_min, flushing_time
     real(realk) :: phase_cntrs(nphases)
     real(realk) :: tot_intloop, tot_intloop_max, tot_intloop_min
     real(realk) :: time_yv_fg_tot_max,   time_yv_fg_tot_min,   time_yv_fg_tot,    time_yvfg
     real(realk) :: time_yo_fg_tot_max,   time_yo_fg_tot_min,   time_yo_fg_tot,    time_yofg 
     real(realk) :: time_xo_fg_tot_max,   time_xo_fg_tot_min,   time_xo_fg_tot,    time_xofg 
     real(realk) :: time_xv_fg_tot_max,   time_xv_fg_tot_min,   time_xv_fg_tot,    time_xvfg 
     real(realk) :: time_uigcj_tot_max,   time_uigcj_tot_min,   time_uigcj_tot,    time_uigcj
     real(realk) :: time_xo_fa_tot_max,   time_xo_fa_tot_min,   time_xo_fa_tot,    time_xofa
     real(realk) :: time_xv_fa_tot_max,   time_xv_fa_tot_min,   time_xv_fa_tot,    time_xvfa 
     real(realk) :: time_yv_fa_tot_max,   time_yv_fa_tot_min,   time_yv_fa_tot,    time_yvfa
     real(realk) :: time_cont1_tot_max,   time_cont1_tot_min,   time_cont1_tot,    time_cont1 
     real(realk) :: time_int1_tot_max,    time_int1_tot_min,    time_int1_tot ,    time_int1
     real(realk) :: time_cont2_tot_max,   time_cont2_tot_min,   time_cont2_tot,    time_cont2
     real(realk) :: time_cont3_tot_max,   time_cont3_tot_min,   time_cont3_tot,    time_cont3
     real(realk) :: time_cont4_tot_max,   time_cont4_tot_min,   time_cont4_tot,    time_cont4
     real(realk) :: time_cont5_tot_max,   time_cont5_tot_min,   time_cont5_tot,    time_cont5
     real(realk) :: time_cont6_tot_max,   time_cont6_tot_min,   time_cont6_tot,    time_cont6 
     real(realk) :: time_cont7_tot_max,   time_cont7_tot_min,   time_cont7_tot,    time_cont7 
     real(realk) :: time_cont8_tot_max,   time_cont8_tot_min,   time_cont8_tot,    time_cont8 
     real(realk) :: time_convGBI_tot_max, time_convGBI_tot_min, time_convGBI_tot,  time_convGBI
     real(realk) :: time_cont9_tot_max,   time_cont9_tot_min,   time_cont9_tot,    time_cont9
     real(realk) :: time_cont10_tot_max,  time_cont10_tot_min,  time_cont10_tot,   time_cont10
     real(realk) :: time_int2_tot_max,    time_int2_tot_min,    time_int2_tot,     time_int2
     real(realk) :: time_cont11_tot_max,  time_cont11_tot_min,  time_cont11_tot,   time_cont11
     real(realk) :: time_cont12_tot_max,  time_cont12_tot_min,  time_cont12_tot,   time_cont12
     real(realk) :: time_cont13_tot_max,  time_cont13_tot_min,  time_cont13_tot,   time_cont13
     real(realk) :: time_cont14_tot_max,  time_cont14_tot_min,  time_cont14_tot,   time_cont14
     real(realk) :: time_cont15_tot_max,  time_cont15_tot_min,  time_cont15_tot,   time_cont15
     real(realk) :: time_cont16_tot_max,  time_cont16_tot_min,  time_cont16_tot,   time_cont16
     real(realk) :: time_cont17_tot_max,  time_cont17_tot_min,  time_cont17_tot,   time_cont17
     real(realk) :: time_cont18_tot_max,  time_cont18_tot_min,  time_cont18_tot,   time_cont18
     real(realk) :: time_cont19_tot_max,  time_cont19_tot_min,  time_cont19_tot,   time_cont19
     real(realk) :: time_cont20_tot_max,  time_cont20_tot_min,  time_cont20_tot,   time_cont20
     real(realk) :: time_cont21_tot_max,  time_cont21_tot_min,  time_cont21_tot,   time_cont21
     real(realk) :: time_cont22_tot_max,  time_cont22_tot_min,  time_cont22_tot,   time_cont22
     real(realk) :: time_cont23_tot_max,  time_cont23_tot_min,  time_cont23_tot,   time_cont23
     real(realk) :: time_w_min, time_c_min, time_i_min, time_w_max, time_c_max, time_i_max
     integer :: starts(4),residual_nr
     integer(kind=long) :: xyz,zyx1,zyx2
     logical :: debug
     character(tensor_MSG_LEN) :: msg
     real(realk), parameter :: frac = 0.8E0_realk
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef VAR_OMP
     integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
     character(4) :: def_atype

     !init timing variables
     call time_start_phase(PHASE_WORK, twall = twall)
#ifdef VAR_MPI

     time_init_work       = 0.0E0_realk
     time_init_comm       = 0.0E0_realk
     tot_intloop          = 0.0E0_realk
     time_yv_fg_tot       = 0.0E0_realk
     time_yo_fg_tot       = 0.0E0_realk
     time_xo_fg_tot       = 0.0E0_realk
     time_xv_fg_tot       = 0.0E0_realk
     time_uigcj_tot       = 0.0E0_realk
     time_xo_fa_tot       = 0.0E0_realk
     time_xv_fa_tot       = 0.0E0_realk
     time_yv_fa_tot       = 0.0E0_realk
     time_cont1_tot       = 0.0E0_realk
     time_int1_tot        = 0.0E0_realk
     time_cont2_tot       = 0.0E0_realk
     time_cont3_tot       = 0.0E0_realk
     time_cont4_tot       = 0.0E0_realk
     time_cont5_tot       = 0.0E0_realk
     time_cont6_tot       = 0.0E0_realk
     time_cont7_tot       = 0.0E0_realk
     time_cont8_tot       = 0.0E0_realk
     time_convGBI_tot     = 0.0E0_realk
     time_cont9_tot       = 0.0E0_realk
     time_cont10_tot      = 0.0E0_realk
     time_int2_tot        = 0.0E0_realk
     time_cont11_tot      = 0.0E0_realk
     time_cont12_tot      = 0.0E0_realk
     time_cont13_tot      = 0.0E0_realk
     time_cont14_tot      = 0.0E0_realk
     time_cont15_tot      = 0.0E0_realk
     time_cont16_tot      = 0.0E0_realk
     time_cont17_tot      = 0.0E0_realk
     time_cont18_tot      = 0.0E0_realk
     time_cont19_tot      = 0.0E0_realk
     time_cont20_tot      = 0.0E0_realk
     time_cont21_tot      = 0.0E0_realk
     time_cont22_tot      = 0.0E0_realk
     time_cont23_tot      = 0.0E0_realk
     unlock_time          = time_lsmpi_win_unlock 
     waiting_time         = time_lsmpi_wait
     flushing_time        = time_lsmpi_win_flush

     ! Set MPI related info
     ! ********************
     me                       = infpar%lg_mynum
     nnod                     = infpar%lg_nodtot
     lock_outside             = DECinfo%CCSD_NO_DEBUG_COMM
     mode                     = MPI_MODE_NOCHECK
     master                   = (me == 0)


     ! Set default values for the path throug the routine
     ! **************************************************
     restart                  = .false.
     if(present(rest))restart = rest
     scheme                   = 0
     dynamic_load             = DECinfo%dyn_load
     startt                   = 0.0E0_realk
     stopp                    = 0.0E0_realk
     print_debug              = (DECinfo%PL>3)
     debug                    = .false.

     ! Set some shorthand notations
     ! ****************************
     nb2                      = nb*nb
     nb3                      = nb2*nb
     nb4                      = int((i8*nb3)*nb,kind=long)
     nv2                      = nv*nv
     no2                      = no*no
     no3                      = no2*no
     no4                      = int((i8*no3)*no,kind=long)
     b2v                      = nb2*nv
     o2v                      = no2*nv
     v2o                      = nv2*no
     o2v2                     = int((i8*nv2)*no2,kind=long)
     nor                      = no*(no+1)/2
     nvr                      = nv*(nv+1)/2
     vs                       = t2%tdim(1)
     os                       = t2%tdim(3)
     residual_nr              = RN_YET_ANOTHER_RES
     bs                       = vs
     nors                     = get_split_scheme_0(nor)
     nvrs                     = get_split_scheme_0(nvr)


     ! Set integral info
     ! *****************
     !R = Regular Basis set on the 1th center R = Regular Basis set on the 2th center R = Regular Basis set on the 3th center R = Regular Basis set on the 4th center 
     !C = Coulomb operator
     INTSPEC = ['R','R','R','R','C'] 
     save_cs_screen = mylsitem%setting%SCHEME%CS_SCREEN
     save_ps_screen = mylsitem%setting%SCHEME%PS_SCREEN
     mylsitem%setting%SCHEME%CS_SCREEN = .FALSE.
     mylsitem%setting%SCHEME%PS_SCREEN = .FALSE.
     doscreen = mylsitem%setting%SCHEME%CS_SCREEN.OR.mylsitem%setting%SCHEME%PS_SCREEN


     StartUpSlaves: if(master .and. nnod>1) then
        call time_start_phase(PHASE_COMM)
        call ls_mpibcast(CCSDDATA,infpar%master,infpar%lg_comm)
        ccmodel_copy = ccmodel
        call mpi_communicate_ccsd_calcdata(ccmodel_copy,omega2,t2,govov,xo_f,xv_f,&
           &yo_f,yv_f,MyLsItem,nb,nv,no,iter,local,residual_nr)
        call time_start_phase(PHASE_WORK, at = time_init_comm)
     endif StartUpSlaves

     t2%access_type     = AT_ALL_ACCESS
     omega2%access_type = AT_ALL_ACCESS
     govov%access_type  = AT_ALL_ACCESS

     call tensor_ainit(xo,[nb,no],2,local=local,tdims=[bs,os],atype="TDAR")
     call tensor_ainit(yo,[nb,no],2,local=local,tdims=[bs,os],atype="TDAR")
     call tensor_ainit(xv,[nb,nv],2,local=local,tdims=[bs,vs],atype="TDAR")
     call tensor_ainit(yv,[nb,nv],2,local=local,tdims=[bs,vs],atype="TDAR")

     call tensor_convert(xo_f,xo)
     call tensor_convert(yo_f,yo)
     call tensor_convert(xv_f,xv)
     call tensor_convert(yv_f,yv)

     ! Memory info - should synchronize the nodes
     ! ***********
     call get_currently_available_memory(MemFree)
     MemFreeMin = MemFree
     call lsmpi_reduce_realk_min(MemFreeMin,infpar%master,infpar%lg_comm)
     ! Estimate free mem to be the min times the number of nodes
     MemFree = MemFreeMin * nnod

     if(print_debug)then
        call print_norm(xo," NORM(xo)    :",print_on_rank = 0)
        call print_norm(xv," NORM(xv)    :",print_on_rank = 0)
        call print_norm(yo," NORM(yo)    :",print_on_rank = 0)
        call print_norm(yv," NORM(yv)    :",print_on_rank = 0)
     endif

     if(master) then
        !==================================================
        !                  Batch construction             !
        !==================================================
#ifdef VAR_ICHOR
        iAO = 4 
        call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
#else
        call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
#endif

        ! Get free memory and determine maximum batch sizes
        ! -------------------------------------------------
        call get_max_batch_sizes(scheme,nb,bs,nv,vs,no,os,MaxActualDimAlpha,MaxActualDimGamma,&
        &MinAObatch,DECinfo%manual_batchsizes,iter,MemFreeMin,.true.,els2add,local,.false.,mylsitem%setting,intspec)

        if(scheme /= 0 ) call lsquit("ERROR(yet_another_ccsd_residual): for the collective memory only scheme 0 possible",-1)
        
     endif

     !all communication for MPI prior to the loop
     call time_start_phase(PHASE_COMM, at = time_init_work )
     call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
     call ls_mpi_buffer(scheme,infpar%master)
     call ls_mpi_buffer(print_debug,infpar%master)
     call ls_mpi_buffer(dynamic_load,infpar%master)
     call ls_mpi_buffer(restart,infpar%master)
     call ls_mpi_buffer(MaxActualDimAlpha,infpar%master)
     call ls_mpi_buffer(MaxActualDimGamma,infpar%master)
     call ls_mpi_buffer(lock_outside,infpar%master)
     call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
     call time_start_phase(PHASE_WORK, at = time_init_comm)

     nbatchesAlpha = nb / MaxActualDimAlpha
     if( mod( nb, MaxActualDimAlpha ) > 0 ) nbatchesAlpha = nbatchesAlpha + 1
     nbatchesGamma = nb / MaxActualDimGamma
     if( mod( nb, MaxActualDimGamma ) > 0 ) nbatchesGamma = nbatchesGamma + 1


     hstatus = 80
     CALL MPI_GET_PROCESSOR_NAME(hname,hstatus,ierr)

     call tensor_lock_wins(omega2,'s',all_nodes = alloc_in_dummy)
     call tensor_zero(omega2)


     ! ************************************************
     ! *  Allocate matrices used in the batched loop  *
     ! ************************************************

     ! Use the dense amplitudes
     ! ------------------------
     !get the t+ and t- for the Kobayshi-like B2 term
     call tensor_ainit(tpl,[nor,nvr],2,local=local,tdims=[nors,nvrs],atype="TDAR")
     call tensor_ainit(tmi,[nor,nvr],2,local=local,tdims=[nors,nvrs],atype="TDAR")
     
     !call get_tpl_and_tmi(t2,tpl,tmi)
     call tensor_zero(tpl)
     call tensor_zero(tmi)

     if(print_debug)then
        call print_norm(tpl," NORM(tpl)   :",print_on_rank=0)
        call print_norm(tmi," NORM(tmi)   :",print_on_rank=0)
     endif

     call tensor_ainit( u2, [nv,nv,no,no], 4, local=local, atype='TDAR', tdims=[vs,vs,os,os] )
     call tensor_add( u2,  2.0E0_realk, t2, a = 0.0E0_realk, order=[2,1,3,4] )
     call tensor_add( u2, -1.0E0_realk, t2, order=[2,1,4,3] )

     if( alloc_in_dummy )then
        call tensor_lock_wins(u2,'s',all_nodes = .true.)
     endif

     if( print_debug )then
        call print_norm(u2," NORM(u2)    :",print_on_rank=0)
     endif

     print *,"alloc 1"
     call lsmpi_barrier(infpar%lg_comm)
     call tensor_ainit(Had,[nv,nb],2,local=local,tdims=[vs,bs],atype="TDAR")
     print *,"alloc 2"
     call lsmpi_barrier(infpar%lg_comm)
     call tensor_ainit(Gbi,[nb,no],2,local=local,tdims=[bs,os],atype="TDAR")
     print *,"alloc 3"
     call lsmpi_barrier(infpar%lg_comm)
     call tensor_ainit(Gbi_local,[nb,no],2,local=local,atype="REAR")

     call tensor_zero( Had )
     call tensor_zero( Gbi )
     call tensor_zero( Gbi_local )

     if( CCmodel > MODEL_CC2 )then

     print *,"alloc 4"
     call lsmpi_barrier(infpar%lg_comm)
        call tensor_ainit(gvvoo, [nv,no,no,nv],4, local=local, atype="TDAR", tdims=[vs,os,os,vs])
     print *,"alloc 5"
     call lsmpi_barrier(infpar%lg_comm)
        call tensor_ainit(gvoov, [nv,no,nv,no],4, local=local, atype="TDAR", tdims=[vs,os,vs,os])
     print *,"alloc 6"
     call lsmpi_barrier(infpar%lg_comm)
        call tensor_ainit(sio4,  [no,no,no,no],4, local=local, atype="TDAR", tdims=[os,os,os,os])
     print *,"alloc 7"
     call lsmpi_barrier(infpar%lg_comm)

        call tensor_zero(gvvoo)
        call tensor_zero(gvoov)
        call tensor_zero(sio4)


        if( alloc_in_dummy )then
           call tensor_lock_wins(gvvoo,'s',all_nodes = .true.)
           call tensor_lock_wins(gvoov,'s',all_nodes = .true.)
           call tensor_lock_wins(sio4, 's',all_nodes = .true.)
        endif

     endif

     ! allocate working arrays depending on the batch sizes
     w0size = get_wsize_for_ccsd_int_direct(0,no,os,nv,vs,nb,bs,MaxActualDimAlpha,&
        &MaxActualDimGamma,0,setting = mylsitem%setting, intspec = intspec)
     print *,"alloc w0",w0size,nb,MaxActualDimAlpha,MaxActualDimGamma
     call lsmpi_barrier(infpar%lg_comm)
     call mem_alloc( w0, w0size , simple = .true. )

     w1size = get_wsize_for_ccsd_int_direct(1,no,os,nv,vs,nb,bs,MaxActualDimAlpha,&
        &MaxActualDimGamma,0,setting = mylsitem%setting, intspec = intspec)
     print *,"alloc w1",w1size
     call lsmpi_barrier(infpar%lg_comm)
     call mem_alloc( w1, w1size , simple = .true. )

     !first touch
     w0%d(1) = 0.0E0_realk
     w1%d(1) = 0.0E0_realk


#ifdef VAR_OMP
     nthreads=OMP_GET_MAX_THREADS()
     if(master.and.DECinfo%PL>2)write(DECinfo%output,*)&
        & 'Starting CCSD residuals - OMP. Number of threads: ', OMP_GET_MAX_THREADS()
#else
     nthreads=1
     if(master.and.DECinfo%PL>2)write(DECinfo%output,*) &
        &'Starting CCSD integral/amplitudes - NO OMP!'
#endif


     if(master.and.DECinfo%PL>2)then
        write(*,'("CCSD time in lsmpi_win_unlock phase A",g10.3)') time_lsmpi_win_unlock - unlock_time
        write(*,'("CCSD time in lsmpi_wait       phase A",g10.3)') time_lsmpi_wait       - waiting_time
        write(*,'("CCSD time in lsmpi_win_flush  phase A",g10.3)') time_lsmpi_win_flush  - flushing_time
     endif
     unlock_time   = time_lsmpi_win_unlock 
     waiting_time  = time_lsmpi_wait
     flushing_time = time_lsmpi_win_flush

     myload = 0

     !TIMING
     call time_start_phase(PHASE_WORK, at = time_init_work, twall = time_intloop )
     time_init = time_init_work + time_init_comm
     if(master.and.DECinfo%PL>2)then
        write( *,'("-------------------------------------------------------------------")')
        write( *,'("CCSD residual init work:",g10.3,"s, comm:",g10.3,"s")')time_init_work,time_init_comm
     endif
     call time_phases_get_current(current_wt=phase_cntrs)

     fullRHS=(nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)

     !**********************************
     ! Begin the loop over gamma batches
     !**********************************
     call time_start_phase( PHASE_WORK, twall = tot_intloop )

     BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
        !short hand notation
        fg = 1 + (gammaB-1)*MaxActualDimGamma
        lg = nb - fg + 1
        if( lg >= MaxActualDimGamma )then
           lg = MaxActualDimGamma
        endif

        if( lg > bs )then
           gs = bs
        else
           gs = lg
        endif

        call time_start_phase(PHASE_WORK, twall = time_yvfg )
        call tensor_ainit(yv_fg, [lg,nv], 2, local=local, atype="TDAR", tdims=[gs,vs])
        call copy_stripe_from_full_matrix(yv_f,w0%d,fg,lg,nb,nv)
        call tensor_convert(w0%d,yv_fg,wrk=w1%d,iwrk=w1%n)
        call time_start_phase(PHASE_WORK, ttot = time_yvfg )
        time_yv_fg_tot = time_yv_fg_tot + time_yvfg

        call time_start_phase(PHASE_WORK, twall = time_yofg )
        call tensor_ainit(yo_fg, [lg,no], 2, local=local, atype="TDAR", tdims=[gs,os])
        call copy_stripe_from_full_matrix(yo_f,w0%d,fg,lg,nb,no)
        call tensor_convert(w0%d,yo_fg,wrk=w1%d,iwrk=w1%n)
        call time_start_phase(PHASE_WORK, ttot = time_yofg )
        time_yo_fg_tot = time_yo_fg_tot + time_yofg

        call time_start_phase(PHASE_WORK, twall = time_xofg )
        call tensor_ainit(xo_fg, [lg,no], 2, local=local, atype="TDAR", tdims=[gs,os])
        call copy_stripe_from_full_matrix(xo_f,w0%d,fg,lg,nb,no)
        call tensor_convert(w0%d,xo_fg,wrk=w1%d,iwrk=w1%n)
        call time_start_phase(PHASE_WORK, ttot = time_xofg )
        time_xo_fg_tot = time_xo_fg_tot + time_xofg

        call time_start_phase(PHASE_WORK, twall = time_xvfg )
        call tensor_ainit(xv_fg, [lg,nv], 2, local=local, atype="TDAR", tdims=[gs,vs])
        call copy_stripe_from_full_matrix(xv_f,w0%d,fg,lg,nb,nv)
        call tensor_convert(w0%d,xv_fg,wrk=w1%d,iwrk=w1%n)
        call time_start_phase(PHASE_WORK, ttot = time_xvfg )
        time_xv_fg_tot = time_xv_fg_tot + time_xvfg

        call time_start_phase(PHASE_WORK, twall = time_uigcj )
        call tensor_ainit( uigcj, [no,lg,nv,no], 4, local=local, atype='TDAR', tdims=[os,gs,vs,os] )
        order4 = [3,1,2,4]
        call tensor_contract( 1.0E0_realk, yv_fg, u2, [2],[1],1,0.0E0_realk,uigcj, order4, force_sync=.true.,wrk=w1%d,iwrk=w1%n)
        call time_start_phase(PHASE_WORK, ttot = time_uigcj )
        time_uigcj_tot = time_uigcj_tot + time_uigcj


        !**********************************
        ! Begin the loop over alpha batches
        !**********************************
        BatchAlpha: do alphaB=1,nbatchesAlpha ! AO batches

           !short hand notation
           fa = 1 + (alphaB-1)*MaxActualDimAlpha
           la = nb - fa + 1
           if( la >= MaxActualDimAlpha )then
              la = MaxActualDimAlpha
           endif

           if( la > bs )then
              as = bs
           else
              as = la
           endif

           myload     = myload + la * lg

           if(master.and.DECinfo%PL>2) write (*, '("starting job (",I3,"/",I3,",",I3,"/",I3,")")')&
              &alphaB,nbatchesAlpha,gammaB,nbatchesGamma


           call time_start_phase(PHASE_WORK, twall = time_xofa )
           call tensor_ainit(xo_fa, [la,no], 2, local=local, atype="TDAR", tdims=[as,os])
           call copy_stripe_from_full_matrix(xo_f,w0%d,fa,la,nb,no)
           call tensor_convert(w0%d,xo_fa,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_xofa )
           time_xo_fa_tot = time_xo_fa_tot + time_xofa

           call time_start_phase(PHASE_WORK, twall = time_xvfa )
           call tensor_ainit(xv_fa, [la,nv], 2, local=local, atype="TDAR", tdims=[as,vs])
           call copy_stripe_from_full_matrix(xv_f,w0%d,fa,la,nb,nv)
           call tensor_convert(w0%d,xv_fa,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_xvfa )
           time_xv_fa_tot = time_xv_fa_tot + time_xvfa

           call time_start_phase(PHASE_WORK, twall = time_yvfa )
           call tensor_ainit(yv_fa, [la,nv], 2, local=local, atype="TDAR", tdims=[as,vs])
           call copy_stripe_from_full_matrix(yv_f,w0%d,fa,la,nb,nv)
           call tensor_convert(w0%d,yv_fa,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_yvfa )
           time_yv_fa_tot = time_yv_fa_tot + time_yvfa

           call time_start_phase(PHASE_WORK, twall = time_cont1 )
           call tensor_ainit( int1, [nv,la,no,lg], 4, local=local, atype="TDAR", tdims=[vs,as,os,gs] )
           order4 = [3,4,1,2]
           call tensor_contract( 1.0E0_realk, uigcj, xo_fa, [4],[2],1,0.0E0_realk, int1, order4,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont1 )
           time_cont1_tot = time_cont1_tot + time_cont1

           call time_start_phase(PHASE_WORK, twall = time_int1 )
           call tensor_ainit(Cint, [nb,nb,la,lg], 4, local=local, atype="TDAR", tdims=[bs,bs,as,gs])
           do i = 1, Cint%nlti

              call get_midx(Cint%ti(i)%gt,starts,Cint%ntpm,Cint%mode)

              ndimA  = Cint%ti(i)%d(1)
              ndimB  = Cint%ti(i)%d(2)
              ndimC  = Cint%ti(i)%d(3)
              ndimD  = Cint%ti(i)%d(4)

              startA = 1  + (starts(1)-1)*Cint%tdim(1)
              startB = 1  + (starts(2)-1)*Cint%tdim(2)
              startC = fa + (starts(3)-1)*Cint%tdim(3)
              startD = fg + (starts(4)-1)*Cint%tdim(4)

              call II_GET_ERI_INTEGRALBLOCK_INQUIRE(DECinfo%output,DECinfo%output,Mylsitem%setting,&
                 & startA,startB,startC,startD,ndimA,ndimB,ndimC,ndimD,&
                 & ndimAs,ndimBs,ndimCs,ndimDs,INTSPEC)

              call II_GET_ERI_INTEGRALBLOCK(DECinfo%output,DECinfo%output,Mylsitem%setting,&
                 & startA,startB,startC,startD,ndimA,ndimB,ndimC,ndimD,&
                 & ndimAs,ndimBs,ndimCs,ndimDs,INTSPEC,Cint%ti(i)%t,w1%d)

           enddo
           call time_start_phase(PHASE_WORK, ttot = time_int1 )
           time_int1_tot = time_int1_tot + time_int1

           call time_start_phase(PHASE_WORK, twall = time_cont2 )
           ! I [beta delta alpha gamma] * Lambda^p [beta l] = I [alpha l gamma delta]
           call tensor_ainit( int2, [la,no,lg,nb], 4, local=local, atype="TDAR", tdims=[as,os,gs,bs] )
           order4 = [2,4,3,1]
           call tensor_contract( 1.0E0_realk, Cint, xo, [2],[1],1,0.0E0_realk, int2, order4,force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont2 )
           time_cont2_tot = time_cont2_tot + time_cont2

           !u [b alpha k gamma] * I [alpha k gamma delta] =+ Had [a delta]
           call time_start_phase(PHASE_WORK, twall = time_cont3 )
           order2 = [1,2]
           call tensor_contract( 1.0E0_realk, int1, int2, [2,3,4],[1,2,3],3,1.0E0_realk, Had, order2,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont3 )
           time_cont3_tot = time_cont3_tot + time_cont3

           call tensor_free( int1 )
           call tensor_free( Cint )


           !VVOO
           if ( Ccmodel > MODEL_CC2 ) then

              !I [alpha  i gamma delta] * Lambda^h [delta j]          = I [alpha i gamma j]
              call time_start_phase(PHASE_WORK, twall = time_cont4 )
              call tensor_ainit( int3, [la,no,no,lg], 4, local=local, atype="TDAR", tdims=[as,os,os,gs] )
              order4 = [1,2,4,3]
              call tensor_contract( 1.0E0_realk, int2, yo   , [4],[1],1, 0.0E0_realk, int3, order4,&
                 &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
              call time_start_phase(PHASE_WORK, ttot = time_cont4 )
              time_cont4_tot = time_cont4_tot + time_cont4

              !I [alpha  i j gamma] * Lambda^h [gamma b]            = I [alpha i j b]
              call time_start_phase(PHASE_WORK, twall = time_cont5 )
              call tensor_ainit( int4, [la,no,no,nv], 4, local=local, atype="TDAR", tdims=[as,os,os,vs] )
              order4 = [1,2,3,4]
              call tensor_contract( 1.0E0_realk, int3, yv_fg, [4],[1],1, 0.0E0_realk, int4, order4,&
                 &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
              call time_start_phase(PHASE_WORK, ttot = time_cont5 )
              time_cont5_tot = time_cont5_tot + time_cont5

              call tensor_free( int3 )

              !Lambda^p [alpha a]^T * I [alpha i j b]             =+ gvvoo [a i j b]
              call time_start_phase(PHASE_WORK, twall = time_cont6 )
              call tensor_contract( 1.0E0_realk, xv_fa, int4, [1],[1],1, 1.0E0_realk, gvvoo, order4,&
                 &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
              call time_start_phase(PHASE_WORK, ttot = time_cont6 )
              time_cont6_tot = time_cont6_tot + time_cont6

              call tensor_free( int4 )

           endif


           ! I [alpha l gamma delta] * Lambda^h [delta c] = I[alpha l gamma c]
           call time_start_phase(PHASE_WORK, twall = time_cont7 )
           call tensor_ainit( int1, [la,no,lg,nv], 4, local=local, atype="TDAR", tdims=[as,os,gs,vs] )
           order4 = [1,2,3,4]
           call tensor_contract( 1.0E0_realk, int2, yv, [4],[1],1, 0.0E0_realk, int1, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont7 )
           time_cont7_tot = time_cont7_tot + time_cont7

           call tensor_free( int2 )

           !I [alpha l gamma c] * u [l gamma c j]  =+ Gbi [alpha j]
           call time_start_phase(PHASE_WORK, twall = time_cont8 )
           call tensor_ainit( int3, [la,no], 2, local=local, atype="TDAR", tdims=[as,os] )
           order2 = [1,2]
           call tensor_contract( 1.0E0_realk, int1, uigcj, [2,3,4],[1,2,3],3, 0.0E0_realk, int3, order2,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont8 )
           time_cont8_tot = time_cont8_tot + time_cont8

           !FIXME: find a more elegant solution for this contraction
           call time_start_phase(PHASE_WORK, twall = time_convGBI )
           if( me == 0 )then
              call tensor_convert(int3,w0%d,wrk=w1%d,iwrk=w1%n)
              call ass_D1to2(w0%d,p2,[la,no])
              !$OMP WORKSHARE
              Gbi_local%elm2(fa:fa+la-1,:) = Gbi_local%elm2(fa:fa+la-1,:) + p2(:,:)
              !$OMP END WORKSHARE
           endif
           call time_start_phase(PHASE_WORK, ttot = time_convGBI )
           time_convGBI_tot = time_convGBI_tot + time_convGBI

           call tensor_free( int3 )

           if ( Ccmodel > MODEL_CC2 ) then

              call time_start_phase(PHASE_WORK, twall = time_cont9 )
              call tensor_ainit( int3, [la,no,nv,no], 4, local=local, atype="TDAR", tdims=[as,os,vs,os] )

              order4 = [1,2,3,4]
              !Reorder I [alpha j gamma b]  * Lambda^h [gamma i]          = I [alpha j b i]
              call tensor_contract( 1.0E0_realk, int1, yo_fg, [3],[1],1, 0.0E0_realk, int3, order4,&
                 &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
              call time_start_phase(PHASE_WORK, ttot = time_cont9 )
              time_cont9_tot = time_cont9_tot + time_cont9

              !Lambda^p [alpha a]^T * I [alpha j b i]             =+ gvoov [a j b i]
              call time_start_phase(PHASE_WORK, twall = time_cont10 )
              call tensor_contract( 1.0E0_realk, xv_fa, int3, [1],[1],1, 1.0E0_realk, gvoov, order4,&
                 &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
              call time_start_phase(PHASE_WORK, ttot = time_cont10 )
              time_cont10_tot = time_cont10_tot + time_cont10


              call tensor_free( int3 )
           endif

           call tensor_free( int1 )

           call time_start_phase(PHASE_WORK, twall = time_int2 )
           call tensor_ainit(Cint, [la,nb,lg,nb], 4, local=local, atype="TDAR", tdims=[as,bs,gs,bs])
           do i = 1, Cint%nlti

              call get_midx(Cint%ti(i)%gt,starts,Cint%ntpm,Cint%mode)

              ndimA  = Cint%ti(i)%d(1)
              ndimB  = Cint%ti(i)%d(2)
              ndimC  = Cint%ti(i)%d(3)
              ndimD  = Cint%ti(i)%d(4)

              startA = fa + (starts(1)-1)*Cint%tdim(1)
              startB = 1  + (starts(2)-1)*Cint%tdim(2)
              startC = fg + (starts(3)-1)*Cint%tdim(3)
              startD = 1  + (starts(4)-1)*Cint%tdim(4)

              call II_GET_ERI_INTEGRALBLOCK_INQUIRE(DECinfo%output,DECinfo%output,Mylsitem%setting,&
                 & startA,startB,startC,startD,ndimA,ndimB,ndimC,ndimD,&
                 & ndimAs,ndimBs,ndimCs,ndimDs,INTSPEC)

              call II_GET_ERI_INTEGRALBLOCK(DECinfo%output,DECinfo%output,Mylsitem%setting,&
                 & startA,startB,startC,startD,ndimA,ndimB,ndimC,ndimD,&
                 & ndimAs,ndimBs,ndimCs,ndimDs,INTSPEC,Cint%ti(i)%t,w1%d)

           enddo

           call time_start_phase(PHASE_WORK, ttot = time_int2 )
           time_int2_tot = time_int2_tot + time_int2

           !FIXME: REPLACE WITH SYMM AND ANTISYMM COMBINATIONS
           !Transform c
           call time_start_phase(PHASE_WORK, twall = time_cont11 )
           call tensor_ainit( int1, [lg,la,nv,nb], 4, local=local, atype="TDAR", tdims=[gs,as,vs,bs] )
           order4=[2,1,4,3]
           call tensor_contract( 1.0E0_realk, Cint, yv, [2],[1],1, 0.0E0_realk, int1, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont11 )
           time_cont11_tot = time_cont11_tot + time_cont11

           !transform d
           call time_start_phase(PHASE_WORK, twall = time_cont12 )
           call tensor_ainit( int2, [lg,la,nv,nv], 4, local=local, atype="TDAR", tdims=[gs,as,vs,vs] )
           order4=[1,2,3,4]
           call tensor_contract( 1.0E0_realk, int1, yv, [4],[1],1, 0.0E0_realk, int2, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont12 )
           time_cont12_tot = time_cont12_tot + time_cont12

           call tensor_free( int1 )
           !contract t(cd) and int2(cd)
           call time_start_phase(PHASE_WORK, twall = time_cont13 )
           call tensor_ainit( int1, [lg,la,no,no], 4, local=local, atype="TDAR", tdims=[gs,as,os,os] )
           order4=[1,2,3,4]
           call tensor_contract( 1.0E0_realk, int2, t2, [3,4],[1,2],2, 0.0E0_realk, int1, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont13 )
           time_cont13_tot = time_cont13_tot + time_cont13

           call tensor_free( int2 )

           !transform gamma to b
           call time_start_phase(PHASE_WORK, twall = time_cont14 )
           call tensor_ainit( int2, [la,nv,no,no], 4, local=local, atype="TDAR", tdims=[as,vs,os,os] )
           order4=[2,1,3,4]
           call tensor_contract( 1.0E0_realk, xv_fg, int1, [1],[1],1, 0.0E0_realk, int2, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont14 )
           time_cont14_tot = time_cont14_tot + time_cont14

           !transform alpha to a and add to omega2
           call time_start_phase(PHASE_WORK, twall = time_cont15 )
           order4=[1,2,3,4]
           call tensor_contract( 0.5E0_realk, xv_fa, int2, [1],[1],1, 1.0E0_realk, omega2, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont15 )
           time_cont15_tot = time_cont15_tot + time_cont15

           call tensor_free( int2 )

           !Do B2 term from int1 - transform occ i
           call time_start_phase(PHASE_WORK, twall = time_cont16 )
           call tensor_ainit( int2, [lg,la,no,nb], 4, local=local, atype="TDAR", tdims=[gs,as,os,bs] )
           order4=[2,1,4,3]
           call tensor_contract( 1.0E0_realk, Cint, yo, [2],[1],1, 0.0E0_realk, int2, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont16 )
           time_cont16_tot = time_cont16_tot + time_cont16

           call tensor_free( Cint )

           !transform occ j
           call time_start_phase(PHASE_WORK, twall = time_cont17 )
           order4=[1,2,3,4]
           call tensor_contract( 1.0E0_realk, int2, yo, [4],[1],1, 1.0E0_realk, int1, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont17 )
           time_cont17_tot = time_cont17_tot + time_cont17

           !backtransform t2
           order4=[1,2,3,4]
           call time_start_phase(PHASE_WORK, twall = time_cont18 )
           call tensor_ainit( int4, [nv,nv,no,la], 4, local=local, atype="TDAR", tdims=[vs,vs,os,as] )

           call tensor_contract( 1.0E0_realk, t2, xo_fa, [3],[2],1, 0.0E0_realk, int4, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont18 )
           time_cont18_tot = time_cont18_tot + time_cont18

           call time_start_phase(PHASE_WORK, twall = time_cont19 )
           call tensor_ainit( int3, [nv,nv,la,lg], 4, local=local, atype="TDAR", tdims=[vs,vs,as,gs] )

           call tensor_contract( 1.0E0_realk, int4, xo_fg, [3],[2],1, 0.0E0_realk, int3, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont19 )
           time_cont19_tot = time_cont19_tot + time_cont19

           call tensor_free( int4 )

           order4=[1,2,3,4]
           call time_start_phase(PHASE_WORK, twall = time_cont20 )
           call tensor_contract( 0.5E0_realk, int3, int1, [3,4],[2,1],2, 1.0E0_realk, omega2, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont20 )
           time_cont20_tot = time_cont20_tot + time_cont20

           call tensor_free( int3 )


           !           if( Ccmodel > MODEL_CC2 )then
           !              if(fa<=fg+lg-1)then
           !                 !CHECK WHETHER THE TERM HAS TO BE DONE AT ALL, i.e. when the first
           !                 !element in the alpha batch has a smaller index as the last element in
           !                 !the gamma batch, chose the trafolength as minimum of alpha batch-length
           !                 !and the difference between first element of alpha batch and last element
           !                 !of gamma batch
           !                 call get_a22_and_prepb22_terms_ex(w0%d,w1%d,w2%d,w3%d,tpl%d,tmi%d,no,nv,nb,fa,fg,la,lg,&
           !                    &xo,yo,xv,yv,omega2,sio4,scheme,[w0%n,w1%n,w2%n,w3%n],lock_outside,&
           !                    &time_intloop_B1work, time_intloop_B1comm, scal=0.5E0_realk  )
           !
           !                 !start a new timing phase after these terms
           !                 call time_start_phase(PHASE_WORK)
           !
           !              endif
           !           endif
           !
           !
           !           !(w0%d):I[ delta gamma alpha beta] <- (w1%d):I[ alpha beta gamma delta ]
           !           call array_reorder_4d(1.0E0_realk,w1%d,la,nb,lg,nb,[2,3,1,4],0.0E0_realk,w0%d)
           !           call lsmpi_poke()
           !           ! (w3%d):I[i gamma alpha beta] = Lambda^h[delta i] I[delta gamma alpha beta]
           !           call dgemm('t','n',no,lg*la*nb,nb,1.0E0_realk,yo,nb,w0%d,nb,0.0E0_realk,w2%d,no)
           !           call lsmpi_poke()
           !           ! (w0%d):I[i gamma alpha j] = (w3%d):I[i gamma alpha beta] Lambda^h[beta j]
           !           call dgemm('n','n',no*lg*la,no,nb,1.0E0_realk,w2%d,no*lg*la,yo,nb,0.0E0_realk,w0%d,no*lg*la)
           !           call lsmpi_poke()
           !           if( Ccmodel > MODEL_CC2 )then
           !              select case(scheme)
           !              case(4,3)
           !                 ! (w3%d):I[alpha gamma i j] <- (w0%d):I[i gamma alpha j]
           !                 call add_int_to_sio4(w0%d,w2%d,w3%d,nor,no,nv,nb,fa,fg,la,lg,xo,sio4%elm1)
           !              case(2)
           !                 ! (w3):I[ gamma i j alpha] <- (w0):I[i gamma alpha  j]
           !                 call array_reorder_4d(1.0E0_realk,w0%d,no,lg,la,no,[2,1,4,3],0.0E0_realk,w3%d)
           !                 ! (w2):I[ l i j alpha] <- (w3):Lambda^p [gamma l ]^T I[gamma i j alpha]
           !                 call dgemm('t','n',no,no*no*la,lg,1.0E0_realk,xo(fg),nb,w3%d,lg,0.0E0_realk,w2%d,no)
           !                 ! (w3):I[ k l i j] <- (w2):Lambda^p [alpha l ]^T I[ l i j , alpha]^T
           !                 call dgemm('t','t',no,no*no*no,la,1.0E0_realk,xo(fa),nb,w2%d,no*no*no,0.0E0_realk,w3%d,no)
           !                 if(.not.alloc_in_dummy)call tensor_lock_wins(sio4,'s',mode)
           !                 call tensor_add(sio4,1.0E0_realk,w3%d,order = [1,2,3,4],wrk=w2%d,iwrk=w2%n)
           !                 if( alloc_in_dummy )then
           !                    call lsmpi_win_flush(sio4%wi(1),local=.true.)
           !                 else
           !                    call tensor_unlock_wins(sio4,.true.)
           !                 endif
           !              end select
           !           endif
           !           call lsmpi_poke()
           !


           !Transform second bas to occ in int2 to get T1 transformed integral contrib
           !transform occ j
           call time_start_phase(PHASE_WORK, twall = time_cont21 )
           order4=[1,2,3,4]
           call tensor_contract( 1.0E0_realk, int2, yo, [4],[1],1, 0.0E0_realk, int1, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont21 )
           time_cont21_tot = time_cont21_tot + time_cont21

           !g´transform gamma->b
           call time_start_phase(PHASE_WORK, twall = time_cont22 )
           call tensor_ainit( int4, [la,nv,no,no], 4, local=local, atype="TDAR", tdims=[as,vs,os,os] )
           order4=[2,1,3,4]
           call tensor_contract( 1.0E0_realk, xv_fg, int1, [1],[1],1, 0.0E0_realk, int4, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont22 )
           time_cont22_tot = time_cont22_tot + time_cont22

           call time_start_phase(PHASE_WORK, twall = time_cont23 )
           order4=[1,2,3,4]
           call tensor_contract( 0.5E0_realk, xv_fa, int4, [1],[1],1, 1.0E0_realk, omega2, order4,&
              &force_sync=.true.,wrk=w1%d,iwrk=w1%n)
           call time_start_phase(PHASE_WORK, ttot = time_cont23 )
           time_cont23_tot = time_cont23_tot + time_cont23


           call tensor_free( int1 )
           call tensor_free( int2 )
           call tensor_free( int4 )


           call tensor_free( xo_fa )
           call tensor_free( xv_fa )
           call tensor_free( yv_fa )

        end do BatchAlpha

        call tensor_free( uigcj )
        call tensor_free( yv_fg )
        call tensor_free( yo_fg )
        call tensor_free( xo_fg )
        call tensor_free( xv_fg )

     end do BatchGamma

     call time_start_phase( PHASE_WORK, ttot = tot_intloop )

     mylsitem%setting%SCHEME%CS_SCREEN = save_cs_screen
     mylsitem%setting%SCHEME%PS_SCREEN = save_ps_screen

     ! free arrays only needed in the batched loops
     call time_start_phase(PHASE_COMM, at = time_intloop_work )

     if( CCmodel > MODEL_CC2 )then
        call tensor_unlock_wins(gvvoo,  all_nodes = alloc_in_dummy, check =.not.alloc_in_dummy)
        call tensor_unlock_wins(gvoov,  all_nodes = alloc_in_dummy, check =.not.alloc_in_dummy)
        call tensor_unlock_wins(sio4,   all_nodes = alloc_in_dummy, check =.not.alloc_in_dummy)
     endif

     call tensor_sync_replicated(Gbi_local)
     call tensor_convert(Gbi_local%elm1,Gbi)
     call tensor_free(Gbi_local)

     call time_start_phase( PHASE_WORK )
     call time_phases_get_diff(current_wt=phase_cntrs)

     call tensor_free(tpl)
     call tensor_free(tmi)
     ! free working matrices and adapt to new requirements
     call mem_dealloc(w0)
     call mem_dealloc(w1)

     if(DECinfo%PL>2)then

        tot_intloop_min      =  tot_intloop 
        time_yv_fg_tot_min   =  time_yv_fg_tot   
        time_yo_fg_tot_min   =  time_yo_fg_tot   
        time_xo_fg_tot_min   =  time_xo_fg_tot   
        time_xv_fg_tot_min   =  time_xv_fg_tot   
        time_uigcj_tot_min   =  time_uigcj_tot   
        time_xo_fa_tot_min   =  time_xo_fa_tot   
        time_xv_fa_tot_min   =  time_xv_fa_tot   
        time_yv_fa_tot_min   =  time_yv_fa_tot   
        time_cont1_tot_min   =  time_cont1_tot   
        time_int1_tot_min    =  time_int1_tot   
        time_cont2_tot_min   =  time_cont2_tot   
        time_cont3_tot_min   =  time_cont3_tot   
        time_cont4_tot_min   =  time_cont4_tot   
        time_cont5_tot_min   =  time_cont5_tot   
        time_cont6_tot_min   =  time_cont6_tot   
        time_cont7_tot_min   =  time_cont7_tot   
        time_cont8_tot_min   =  time_cont8_tot   
        time_convGBI_tot_min =  time_convGBI_tot 
        time_cont9_tot_min   =  time_cont9_tot 
        time_cont10_tot_min  =  time_cont10_tot  
        time_int2_tot_min    =  time_int2_tot    
        time_cont11_tot_min  =  time_cont11_tot
        time_cont12_tot_min  =  time_cont12_tot
        time_cont13_tot_min  =  time_cont13_tot
        time_cont14_tot_min  =  time_cont14_tot
        time_cont15_tot_min  =  time_cont15_tot
        time_cont16_tot_min  =  time_cont16_tot
        time_cont17_tot_min  =  time_cont17_tot
        time_cont18_tot_min  =  time_cont18_tot
        time_cont19_tot_min  =  time_cont19_tot
        time_cont20_tot_min  =  time_cont20_tot
        time_cont21_tot_min  =  time_cont21_tot
        time_cont22_tot_min  =  time_cont22_tot
        time_cont23_tot_min  =  time_cont23_tot

        tot_intloop_max      =  tot_intloop 
        time_yv_fg_tot_max   =  time_yv_fg_tot   
        time_yo_fg_tot_max   =  time_yo_fg_tot   
        time_xo_fg_tot_max   =  time_xo_fg_tot   
        time_xv_fg_tot_max   =  time_xv_fg_tot   
        time_uigcj_tot_max   =  time_uigcj_tot   
        time_xo_fa_tot_max   =  time_xo_fa_tot   
        time_xv_fa_tot_max   =  time_xv_fa_tot   
        time_yv_fa_tot_max   =  time_yv_fa_tot   
        time_cont1_tot_max   =  time_cont1_tot   
        time_int1_tot_max    =  time_int1_tot   
        time_cont2_tot_max   =  time_cont2_tot   
        time_cont3_tot_max   =  time_cont3_tot   
        time_cont4_tot_max   =  time_cont4_tot   
        time_cont5_tot_max   =  time_cont5_tot   
        time_cont6_tot_max   =  time_cont6_tot   
        time_cont7_tot_max   =  time_cont7_tot   
        time_cont8_tot_max   =  time_cont8_tot   
        time_convGBI_tot_max =  time_convGBI_tot 
        time_cont9_tot_max   =  time_cont9_tot 
        time_cont10_tot_max  =  time_cont10_tot  
        time_int2_tot_max    =  time_int2_tot    
        time_cont11_tot_max  =  time_cont11_tot
        time_cont12_tot_max  =  time_cont12_tot
        time_cont13_tot_max  =  time_cont13_tot
        time_cont14_tot_max  =  time_cont14_tot
        time_cont15_tot_max  =  time_cont15_tot
        time_cont16_tot_max  =  time_cont16_tot
        time_cont17_tot_max  =  time_cont17_tot
        time_cont18_tot_max  =  time_cont18_tot
        time_cont19_tot_max  =  time_cont19_tot
        time_cont20_tot_max  =  time_cont20_tot
        time_cont21_tot_max  =  time_cont21_tot
        time_cont22_tot_max  =  time_cont22_tot
        time_cont23_tot_max  =  time_cont23_tot


        call lsmpi_reduce_realk_min( tot_intloop_min      , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_yv_fg_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_yo_fg_tot_min   , infpar%master, infpar%lg_comm  )
        call lsmpi_reduce_realk_min( time_xo_fg_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_xv_fg_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_uigcj_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_xo_fa_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_xv_fa_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_yv_fa_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont1_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_int1_tot_min    , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont2_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont3_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont4_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont5_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont6_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont7_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont8_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_convGBI_tot_min , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont9_tot_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont10_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_int2_tot_min    , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont11_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont12_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont13_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont14_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont15_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont16_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont17_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont18_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont19_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont20_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont21_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont22_tot_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_cont23_tot_min  , infpar%master, infpar%lg_comm )

        call lsmpi_reduce_realk_max( tot_intloop_max      , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_yv_fg_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_yo_fg_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_xo_fg_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_xv_fg_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_uigcj_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_xo_fa_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_xv_fa_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_yv_fa_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont1_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_int1_tot_max    , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont2_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont3_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont4_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont5_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont6_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont7_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont8_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_convGBI_tot_max , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont9_tot_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont10_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_int2_tot_max    , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont11_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont12_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont13_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont14_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont15_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont16_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont17_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont18_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont19_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont20_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont21_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont22_tot_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_cont23_tot_max  , infpar%master, infpar%lg_comm )

        call lsmpi_local_reduction( tot_intloop      , infpar%master )
        call lsmpi_local_reduction( time_yv_fg_tot   , infpar%master )
        call lsmpi_local_reduction( time_yo_fg_tot   , infpar%master )
        call lsmpi_local_reduction( time_xo_fg_tot   , infpar%master )
        call lsmpi_local_reduction( time_xv_fg_tot   , infpar%master )
        call lsmpi_local_reduction( time_uigcj_tot   , infpar%master )
        call lsmpi_local_reduction( time_xo_fa_tot   , infpar%master )
        call lsmpi_local_reduction( time_xv_fa_tot   , infpar%master )
        call lsmpi_local_reduction( time_yv_fa_tot   , infpar%master )
        call lsmpi_local_reduction( time_cont1_tot   , infpar%master )
        call lsmpi_local_reduction( time_int1_tot    , infpar%master )
        call lsmpi_local_reduction( time_cont2_tot   , infpar%master )
        call lsmpi_local_reduction( time_cont3_tot   , infpar%master )
        call lsmpi_local_reduction( time_cont4_tot   , infpar%master )
        call lsmpi_local_reduction( time_cont5_tot   , infpar%master )
        call lsmpi_local_reduction( time_cont6_tot   , infpar%master )
        call lsmpi_local_reduction( time_cont7_tot   , infpar%master )
        call lsmpi_local_reduction( time_cont8_tot   , infpar%master )
        call lsmpi_local_reduction( time_convGBI_tot , infpar%master )
        call lsmpi_local_reduction( time_cont9_tot   , infpar%master )
        call lsmpi_local_reduction( time_cont10_tot  , infpar%master )
        call lsmpi_local_reduction( time_int2_tot    , infpar%master )
        call lsmpi_local_reduction( time_cont11_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont12_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont13_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont14_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont15_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont16_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont17_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont18_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont19_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont20_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont21_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont22_tot  , infpar%master )
        call lsmpi_local_reduction( time_cont23_tot  , infpar%master )

        unlock_time   = time_lsmpi_win_unlock - unlock_time
        waiting_time  = time_lsmpi_wait       - waiting_time
        flushing_time = time_lsmpi_win_flush  - flushing_time

        unlock_time_min    = unlock_time
        waiting_time_min   = waiting_time      
        flushing_time_min  = flushing_time 

        unlock_time_max    = unlock_time
        waiting_time_max   = waiting_time      
        flushing_time_max  = flushing_time 

        call lsmpi_reduce_realk_min( unlock_time_min   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( waiting_time_min  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( flushing_time_min , infpar%master, infpar%lg_comm )

        call lsmpi_reduce_realk_max( unlock_time_max   , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( waiting_time_max  , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( flushing_time_max , infpar%master, infpar%lg_comm )

        call lsmpi_local_reduction( unlock_time   , infpar%master )
        call lsmpi_local_reduction( waiting_time  , infpar%master )
        call lsmpi_local_reduction( flushing_time , infpar%master )

        time_w_min = phase_cntrs( PHASE_WORK_IDX )
        time_c_min = phase_cntrs( PHASE_COMM_IDX )
        time_i_min = phase_cntrs( PHASE_IDLE_IDX )

        time_w_max = phase_cntrs( PHASE_WORK_IDX )
        time_c_max = phase_cntrs( PHASE_COMM_IDX )
        time_i_max = phase_cntrs( PHASE_IDLE_IDX )

        call lsmpi_reduce_realk_min( time_w_min , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_c_min , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_min( time_i_min , infpar%master, infpar%lg_comm )

        call lsmpi_reduce_realk_max( time_w_max , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_c_max , infpar%master, infpar%lg_comm )
        call lsmpi_reduce_realk_max( time_i_max , infpar%master, infpar%lg_comm )

        call lsmpi_local_reduction(phase_cntrs,nphases,infpar%master)

        if(master)then
           write(*,'("CCSD total time               phase B",g10.3,g10.3,g10.3)')&
              & tot_intloop_max     ,      tot_intloop      /dble(nnod),tot_intloop_min 
           write(*,'("CCSD time_yv_fg_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_yv_fg_tot_max  ,  time_yv_fg_tot   /dble(nnod),time_yv_fg_tot_min  ,  time_yv_fg_tot   / tot_intloop
           write(*,'("CCSD time_yo_fg_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_yo_fg_tot_max  ,  time_yo_fg_tot   /dble(nnod),time_yo_fg_tot_min  ,  time_yo_fg_tot   / tot_intloop
           write(*,'("CCSD time_xo_fg_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_xo_fg_tot_max  ,  time_xo_fg_tot   /dble(nnod),time_xo_fg_tot_min  ,  time_xo_fg_tot   / tot_intloop
           write(*,'("CCSD time_xv_fg_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_xv_fg_tot_max  ,  time_xv_fg_tot   /dble(nnod),time_xv_fg_tot_min  ,  time_xv_fg_tot   / tot_intloop
           write(*,'("CCSD time_uigcj_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_uigcj_tot_max  ,  time_uigcj_tot   /dble(nnod),time_uigcj_tot_min  ,  time_uigcj_tot   / tot_intloop
           write(*,'("CCSD time_xo_fa_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_xo_fa_tot_max  ,  time_xo_fa_tot   /dble(nnod),time_xo_fa_tot_min  ,  time_xo_fa_tot   / tot_intloop
           write(*,'("CCSD time_xv_fa_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_xv_fa_tot_max  ,  time_xv_fa_tot   /dble(nnod),time_xv_fa_tot_min  ,  time_xv_fa_tot   / tot_intloop
           write(*,'("CCSD time_yv_fa_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_yv_fa_tot_max  ,  time_yv_fa_tot   /dble(nnod),time_yv_fa_tot_min  ,  time_yv_fa_tot   / tot_intloop
           write(*,'("CCSD time_cont1_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont1_tot_max  ,  time_cont1_tot   /dble(nnod),time_cont1_tot_min  ,  time_cont1_tot   / tot_intloop
           write(*,'("CCSD time_int1_tot            phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_int1_tot_max   ,  time_int1_tot    /dble(nnod),time_int1_tot_min   ,  time_int1_tot    / tot_intloop
           write(*,'("CCSD time_cont2_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont2_tot_max  ,  time_cont2_tot   /dble(nnod),time_cont2_tot_min  ,  time_cont2_tot   / tot_intloop
           write(*,'("CCSD time_cont3_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont3_tot_max  ,  time_cont3_tot   /dble(nnod),time_cont3_tot_min  ,  time_cont3_tot   / tot_intloop
           write(*,'("CCSD time_cont4_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont4_tot_max  ,  time_cont4_tot   /dble(nnod),time_cont4_tot_min  ,  time_cont4_tot   / tot_intloop
           write(*,'("CCSD time_cont5_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont5_tot_max  ,  time_cont5_tot   /dble(nnod),time_cont5_tot_min  ,  time_cont5_tot   / tot_intloop
           write(*,'("CCSD time_cont6_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont6_tot_max  ,  time_cont6_tot   /dble(nnod),time_cont6_tot_min  ,  time_cont6_tot   / tot_intloop
           write(*,'("CCSD time_cont7_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont7_tot_max  ,  time_cont7_tot   /dble(nnod),time_cont7_tot_min  ,  time_cont7_tot   / tot_intloop
           write(*,'("CCSD time_cont8_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont8_tot_max  ,  time_cont8_tot   /dble(nnod),time_cont8_tot_min  ,  time_cont8_tot   / tot_intloop
           write(*,'("CCSD time_convGBI_tot         phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_convGBI_tot_max,  time_convGBI_tot /dble(nnod),time_convGBI_tot_min,  time_convGBI_tot / tot_intloop
           write(*,'("CCSD time_cont9_tot           phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont9_tot_max  ,  time_cont9_tot   /dble(nnod),time_cont9_tot_min  ,  time_cont9_tot   / tot_intloop
           write(*,'("CCSD time_cont10_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont10_tot_max ,  time_cont10_tot  /dble(nnod),time_cont10_tot_min ,  time_cont10_tot  / tot_intloop
           write(*,'("CCSD time_int2_tot            phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_int2_tot_max   ,  time_int2_tot    /dble(nnod),time_int2_tot_min   ,  time_int2_tot    / tot_intloop
           write(*,'("CCSD time_cont11_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont11_tot_max ,  time_cont11_tot  /dble(nnod),time_cont11_tot_min ,  time_cont11_tot  / tot_intloop
           write(*,'("CCSD time_cont12_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont12_tot_max ,  time_cont12_tot  /dble(nnod),time_cont12_tot_min ,  time_cont12_tot  / tot_intloop
           write(*,'("CCSD time_cont13_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont13_tot_max ,  time_cont13_tot  /dble(nnod),time_cont13_tot_min ,  time_cont13_tot  / tot_intloop
           write(*,'("CCSD time_cont14_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont14_tot_max ,  time_cont14_tot  /dble(nnod),time_cont14_tot_min ,  time_cont14_tot  / tot_intloop
           write(*,'("CCSD time_cont15_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont15_tot_max ,  time_cont15_tot  /dble(nnod),time_cont15_tot_min ,  time_cont15_tot  / tot_intloop
           write(*,'("CCSD time_cont16_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont16_tot_max ,  time_cont16_tot  /dble(nnod),time_cont16_tot_min ,  time_cont16_tot  / tot_intloop
           write(*,'("CCSD time_cont17_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont17_tot_max ,  time_cont17_tot  /dble(nnod),time_cont17_tot_min ,  time_cont17_tot  / tot_intloop
           write(*,'("CCSD time_cont18_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont18_tot_max ,  time_cont18_tot  /dble(nnod),time_cont18_tot_min ,  time_cont18_tot  / tot_intloop
           write(*,'("CCSD time_cont19_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont19_tot_max ,  time_cont19_tot  /dble(nnod),time_cont19_tot_min ,  time_cont19_tot  / tot_intloop
           write(*,'("CCSD time_cont20_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont20_tot_max ,  time_cont20_tot  /dble(nnod),time_cont20_tot_min ,  time_cont20_tot  / tot_intloop
           write(*,'("CCSD time_cont21_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont21_tot_max ,  time_cont21_tot  /dble(nnod),time_cont21_tot_min ,  time_cont21_tot  / tot_intloop
           write(*,'("CCSD time_cont22_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont22_tot_max ,  time_cont22_tot  /dble(nnod),time_cont22_tot_min ,  time_cont22_tot  / tot_intloop
           write(*,'("CCSD time_cont23_tot          phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_cont23_tot_max ,  time_cont23_tot  /dble(nnod),time_cont23_tot_min ,  time_cont23_tot  / tot_intloop
           write(*,'("CCSD time in lsmpi_win_unlock phase B",g10.3,g10.3,g10.3,g10.3)')&
              & unlock_time_max, unlock_time/dble(nnod),unlock_time_min,          unlock_time   / tot_intloop
           write(*,'("CCSD time in lsmpi_wait       phase B",g10.3,g10.3,g10.3,g10.3)')&
              & waiting_time_max,   waiting_time  /dble(nnod),waiting_time_min,   waiting_time  / tot_intloop
           write(*,'("CCSD time in lsmpi_win_flush  phase B",g10.3,g10.3,g10.3,g10.3)')&
              & flushing_time_max,  flushing_time /dble(nnod),flushing_time_min,  flushing_time / tot_intloop
           write(*,'("CCSD time WORK                phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_w_max,phase_cntrs(PHASE_WORK_IDX)/dble(nnod),time_w_min,phase_cntrs(PHASE_WORK_IDX)/tot_intloop
           write(*,'("CCSD time COMM                phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_c_max,phase_cntrs(PHASE_COMM_IDX)/dble(nnod),time_c_min,phase_cntrs(PHASE_COMM_IDX)/tot_intloop
           write(*,'("CCSD time IDLE                phase B",g10.3,g10.3,g10.3,g10.3)')&
              & time_i_max,phase_cntrs(PHASE_IDLE_IDX)/dble(nnod),time_i_min,phase_cntrs(PHASE_IDLE_IDX)/tot_intloop
        endif
     endif
     unlock_time   = time_lsmpi_win_unlock 
     waiting_time  = time_lsmpi_wait
     flushing_time = time_lsmpi_win_flush

     call ccsd_debug_print(ccmodel,1,master,local,scheme,print_debug,o2v2,w1,&
        &omega2,govov,gvvoo,gvoov)

     if(Ccmodel>MODEL_CC2)then


        !get B2.2 contributions - CURRENTLY NOT NECESSARY
        !**********************
        !call get_B22_contrib_mo(sio4,t2,w1%d,w2%d,no,nv,omega2,scheme,lock_outside,&
        !   &time_Bcnd_work,time_Bcnd_comm)

        call tensor_free(sio4)

        !Get the C2 and D2 terms
        !***********************
        call get_cnd_terms_mo_2(w1%d,w2%d,w3%d,t2,u2,govov,gvoov,gvvoo,no,nv,omega2,&
           &scheme,lock_outside,local)

        call ccsd_debug_print(ccmodel,3,master,local,scheme,print_debug,o2v2,w1,&
           &omega2,govov,gvvoo,gvoov)

        call tensor_free( gvoov )
        call tensor_free( gvvoo )

     endif

     call tensor_unlock_wins(omega2, all_nodes = alloc_in_dummy, check =.not.alloc_in_dummy)

     if(alloc_in_dummy)then
        call tensor_unlock_wins(u2,     all_nodes=.true.)
     endif

     govov%access_type  = AT_MASTER_ACCESS
     omega2%access_type = AT_MASTER_ACCESS
     t2%access_type     = AT_MASTER_ACCESS
     u2%access_type     = AT_MASTER_ACCESS
     Gbi%access_type    = AT_MASTER_ACCESS
     Had%access_type    = AT_MASTER_ACCESS
     xo%access_type     = AT_MASTER_ACCESS
     xv%access_type     = AT_MASTER_ACCESS
     yo%access_type     = AT_MASTER_ACCESS
     yv%access_type     = AT_MASTER_ACCESS


     if(.not. master) then

        return

     endif

#ifdef VAR_LSDEBUG
     if(print_debug)then
        call print_norm(Gbi,    " NORM(Gbi)       :",print_on_rank=0)
        call print_norm(Had,    " NORM(Had)       :",print_on_rank=0)
        call print_norm(omega2, " NORM(omega2 s-o):",print_on_rank=0)
        call print_norm(govov,  " NORM(govov s-o) :",print_on_rank=0)
     endif
#endif


     call tensor_minit(iFock,[nb,nb],2,atype="REAR",local=local)
     call tensor_zero(iFock)

     call mem_alloc(w0, max(max(i8*nb*nb,i8*nb*no),i8*nb*nv), simple=.true.)

     !calculate inactive fock matrix in ao basis
     call dgemm('n','t',nb,nb,no,1.0E0_realk,yo_f,nb,xo_f,nb,0.0E0_realk,w0%d,nb)

     call II_get_fock_mat_full(DECinfo%output,DECinfo%output,MyLsItem%setting,nb,w0%d,.false.,iFock%elm1)
     if(DECinfo%DFTreference)then
        call lsquit("ERROR(yet_another_ccsd_residual): DFT ref not implemented",-1)
        !call II_get_xc_fock_mat_full(DECinfo%output,DECinfo%output,MyLsItem%setting,nb,Dens%elms,.false.,iFock%elms)
     endif

     !use dens as temporay array 
     call ii_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsItem%setting,w0%d,nb,nb,AORdefault,AORdefault)
     ! Add one- and two-electron contributions to Fock matrix
     call daxpy(nb2,1.0E0_realk,w0%d,1,iFock%elm1,1)
     !Free the density matrix

     ! KK: Add long-range Fock correction
     call daxpy(nb2,1.0E0_realk,deltafock,1,iFock%elm1,1)

     call tensor_sync_replicated(iFock)

     if(print_debug)then
        call print_norm(deltafock,int((i8*nb)*nb,kind=8), " NORM(deltafock):")
        call print_norm(iFock, " NORM(iFock)    :")
     endif



     !Transform inactive Fock matrix into the different mo subspaces
     if (Ccmodel>MODEL_CC2) then
        ! -> Foo
        call dgemm('t','n',no,nb,nb,1.0E0_realk,xo_f,nb,iFock%elm1,nb,0.0E0_realk,w0%d,no)
        call dgemm('n','n',no,no,nb,1.0E0_realk,w0%d,no,yo_f,nb,0.0E0_realk,ppfock,no)
        ! -> Fov
        call dgemm('n','n',no,nv,nb,1.0E0_realk,w0%d,no,yv_f,nb,0.0E0_realk,pqfock,no)
        ! -> Fvo
        call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv_f,nb,iFock%elm1,nb,0.0E0_realk,w0%d,nv)
        call dgemm('n','n',nv,no,nb,1.0E0_realk,w0%d,nv,yo_f,nb,0.0E0_realk,qpfock,nv)
        ! -> Fvv
        call dgemm('n','n',nv,nv,nb,1.0E0_realk,w0%d,nv,yv_f,nb,0.0E0_realk,qqfock,nv)
     else
        ! -> Foo
        call dgemm('t','n',no,nb,nb,1.0E0_realk,xo_f,nb,fock,nb,0.0E0_realk,w0%d,no)
        call dgemm('n','n',no,no,nb,1.0E0_realk,w0%d,no,yo_f,nb,0.0E0_realk,ppfock,no)
        ! -> Fov
        call dgemm('t','n',no,nb,nb,1.0E0_realk,xo_f,nb,iFock%elm1,nb,0.0E0_realk,w0%d,no)
        call dgemm('n','n',no,nv,nb,1.0E0_realk,w0%d,no,yv_f,nb,0.0E0_realk,pqfock,no)
        ! -> Fvo
        call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv_f,nb,iFock%elm1,nb,0.0E0_realk,w0%d,nv)
        call dgemm('n','n',nv,no,nb,1.0E0_realk,w0%d,nv,yo_f,nb,0.0E0_realk,qpfock,nv)
        ! -> Fvv
        call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv_f,nb,fock,nb,0.0E0_realk,w0%d,nv)
        call dgemm('n','n',nv,nv,nb,1.0E0_realk,w0%d,nv,yv_f,nb,0.0E0_realk,qqfock,nv)
     endif



     if(print_debug)then
        call print_norm(ppfock,int((i8*no)*no,kind=8)," NORM(ppfock)   :")
        call print_norm(pqfock,int((i8*no)*nv,kind=8)," NORM(pqfock)   :")
        call print_norm(qpfock,int((i8*no)*nv,kind=8)," NORM(qpfock)   :")
        call print_norm(qqfock,int((i8*nv)*nv,kind=8)," NORM(qqfock)   :")
     endif

     !Free the AO fock matrix
     call tensor_free(iFock)



     !CCD can be achieved by not using singles residual updates here
     if(.not. DECinfo%CCDhack)then

        !GET SINGLES CONTRIBUTIONS
        !*************************

        !FIXME:this is a workaround, use complete pdm also in the trafos
        call tensor_minit(int1, [nv,no], 2, atype="TDAR", tdims=[vs,os],local=local)
        call tensor_minit(int2, [no,nv], 2, atype="TDAR", tdims=[os,vs],local=local)
        call tensor_convert(pqfock,int2,wrk=w0%d,iwrk=w0%n)

        ! u [a i k c] * F[k c] =+ Omega [a i]
        order2=[1,2]
        call tensor_contract(1.0E0_realk,u2,int2,[1,4],[2,1],2,0.0E0_realk,int1,order2,force_sync=.true.)

        call tensor_free(int2)

        !calculate singles G term
        ! Lambda^p [alpha a]^T Gbi [alpha i] =+ Omega [a i]
        call tensor_contract(1.0E0_realk,xv,Gbi,[1],[1],1,1.0E0_realk,int1,order2,force_sync=.true.)

        !calculate singles H term
        ! (-1) Had [a delta] * Lambda^h [delta i] =+ Omega[a i]
        call tensor_contract(-1.0E0_realk,Had,yo,[2],[1],1,1.0E0_realk,int1,order2,force_sync=.true.)

        call tensor_convert(int1,omega1,wrk = w0%d, iwrk = w0%n)

        call tensor_free(int1)

        !calculate singles J term
        ! F [a i] = Omega [a i]
        call daxpy(no*nv,1.0E0_realk,qpfock,1,omega1,1)

     endif


     !GET DOUBLES E2 TERM - AND INTRODUCE PERMUTATIONAL SYMMMETRY
     !***********************************************************
     !Prepare the E2 term by transforming Had and Cbi and move them in
     !PDM, here rendundand work is performed, but only O^3, so cheap
     call tensor_minit(int1,[nv,nv],2,tdims=[vs,vs],atype="TDAR",local=local)
     call tensor_minit(int2,[no,no],2,tdims=[os,os],atype="TDAR",local=local)
     call tensor_convert(qqfock,int1,wrk=w0%d,iwrk=w0%n)
     call tensor_convert(ppfock,int2,wrk=w0%d,iwrk=w0%n)

     order2=[1,2]
     if (Ccmodel>MODEL_CC2)then
        call tensor_contract(-1.0E0_realk,Had,yv,[2],[1],1,1.0E0_realk,int1,order2,force_sync=.true.)
        call tensor_contract(1.0E0_realk,xo,Gbi,[1],[1],1,1.0E0_realk,int2,order2,force_sync=.true.)
     endif

     order4 = [1,4,2,3]
     call tensor_contract( 1.0E0_realk,t2,int1,[2],[2],1,1.0E0_realk,omega2,order4,force_sync=.true.)
     order4 = [1,2,3,4]
     call tensor_contract(-1.0E0_realk,t2,int2,[4],[1],1,1.0E0_realk,omega2,order4,force_sync=.true.)

     call tensor_minit(int3,omega2%dims,4,tdims=omega2%tdim,atype="TDAR",local=local)

     call tensor_free(int1)
     call tensor_free(int2)

     !INTRODUCE PERMUTATION
     order4 = [2,1,4,3]
     call tensor_cp_data(omega2,int3, order = order4 )
     call tensor_add(omega2,1.0E0_realk,int3)

     call tensor_free(int3)

     !call calculate_E2_and_permute(ccmodel,ppfock,qqfock,w1%d,t2,xo,yv,Gbi,Had,no,nv,nb,&
     !   &omega2,o2v2,scheme,print_debug,lock_outside,time_Esing_work,time_Esing_comm)

     call tensor_free(Had)
     call tensor_free(Gbi)
     call tensor_free(xo)
     call tensor_free(yo)
     call tensor_free(xv)
     call tensor_free(yv)
     call tensor_free(u2)

     call mem_dealloc(w0)

     if(print_debug)then
        call print_norm(omega1,int((i8*no)*nv,kind=8)," NORM(omega1):")
        call print_norm(omega2,                       " NORM(omega2):")
     endif

#else
     call lsquit("ERROR(yet_another_ccsd_residual): MPI only",-1)
#endif

  end subroutine yet_another_ccsd_residual

  subroutine calculate_E2_and_permute(ccmodel,ppf,qqf,w1,t2,xo,yv,Gbi,Had,no,nv,nb,&
        &omega2,o2v2,s,pd,lock_outside,tw,tc)
     implicit none
     !> CC model
     integer,intent(in) :: ccmodel
     integer(kind=8),intent(in)::o2v2
     real(realk),intent(inout)::ppf(:)
     real(realk),intent(inout)::qqf(:)
     real(realk) :: w1(o2v2)
     type(tensor),intent(inout) :: t2
     real(realk),pointer:: xo(:)
     real(realk),pointer:: yv(:)
     real(realk),pointer:: Gbi(:)
     real(realk),pointer:: Had(:)
     integer, intent(in) :: no,nv,nb
     type(tensor),intent(inout) :: omega2
     integer, intent(in) :: s
     logical, intent(in) :: pd
     logical,intent(inout) :: lock_outside
     !> timing
     real(realk),intent(inout) :: tw,tc
     integer :: no2,nv2,v2o,o2v
     logical :: master 
     real(realk),pointer :: w2(:),w3(:)
     integer :: i,ml
     integer(kind=ls_mpik) :: me,nnod,nod
     integer :: ml1,fai1,l1,tl1,lai1
     integer :: ml2,fai2,l2,tl2,lai2
     integer :: fri,tri,ccmodel_copy
     character(tensor_MSG_LEN) :: msg
     real(realk) :: nrm
     integer(kind=8) :: w3size
     integer(kind=ls_mpik) :: mode
     logical :: lock_safe,traf1,traf2,trafi
     type(tensor) :: E1,E2, Pijab_om2
     integer :: os, vs, ord(4)

     call time_start_phase(PHASE_WORK)

     master       = .true.
     nrm          = 0.0E0_realk
     no2          = no*no
     nv2          = nv*nv
     v2o          = nv*nv*no
     o2v          = no*no*nv
     me           = 0
     nnod         = 1

#ifdef VAR_MPI
     master=(infpar%lg_mynum==infpar%master)
     if((s==2).and.master)then
        ccmodel_copy = ccmodel
        call time_start_phase(PHASE_COMM, at = tw)
        call share_E2_with_slaves(ccmodel_copy,ppf,qqf,t2,xo,yv,Gbi,Had,no,nv,nb,omega2,s,lock_outside)
        call time_start_phase(PHASE_WORK, at = tc)
     endif
#endif


     if(s==4.or.s==3.or.s==0)then
        !calculate first part of doubles E term and its permutation
        ! F [k j] + Lambda^p [alpha k]^T * Gbi [alpha j] = G' [k j]
        call dcopy(no2,ppf,1,w1,1)
        if (Ccmodel>MODEL_CC2) call dgemm('t','n',no,no,nb,1.0E0_realk,xo,nb,Gbi,nb,1.0E0_realk,w1,no)
        ! (-1) t [a b i k] * G' [k j] =+ Omega [a b i j]
        call dgemm('n','n',v2o,no,no,-1.0E0_realk,t2%elm1,v2o,w1,no,1.0E0_realk,omega2%elm1,v2o)

        !calculate second part of doubles E term
        ! F [b c] - Had [a delta] * Lambda^h [delta c] = H' [b c]
        call dcopy(nv2,qqf,1,w1,1)
        if (Ccmodel>MODEL_CC2) call dgemm('n','n',nv,nv,nb,-1.0E0_realk,Had,nv,yv,nb,1.0E0_realk,w1,nv)
        ! H'[a c] * t [c b i j] =+ Omega [a b i j]
        call dgemm('n','n',nv,o2v,nv,1.0E0_realk,w1,nv,t2%elm1,nv,1.0E0_realk,omega2%elm1,nv)

        if(pd) then 
           call print_norm(omega2," NORM(omega2 before permut):")
        endif

        !INTRODUCE PERMUTATION
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
        call assign_in_subblocks(w1,'=',omega2%elm1,o2v2)
#else
        !$OMP WORKSHARE
        w1(1_long:o2v2) = omega2%elm1(1_long:o2v2)
        !$OMP END WORKSHARE
#endif

        call array_reorder_4d(1.0E0_realk,w1,nv,nv,no,no,[2,1,4,3],1.0E0_realk,omega2%elm1)


#ifdef VAR_MPI
        !THE INTENSIVE SCHEMES
     else if(s==2)then
        omega2%access_type = AT_ALL_ACCESS
        t2%access_type     = AT_ALL_ACCESS
        nnod             = infpar%lg_nodtot
        me               = infpar%lg_mynum
        mode             = int(MPI_MODE_NOCHECK,kind=ls_mpik)
        vs               = t2%tdim(1)
        os               = t2%tdim(3)

        call tensor_lock_local_wins(omega2,'e',mode)

        !Prepare the E2 term by transforming Had and Cbi and move them in
        !PDM, here rendundand work is performed, but only O^3, so cheap
        call tensor_ainit(E1,[nv,nv],2,tdims=[vs,vs],atype="TDPD")
        E1%itype = TT_TILED_DIST
        call dcopy(nv2,qqf,1,E1%elm1,1)
        call tensor_lock_local_wins(E1,'e',mode)
        if (Ccmodel>MODEL_CC2) call dgemm('n','n',nv,nv,nb,-1.0E0_realk,Had,nv,yv,nb,1.0E0_realk,E1%elm1,nv)
        call tensor_mv_dense2tiled(E1,.true.)

        call tensor_ainit(E2,[no,no],2,tdims=[os,os],atype="TDPD")
        call tensor_lock_local_wins(E2,'e',mode)
        E2%itype = TT_TILED_DIST
        call dcopy(no2,ppf,1,E2%elm1,1)
        if (Ccmodel>MODEL_CC2) call dgemm('t','n',no,no,nb,1.0E0_realk,xo,nb,Gbi,nb,1.0E0_realk,E2%elm1,no)
        call tensor_mv_dense2tiled(E2,.true.)
        call tensor_unlock_local_wins(E1)
        call tensor_unlock_local_wins(E2)


        ord = [1,4,2,3]
        call tensor_contract( 1.0E0_realk,t2,E1,[2],[2],1,1.0E0_realk,omega2,ord)
        ord = [1,2,3,4]
        call tensor_contract(-1.0E0_realk,t2,E2,[4],[1],1,1.0E0_realk,omega2,ord)

        call tensor_ainit(Pijab_om2,omega2%dims,4,tdims=omega2%tdim,atype="TDAR")

        call tensor_lock_local_wins(Pijab_om2,'e',mode)
        call tensor_unlock_local_wins(omega2)

        call tensor_free(E1)
        call tensor_free(E2)

        !INTRODUCE PERMUTATION
        ord = [2,1,4,3]
        call tensor_add(Pijab_om2,1.0E0_realk,omega2, a = 0.0E0_realk, order = ord )
        call tensor_unlock_local_wins(Pijab_om2)
        call tensor_add(omega2,1.0E0_realk,Pijab_om2)

        call lsmpi_barrier(infpar%lg_comm)
        call tensor_free(Pijab_om2)

        omega2%access_type = AT_MASTER_ACCESS
        t2%access_type     = AT_MASTER_ACCESS

#endif
     endif

     call time_start_phase(PHASE_WORK, at = tw)
  end subroutine calculate_E2_and_permute


  subroutine check_job(s,fr,dyn,a,g,na,ng,static,win,prnt)
    implicit none
    logical,intent(in) :: dyn,prnt
    logical,intent(inout) :: fr
    integer,intent(in) :: s,g,na,ng
    integer,intent(inout) :: a
    integer :: static(:)
    integer(kind=ls_mpik) :: win
    real(realk) :: mpi_buf
    integer :: el 
    integer(kind=ls_mpik) :: i, job
#ifdef VAR_MPI
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
           el=1
#ifdef VAR_HAVE_MPI3
           call lsmpi_get_acc(el,a,infpar%master,g,win)
           call lsmpi_win_flush(win,rank = infpar%master, local=.true.)
#else
           call lsmpi_win_lock(infpar%master,win,'e')
           call lsmpi_get_acc(el,a,infpar%master,g,win)
           call lsmpi_win_unlock(infpar%master,win)
#endif
         endif
       endif
       if(fr) fr=.false.
#endif
       a=a+1
       if(a>na)return
#ifdef VAR_MPI
       if(prnt) write (*, '("Rank ",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")') infpar%mynum,&
       &a,na,g,ng
#else
       if(prnt) write (*, '("starting job (",I3,"/",I3,",",I3,"/",I3,")")')a,&
       &na,g,ng
#endif
       call lsmpi_poke()
  end subroutine check_job

  !> \brief Routine to get the c and the d terms from t1 tranformed integrals
  !using a simple mpi-parallelization
  !> \author Patrick Ettenhuber
  !> \Date January 2013 
  subroutine get_cnd_terms_mo_2(w1,w2,w3,t2,u2,govov,gvoov,gvvoo,&
        &no,nv,omega2,s,lock_outside,local)
     implicit none
     !> input some empty workspace of zise v^2*o^2 
     real(realk), intent(inout) :: w1(:)
     real(realk),pointer :: w2(:),w3(:)
     !> the t1-transformed integrals
     type(tensor), intent(inout) :: govov,gvvoo,gvoov
     !> number of occupied orbitals 
     integer, intent(in) :: no
     !> nuber of virtual orbitals
     integer, intent(in) :: nv
     !> ampitudes on input ordered as abij
     !real(realk), intent(in) :: t2(:)
     type(tensor), intent(inout) :: t2
     !> u on input u{aibj}=2t{aibj}-t{ajbi} ordered as abij
     type(tensor), intent(inout) :: u2
     !> the residual to add the contribution
     type(tensor), intent(inout) :: omega2
     !> integer specifying the scheme
     integer, intent(in) :: s
     !> specifiaction if lock stuff
     logical, intent(in) :: lock_outside,local

     !INTERNAL VARIABLES:
     type(tensor) :: Dvoov, Lovov, Coovv, O_pre
     integer :: fdim1(4), sdim1(4), fdim2(4), sdim2(4),ord(4)
     integer :: os, vs
     integer(kind=ls_mpik) :: me, nnod, mode
     integer(kind=8) :: o2v2
     logical :: master
     character(4) :: atype

     call time_start_phase(PHASE_WORK)

     me     = 0_ls_mpik
     nnod   = 1_ls_mpik
#ifdef VAR_MPI
     nnod   = infpar%lg_nodtot
     me     = infpar%lg_mynum
     mode   = MPI_MODE_NOCHECK
     o2v2   = int((i8*no)*no*nv*nv,kind=8)
     master = ( me == 0_ls_mpik )
     os     = govov%tdim(1)
     vs     = govov%tdim(2)
     atype  = "TDAR"

     if(.not.alloc_in_dummy)call tensor_unlock_wins(omega2,.true.)
     !call tensor_lock_local_wins(omega2,'e',mode)

     !Cterm
     fdim1 = [no,no,nv,nv]
     sdim1 = [os,os,vs,vs]
     call tensor_ainit(Coovv,fdim1,4,tdims=sdim1,atype=atype,local=local)
     call tensor_lock_local_wins(Coovv,'e',mode)

     !Build C intermediate
     ord = [2,3,1,4]
     call tensor_add(Coovv,1.0E0_realk,gvvoo, a = 0.0E0_realk, order = ord)
     ord = [3,2,1,4]
     call tensor_contract(-0.5E0_realk,t2,govov,[2,3],[2,3],2,1.0E0_realk,Coovv,ord)

     !Inser synchronizatipn point
     fdim1 = [nv,nv,no,no]
     sdim1 = [vs,vs,os,os]
     call tensor_ainit(O_pre,fdim1,4,tdims=sdim1,atype=atype,local=local,fo = omega2%offset)
     call tensor_lock_local_wins(O_pre,'e',mode)

     !now allow for access to the completed tiles
     call tensor_unlock_local_wins(Coovv)

     ord = [4,1,2,3]
     call tensor_contract(-1.0E0_realk,t2,Coovv,[2,3],[4,1],2,0.0E0_realk,O_pre,ord)
     
     !synchronize
     call tensor_free(Coovv)

     call tensor_unlock_local_wins(O_pre)

     !add in permutations (1+0.5P_ij)
     call tensor_add(omega2,1.0E0_realk,O_pre)
     ord = [1,2,4,3]
     call tensor_add(omega2,0.5E0_realk,O_pre,order=ord)

     !synchronize
     call tensor_free(O_pre)

     !Dterm
     !Calculate intermediates needed in D2 term
     fdim1  = [nv,no,no,nv]
     sdim1  = [vs,os,os,vs]
     fdim2  = govov%dims
     sdim2  = govov%tdim
     call tensor_ainit(Dvoov,fdim1,4,tdims=sdim1,atype=atype,local=local)
     call tensor_ainit(Lovov,fdim2,4,tdims=sdim2,atype=atype,local=local)
     call tensor_lock_local_wins(Dvoov,'e',mode)
     call tensor_lock_local_wins(Lovov,'e',mode)

     !careful gvvoo is ordered as (aijb) and gvoov is ordered as (ajbi) 
     ord = [1,4,2,3]
     call tensor_add(Dvoov, 2.0E0_realk,gvoov, a = 0.0E0_realk, order = ord)
     ord = [1,3,2,4]
     call tensor_add(Dvoov,-1.0E0_realk,gvvoo, order = ord)

     ord = [1,4,3,2]
     call tensor_add(Lovov, 2.0E0_realk,govov, a = 0.0E0_realk)
     call tensor_add(Lovov,-1.0E0_realk,govov,order = ord )
     call tensor_unlock_local_wins(Lovov)

     !u2 is saved as (baij) 
     ord = [1,2,3,4]
     call tensor_contract(0.5E0_realk,u2,Lovov,[4,1],[1,2],2,1.0E0_realk,Dvoov,ord)

     !Inser synchronizatipn point
     fdim1 = [nv,nv,no,no]
     sdim1 = [vs,vs,os,os]
     call tensor_ainit(O_pre,fdim1,4,tdims=sdim1,atype=atype,local=local,fo = omega2%offset)
     call tensor_lock_local_wins(O_pre,'e',mode)

     call tensor_unlock_local_wins(Dvoov)

     !synchronization point
     call tensor_free(Lovov)

     !u2 is saved as (baij) 
     ord = [3,1,4,2]
     call tensor_contract(0.5E0_realk,u2,Dvoov,[1,4],[4,3],2,0.0E0_realk,O_pre,ord)

     !synchronization point
     call tensor_free(Dvoov)

     call tensor_unlock_local_wins(O_pre)

     !add in permutations P_ij^ab (1+0.5P_ij)
     call tensor_add(omega2,1.0E0_realk,O_pre)


     !call tensor_unlock_wins(omega2,.true.)
     call tensor_free(O_pre)

#else
     call lsquit("ERROR(get_cnd_terms_mo_2): MPI only routine",-1)
#endif

  end subroutine get_cnd_terms_mo_2
  
  !> \brief Routine to get the c and the d terms from t1 tranformed integrals
  !using a simple mpi-parallelization, only for schemes 4 and 3, 2 is also
  !possible, but some functions seem to be buggy.
  !> \author Patrick Ettenhuber
  !> \Date January 2013 
  subroutine get_cnd_terms_mo_3n4(w1,w2,w3,t2,u2,govov,gvoov,gvvoo,&
        &no,nv,omega2,s,lock_outside,els2add,tw,tc)
     implicit none
     !> input some empty workspace of zise v^2*o^2 
     real(realk), intent(inout) :: w1(:)
     real(realk),pointer :: w2(:),w3(:)
     !> the t1-transformed integrals
     type(tensor), intent(inout) :: govov,gvvoo,gvoov
     !> number of occupied orbitals 
     integer, intent(in) :: no
     !> nuber of virtual orbitals
     integer, intent(in) :: nv
     !> ampitudes on input ordered as abij
     !real(realk), intent(in) :: t2(:)
     type(tensor), intent(inout) :: t2
     !> u on input u{aibj}=2t{aibj}-t{ajbi} ordered as abij
     type(tensor), intent(inout) :: u2
     !> the residual to add the contribution
     type(tensor), intent(inout) :: omega2
     !> integer specifying the scheme
     integer, intent(in) :: s
     !> specifiaction if lock stuff
     logical, intent(in) :: lock_outside
     !> specify how many elements can be added to w3 buffer
     integer(kind=8),intent(in) :: els2add
     !> timing information
     real(realk),intent(inout) :: tw, tc
     integer(kind=ls_mpik) :: nnod, me
     integer :: tl,fai,lai,i,faif,lead
     integer :: l,ml
     integer(kind=ls_mpik) :: nod,mode
     real(realk) :: nrm1,nrm2,nrm3,nrm4,t_comm_b
     integer :: a,b,j,fri,tri
     integer(kind=8) :: o2v2,tlov,w1size,w2size,w3size
     character(tensor_MSG_LEN) :: msg
     real(realk) :: MemFree, startt, stopp
     logical :: traf,trafi,master

     call time_start_phase(PHASE_WORK)

     me     = 0_ls_mpik
     nnod   = 1_ls_mpik
#ifdef VAR_MPI
     nnod   = infpar%lg_nodtot
     me     = infpar%lg_mynum
     mode   = MPI_MODE_NOCHECK
#endif
     o2v2   = int((i8*no)*no*nv*nv,kind=8)
     w1size = o2v2
     master = ( me == 0_ls_mpik )

     !Setting transformation variables for each rank
     !**********************************************
     call mo_work_dist(nv*no,fai,tl,traf)

     tlov  = int((i8*tl)*no*nv,kind=8)

     if(DECinfo%PL>3.and.me==0)then
        write(DECinfo%output,'("Trafolength in striped CD:",I5)')tl
     endif

     ! Go through that only if the transformation has to be performed
     if(traf.or.me==0)then

        if(s==4)then
           faif = fai
           lead = no * nv
           w2size = o2v2
           w3size = o2v2
        else if(s==3.or.s==2)then
           faif = 1
           lead = tl
           !use w3 as buffer which is allocated largest possible
           w2size  = tlov
           w3size  = min(o2v2,tlov + els2add)
        else
           call lsquit("ERROR(get_cnd_terms_mo_3n4):no valid scheme",-1)
        endif

        if(me==0.and.DECinfo%PL>3)then
           print *,"w2size(2)",w2size
           print *,"w3size(2)",w3size
        endif

        call mem_alloc(w2,w2size)
        call mem_alloc(w3,w3size)


        !calculate doubles C term
        !*************************

        !Reorder gvvoo [a c k i] -> goovv [a i c k]
        if(s==4)then
           call array_reorder_4d(1.0E0_realk,gvvoo%elm1,nv,no,no,nv,[1,3,4,2],0.0E0_realk,w2)
        else if(s==3)then
           call array_reorder_4d(1.0E0_realk,gvvoo%elm1,nv,no,no,nv,[1,3,4,2],0.0E0_realk,w1)
           do i=1,tl
              call dcopy(no*nv,w1(fai+i-1),no*nv,w2(i),tl)
           enddo
        endif

        !Reorder t [a d l i] -> t [a i d l]
        if(s==4)then
           call array_reorder_4d(1.0E0_realk,t2%elm1,nv,nv,no,no,[1,4,2,3],0.0E0_realk,w3)
           !w3 = 0.0E0_realk
           !call array_reorder_4d(1.0E0_realk,t2%elm1,nv,nv,no,no,[1,4,2,3],0.0E0_realk,w1)
           !do i=1,tl
           !   call dcopy(no*nv,w1(fai+i-1),no*nv,w3(fai+i-1),no*nv)
           !enddo
        else if(s==3)then
           call array_reorder_4d(1.0E0_realk,t2%elm1,nv,nv,no,no,[1,4,2,3],0.0E0_realk,w1)
           do i=1,tl
              call dcopy(no*nv,w1(fai+i-1),no*nv,w3(i),tl)
           enddo
        endif

        !stop 0
        !SCHEME 4 AND 3 because of w1 being buffer before
        !Reorder govov [k d l c] -> govov [d l c k]
        if(s==3.or.s==4)then
           call array_reorder_4d(1.0E0_realk,govov%elm1,no,nv,no,nv,[2,3,4,1],0.0E0_realk,w1)
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       CENTRAL GEMM 1         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !(-0.5) * t [a i d l] * govov [d l c k] + goovv [a i c k] = C [a i c k]
        call dgemm('n','n',tl,no*nv,no*nv,-0.5E0_realk,w3(faif),lead,w1,no*nv,1.0E0_realk,w2(faif),lead)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       CENTRAL GEMM 2         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !(-1) * C [a i c k] * t [c k b j] = preOmC [a i b j]
        if(s==4)then
           w1=0.0E0_realk
           call dgemm('n','t',tl,no*nv,no*nv,-1.0E0_realk,w2(faif),lead,w3,no*nv,0.0E0_realk,w1(fai),no*nv)
        else if(s==3)then
           call array_reorder_4d(1.0E0_realk,t2%elm1,nv,nv,no,no,[1,4,2,3],0.0E0_realk,w1)
           call dgemm('n','t',tl,no*nv,no*nv,-1.0E0_realk,w2(faif),lead,w1,no*nv,0.0E0_realk,w3,lead)
           w1=0.0E0_realk
           do i=1,tl
              call dcopy(no*nv,w3(i),tl,w1(fai+i-1),no*nv)
           enddo
           !print *,infpar%lg_mynum,"DGEMM2 -- out",norm2(w1),norm2(w3)
        endif


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       Omega update           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(s==3.or.s==4)then
           !contribution 1: 0.5*preOmC [a i b j] -> =+ Omega [a b i j]
           call array_reorder_4d(0.5E0_realk,w1,nv,no,nv,no,[1,3,2,4],1.0E0_realk,omega2%elm1)
           !contribution 3: preOmC [a j b i] -> =+ Omega [a b i j]
           call array_reorder_4d(1.0E0_realk,w1,nv,no,nv,no,[1,3,4,2],1.0E0_realk,omega2%elm1)
        endif


        !print *,infpar%lg_mynum,"done C"
        !call lsmpi_barrier(infpar%lg_comm)
        !call print_norm(omega2)



        !calculate doubles D term
        !************************
        !(-1) * gvvoo [a c k i] -> + 2*gvoov[a i k c] = L [a i k c]

        if(s==4)then
           call array_reorder_4d(2.0E0_realk,gvoov%elm1,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w2)
           call array_reorder_4d(-1.0E0_realk,gvvoo%elm1,nv,no,no,nv,[1,3,2,4],1.0E0_realk,w2)
        else if(s==3)then
           call array_reorder_4d(2.0E0_realk,gvoov%elm1,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w1)
           call array_reorder_4d(-1.0E0_realk,gvvoo%elm1,nv,no,no,nv,[1,3,2,4],1.0E0_realk,w1)
           do i=1,tl
              call dcopy(no*nv,w1(fai+i-1),no*nv,w2(i),tl)
           enddo
        endif


        !Transpose u [d a i l] -> u [a i l d]
        if(s==4)then
           call array_reorder_4d(1.0E0_realk,u2%elm1,nv,nv,no,no,[2,3,4,1],0.0E0_realk,w3)
        else if(s==3)then
           call array_reorder_4d(1.0E0_realk,u2%elm1,nv,nv,no,no,[2,3,4,1],0.0E0_realk,w1)
           do i=1,tl
              call dcopy(no*nv,w1(fai+i-1),no*nv,w3(i),tl)
           enddo
        endif

        !SCHEME 3 AND 4, because of the reordering using w1
        !(-1) * govov [l c k d] + 2*govov[l d k c] = L [l d k c]
        if(s==3.or.s==4)then
           call array_reorder_4d(2.0E0_realk,govov%elm1,no,nv,no,nv,[1,2,3,4],0.0E0_realk,w1)
           call array_reorder_4d(-1.0E0_realk,govov%elm1,no,nv,no,nv,[1,4,3,2],1.0E0_realk,w1)
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       CENTRAL GEMM 1         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! (0.5) * u [a i l d] * L [l d k c] + L [a i k c] = D [a i k c]
        call dgemm('n','n',tl,nv*no,nv*no,0.5E0_realk,w3(faif),lead,w1,nv*no,1.0E0_realk,w2(faif),lead)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       CENTRAL GEMM 2         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! (0.5)*D[a i k c] * u [b j k c]^T  = preOmD [a i b j]
        if(s==4)then
           w1=0.0E0_realk
           call dgemm('n','t',tl,nv*no,nv*no,0.5E0_realk,w2(faif),lead,w3,nv*no,0.0E0_realk,w1(fai),nv*no)
        else if(s==3)then
           call array_reorder_4d(1.0E0_realk,u2%elm1,nv,nv,no,no,[2,3,4,1],0.0E0_realk,w1)
           call dgemm('n','t',tl,nv*no,nv*no,0.5E0_realk,w2(faif),lead,w1,nv*no,0.0E0_realk,w3,lead)
           w1=0.0E0_realk
           do i=1,tl
              call dcopy(no*nv,w3(i),tl,w1(fai+i-1),no*nv)
           enddo
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       Omega update           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! preOmD [a i b j] -> =+ Omega [a b i j]
        if(s==4.or.s==3)then
           call array_reorder_4d(1.0E0_realk,w1,nv,no,nv,no,[1,3,2,4],1.0E0_realk,omega2%elm1)
        endif

        call mem_dealloc(w2)
        call mem_dealloc(w3)

     endif

     !print *,infpar%lg_mynum,"done D"
     !call lsmpi_barrier(infpar%lg_comm)
     !call print_norm(omega2)

     call time_start_phase(PHASE_WORK, at = tw )
  end subroutine get_cnd_terms_mo_3n4




  !> Purpose: wrapper for batch size determination routines
  !           in CCSD and MO-CCSD algorithms.
  !
  !> Author:  Pablo Baudin 
  !> Date:    March 2014
  subroutine wrapper_get_ccsd_batch_sizes(MyFragment,bat,mpi_split,ntasks)

    implicit none

    !> Atomic fragment
    type(decfrag), intent(inout) :: MyFragment
    !> AO batch information
    type(mp2_batch_construction), intent(inout) :: bat
    !> return number of tasks (only for MO-CCSD)
    integer, intent(inout) :: ntasks

    real(realk) :: MemFree
    integer :: scheme, nbas, nocc, nvir, MinAObatch, iter
    integer :: dimMO, nMObatch, ntot,os,vs
    integer(kind=8) :: dummy
    logical :: mo_ccsd, local_moccsd, mpi_split
    integer :: bs

    ! For fragment with local orbitals where we really want to use the fragment-adapted orbitals
    ! we need to set nocc and nvirt equal to the fragment-adapted dimensions
    nocc = MyFragment%noccAOS
    nvir = MyFragment%nunoccAOS 
    bs   = get_split_scheme_0(MyFragment%nbasis)

    ! For MO-CCSD part
    ntot    = nocc + nvir
    nbas    = MyFragment%nbasis
    mo_ccsd = .true.
    if (DECinfo%NO_MO_CCSD.or.(nbas>400)) mo_ccsd = .false.

    if (DECinfo%force_scheme) then
      scheme=DECinfo%en_mem
      if (scheme<5) then
        DECinfo%NO_MO_CCSD = .true.
        mo_ccsd            = .false.
      else if (scheme>=5) then 
        mo_ccsd            = .true.
        if (DECinfo%NO_MO_CCSD) call lsquit('ERROR(CCSD): Inconsistent input, CCSD schemes &
           & 5 and 6 require the MO based algorithm. (Remove NO_MO_CCSD keyword)', DECinfo%output)
      end if
    end if

        
#ifdef MOD_UNRELEASED
    ! The two if statments are necessary as mo_ccsd might become false
    ! after the first statement (if not enought memory).
    if (mo_ccsd) then
      call get_MO_and_AO_batches_size(mo_ccsd,local_moccsd,ntot,nbas,nocc,nvir, &
           & dimMO,nMObatch,bat%MaxAllowedDimAlpha,bat%MaxAllowedDimGamma, &
           & MyFragment%MyLsItem,mpi_split)
      ntasks = nMObatch*(nMObatch+1)/2
    end if
#endif

    if (.not.mo_ccsd) then 
      iter=1
      call get_symm_tensor_segmenting_simple(nocc,nvir,os,vs)
      call determine_maxBatchOrbitalsize(DECinfo%output,MyFragment%MyLsItem%setting,MinAObatch,'R')
      call get_currently_available_memory(MemFree)
      call get_max_batch_sizes(scheme,MyFragment%nbasis,bs,nvir,vs,nocc,os,bat%MaxAllowedDimAlpha, &
           & bat%MaxAllowedDimGamma,MinAObatch,DECinfo%manual_batchsizes,iter,MemFree, &
           & .true.,dummy,(.not.DECinfo%solver_par),mpi_split,MyFragment%mylsitem%setting,['R','R','R','R','C'])
    end if

  end subroutine wrapper_get_ccsd_batch_sizes

  !> \brief calculate batch sizes automatically-->dirty but better than nothing
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine get_max_batch_sizes(scheme,nb,bs,nv,vs,no,os,nba,nbg,&
  &minbsize,manual,iter,MemFree,first,e2a,local,mpi_split,se,is,nbuf)
    implicit none
    integer, intent(inout) :: scheme
    integer, intent(in)    :: nb,bs,nv,vs,no,os
    integer :: iter
    integer, intent(inout) :: nba,nbg,minbsize
    real(realK),intent(in) :: MemFree
    real(realk)            :: mem_used,frac_of_total_mem,m
    logical,intent(in)     :: manual,first
    type(lssetting),intent(inout) :: se
    Character,intent(in) :: is(5)
    integer(kind=8), intent(inout) :: e2a
    logical, intent(in)    :: local, mpi_split
    integer, intent(out),optional :: nbuf
    integer(kind=8) :: v2o2,thrsize,w0size,w1size,w2size,w3size
    integer :: nnod,magic,nbuffs

    frac_of_total_mem=0.80E0_realk
    nba=minbsize
    nbg=minbsize
    nnod=1
    e2a  = 0
    v2o2 = (i8*no*no)*nv*nv
#ifdef VAR_MPI
    nnod=infpar%lg_nodtot
#endif
    !magic = DECinfo%MPIsplit/5
    magic = 2
    !test for scheme with highest reqirements --> fastest
    scheme=4
    mem_used=get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,4,scheme,.false.,se,is)
    if (mem_used>frac_of_total_mem*MemFree)then
#ifdef VAR_MPI
       !test for scheme with medium requirements
       scheme=3
       mem_used=get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,4,scheme,.false.,se,is)
       if (mem_used>frac_of_total_mem*MemFree)then
          !test for scheme with low requirements
          scheme=2
          mem_used=get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,4,scheme,.false.,se,is)
          if (mem_used>frac_of_total_mem*MemFree)then
             scheme=0
             mem_used=get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,4,scheme,.false.,se,is)
             frac_of_total_mem=0.60E0_realk
             print *,"mem",mem_used,frac_of_total_mem,MemFree
             if (mem_used>frac_of_total_mem*MemFree)then
                write(DECinfo%output,*) "MINIMUM MEMORY REQUIREMENT IS NOT AVAILABLE"
                write(DECinfo%output,'("Fraction of free mem to be used:          ",f8.3," GB")')&
                   &frac_of_total_mem*MemFree
                write(DECinfo%output,'("Memory required in memory saving scheme:  ",f8.3," GB")')mem_used
                mem_used=get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,4,3,.false.,se,is)
                write(DECinfo%output,'("Memory required in intermediate scheme: ",f8.3," GB")')mem_used
                mem_used=get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,4,4,.false.,se,is)
                write(DECinfo%output,'("Memory required in memory wasting scheme: ",f8.3," GB")')mem_used
                call lsquit("ERROR(CCSD): there is just not enough memory&
                   & available",DECinfo%output)
             endif
          endif
       endif
#endif
    endif

    if(DECinfo%force_scheme)then
      scheme=DECinfo%en_mem
      print *,"!!FORCING CCSD!!"
      if(local)then
        if(scheme==3.or.scheme==2.or.scheme==1)then
          print *,"CHOSEN SCHEME DOES NOT WORK WITHOUT PARALLEL SOLVER, USE&
          & MORE THAN ONE NODE"
          call lsquit("ERROR(get_ccsd_residual_integral_driven):invalid scheme",-1)
        endif
      endif
      if(scheme==4)then
        print *,"SCHEME 4: NON PDM-SCHEME WITH HIGH MEMORY REQUIREMENTS"
      else if(scheme==3)then
        print *,"SCHEME 3: WITH MEDIUM MEMORY REQUIREMENTS (PDM)"
      else if(scheme==2)then
        print *,"SCHEME 2: WITH LOW MEMORY REQUIREMENTS (PDM)"
      else if(scheme==1)then
        print *,"SCHEME 1: DMITRY's SCHEME (PDM)"
      else if(scheme==0)then
        print *,"SCHEME 0: COMPLETE PDM"
      else
        print *,"SCHEME ",scheme," DOES NOT EXIST"
        call lsquit("ERROR(get_ccsd_residual_integral_driven):invalid scheme2",-1)
      endif
    endif


    ! Attention this manual block should be only used for debugging, also you
    ! will have to ajust its properties to your system. The block is called
    ! with the keywordk ".manual_batchsizes" in LSDALTON.INP the next line
    ! then contains the alpha and gamma batch sizes separated by space
    if (manual) then
      ! KK and PE hacks -> only for debugging
      ! extended to mimic the behaviour of the mem estimation routine when memory is filled up
      !if((DECinfo%ccsdGbatch==0).and.(DECinfo%ccsdAbatch==0)) then
      !  call get_max_batch_sizes(scheme,nb,nv,no,nba,nbg,minbsize,.false.,iter,MemFree, &
      !       & .false.,e2a,local,mpi_split)
      !else
        nba = DECinfo%ccsdAbatch - iter * 0
        nbg = DECinfo%ccsdGbatch - iter * 0
      !endif
      ! Use value given in input --> the zero can be adjusted to vary batch sizes during the iterations

      m = nbg-0*iter
      if( minbsize <  m ) nbg = m
      if( minbsize >= m ) nbg = minbsize

      m = nba-0*iter
      if( minbsize <  m ) nba = m
      if( minbsize >= m ) nba = minbsize

      if( nbg>=nb )       nbg = nb
      if( nba>=nb )       nba = nb

      mem_used=get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,4,scheme,.false.,se,is)

      if (frac_of_total_mem*MemFree<mem_used) then
        print *, "ATTENTION your chosen batch sizes might be too large!!!"
      endif

    ! This block is routinely used in a calculation
    else

      !optimize buffers
      nbuffs = 0
      if(scheme==0)then
         do nbuffs=11,5
            mem_used = get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,nbuffs,iter,4,scheme,.false.,se,is)
            if (frac_of_total_mem*MemFree>mem_used) exit
         enddo
      endif
      if(present(nbuf))nbuf=nbuffs

      !make batches larger until they do not fit anymore
      !determine gamma batch first
      do while ((frac_of_total_mem*MemFree>mem_used) .and. (nb>=nbg))

        nbg=nbg+1
        mem_used=get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,3,scheme,.false.,se,is)

      enddo

      if (nbg>=nb)then
        nbg = nb
      else if (nbg<=minbsize)then
        nbg = minbsize
      else
        nbg=nbg-1
      endif

      !determine alpha batch

      do while ((frac_of_total_mem*MemFree>mem_used) .and. (nb>=nba))
        nba      = nba+1
        mem_used = get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,3,scheme,.false.,se,is)
      enddo

      if (nba>=nb)then
         nba = nb
      else if (nba<=minbsize)then
         nba = minbsize
      else
         nba=nba-1
      endif

    endif
    mem_used=get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,4,scheme,.false.,se,is)

    ! mpi_split should be true when we want to estimate the workload associated
    ! to a DEC fragment and eventually split the slots. In this case, the next
    ! step must be skiped.
    if (.not.mpi_split.and..not.scheme==0) then
      !if much more slaves than jobs are available, split the jobs to get at least
      !one for all the slaves
      !print *,"JOB SPLITTING WITH THE NUMBER OF NODES HAS BEEN DEACTIVATED"
      if(.not.manual)then
        if((nb/nba)*(nb/nbg)<magic*nnod.and.(nba>minbsize).and.nnod>1)then
          nba=(nb/(magic*nnod))
          if(nba<minbsize)nba=minbsize
        endif
       
        if((nb/nba)*(nb/nbg)<magic*nnod.and.(nba==minbsize).and.nnod>1)then
          do while((nb/nba)*(nb/nbg)<magic*nnod)
            nbg=nbg-1
            if(nbg<1)exit
          enddo
          if(nbg<minbsize)nbg=minbsize
        endif
      endif
    end if

    if(scheme==2)then
      mem_used = get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,iter,2,scheme,.false.,se,is)
      e2a = min(v2o2,int(((frac_of_total_mem*MemFree - mem_used)*1E9_realk*0.5E0_realk/8E0_realk),kind=8))
    endif
  end subroutine get_max_batch_sizes

!> \brief calculate the memory requirement for the matrices in the ccsd routine
!> \author Patrick Ettenhuber
!> \date January 2012
  function get_min_mem_req(no,os,nv,vs,nb,bs,nba,nbg,nbuffs,iter,choice,s,print_stuff,se,is) result (memrq)
    implicit none
    integer, intent(in) :: no,os,nv,vs,nb,bs
    integer, intent(in) :: nba,nbg,nbuffs
    integer, intent(in) :: iter,choice
    real(realk) :: memrq, memin, memout
    logical, intent(in) :: print_stuff
    integer, intent(in) :: s
    type(lssetting),intent(inout) :: se
    Character,intent(in) :: is(5)
    integer :: nor, nvr, fe, ne , nnod, me, i
    integer :: d1(4),ntpm(4),tdim(4),mode,splt,ntiles,tsze
    integer(kind=ls_mpik) :: master
    integer :: l1,ml1,fai1,tl1
    integer :: l2,ml2,fai2,tl2
    integer :: l3,ml3,fai3,tl3
    integer :: l4,ml4,fai4,tl4
    integer :: nloctiles
    integer :: cd , e2
    integer(kind=long) :: w0size, w1size, w2size, w3size
    nor = no*(no+1)/2
    nvr = nv*(nv+1)/2
    master = 0
    mode=4
    d1(1)=nv;d1(2)=nv;d1(3)=no;d1(4)=no
    nnod = 1
    me = 0
#ifdef VAR_MPI
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

    w0size = get_wsize_for_ccsd_int_direct(0,no,os,nv,vs,nb,bs,nba,nbg,nbuffs,s,se,is)
    w1size = get_wsize_for_ccsd_int_direct(1,no,os,nv,vs,nb,bs,nba,nbg,nbuffs,s,se,is)
    w2size = get_wsize_for_ccsd_int_direct(2,no,os,nv,vs,nb,0,nba,nbg,nbuffs,s,se,is)
    w3size = get_wsize_for_ccsd_int_direct(3,no,os,nv,vs,nb,0,nba,nbg,nbuffs,s,se,is)
    !w0
    memin = 1.0E0_realk * w0size
    !w1
    memin = memin + 1.0E0_realk * w1size
    !w2
    memin = memin + 1.0E0_realk * w2size
    !w3
    memin = memin + 1.0E0_realk * w3size

    !calculate minimum memory requirement
    ! u+3*integrals
    select case(s)
    case(4)

      !THROUGHOUT THE ALGORITHM
      !************************

      ! u 2 + Omega 2 +  H +G  
      memrq = 1.0E0_realk*(2_long*no*no*nv*nv+ i8*nb*nv+i8*nb*no)
      !gvoov gvvoo
      memrq=memrq+ 2.0E0_realk*(i8*nv*nv)*no*no


      !INSIDE OF MAIN LOOP
      !*******************

      !uigcj sio4
      memin = memin +1.0E0_realk*((i8*no*no)*nv*nbg+(i8*no*no)*nor)
      !tpl tmi
      memin = memin + (i8*nor)*nvr*2.0E0_realk


      !OUTSIDE OF MAIN LOOP
      !********************

      ! w1 + FO + w2 + w3
      memout = 1.0E0_realk*(max((i8*nv*nv)*no*no,i8*nb*nb)+i8*nb*nb+(2_long*no*no)*nv*nv)
      ! govov
      memout = memout + (1.0E0_realk*no*no)*nv*nv

    case(3)



      !THROUGHOUT THE ALGORITHM
      !************************

      !govov stays in pdm and is dense in second part
      ! u 2 + omega2 + H +G  + keep space for one update batch
      call get_int_dist_info(int((i8*nv)*nv*no*no,kind=8),fe,ne,master)
      memrq = 1.0E0_realk*((2_long*no*no)*nv*nv+ i8*nb*nv+i8*nb*no) + i8*ne
      !gvoov gvvoo (govov, allocd outside)
      memrq=memrq+ 2.0E0_realk*ne 


      !INSIDE OF MAIN LOOP
      !*******************

      !uigcj sio4
      memin = memin + 1.0E0_realk*((i8*no*no)*nv*nbg+(i8*no*no)*nor)
      !tpl tmi
      memin = memin + 1.0E0_realk*nor*(nvr*2_long)

      !OUTSIDE OF MAIN LOOP
      !********************

      ! w1 + FO + w2 + w3 + govov + full gvvoo + full gvoov
      memout = 1.0E0_realk*(max((i8*nv*nv)*no*no,i8*nb*nb) &
           & + max(i8*nb*nb,max(2_long*tl1,i8*tl2)))       &
           & + 2.0E0_realk * (i8*nv**2)*no**2

    case(2)

       !TODO: ADAPT TO ACTUAL REQUIREMENTS

      !THROUGHOUT THE ALGORITHM
      !************************

      call tensor_default_batches(d1,mode,tdim,splt)
      call tensor_get_ntpm(d1,tdim,mode,ntpm,ntiles)
      nloctiles=ceiling(float(ntiles)/float(nnod))
      tsze = 1
      do i = 1, mode
        tsze = tsze * tdim(i)
      enddo
      !govov stays in pdm and is dense in second part
      ! u2 + H +G + space for 2 update tile s
      memrq = 1.0E0_realk*((i8*tsze)*nloctiles+ i8*nb*nv+i8*nb*no + i8*2*tsze)
      !gvoov gvvoo
      memrq=memrq+ 2.0E0_realk*tsze*nloctiles


      !INSIDE OF MAIN LOOP
      !*******************

      !uigcj sio4
      memin = memin + 1.0E0_realk*((i8*no*no)*nv*nbg+(i8*no*no)*nor)
      !tpl tmi
      memin = memin + 1.0E0_realk*(nor*nvr*i8)

      !OUTSIDE OF MAIN LOOP
      !********************

      ! w1 + FO + w2 + w3
      !in cd terms w2 and w3 have tl1, in b2 w2 has tl2
      cd = max(2_long*tl1,i8*tl2)
      ! in e2 term w2 has max(tl2,tl3) and w3 has max(no2,nv2)
      e2 = max(tl3,tl4) + max(no*no,nv*nv)

      memout = 1.0E0_realk*(max((i8*nv*nv)*no*no,(i8*nb*nb))+max(i8*nb*nb,i8*max(cd,e2)))

    case(1)

       print *,"Dmitry, please implement your memory requirements here, such&
       & that a memory estimation can be made and the batch sizes adapted -- PE"

    case(0)

       !RIGHT NOW THIS IS SET TO THE SAME VALUE AS SCHEME 4 JUST /nnod
       ! u 2 + Omega 2 +  H +G  
       memrq = 1.0E0_realk*(2_long*no*no*nv*nv+ i8*nb*nv+i8*nb*no)
       !gvoov gvvoo
       memrq=memrq+ 2.0E0_realk*(i8*nv*nv)*no*no

       memrq = memrq/float(nnod)

       !INSIDE OF MAIN LOOP
       !*******************

       !uigcj sio4
       memin = memin +1.0E0_realk*((i8*no*no)*nv*nbg+(i8*no*no)*nor)/float(nnod)
       !tpl tmi
       memin = memin + (i8*nor)*nvr*2.0E0_realk/float(nnod)

       !OUTSIDE OF MAIN LOOP
       !********************
       memout = 1.0E0_realk*max(max(max(max(i8*nb*nb,i8*nb*nv),i8*nb*no),i8*nv*nv),i8*no*no)

    case default

      print *,"DECinfo%force_scheme",DECinfo%force_scheme,s
      call lsquit("ERROR(get_min_mem_req):requested memory scheme not known,&
      & should not happen",-1)

    end select

    if(print_stuff) then
      write(DECinfo%output,*) "Memory requirements:"
      write(DECinfo%output,*) "Basic  :",(memrq *8.0E0_realk)/(1.024E3_realk**3)
      write(DECinfo%output,*) "Part B :",(memin *8.0E0_realk)/(1.024E3_realk**3)
      write(DECinfo%output,*) "Part C :",(memout*8.0E0_realk)/(1.024E3_realk**3)
    endif


    !SELECTOR OF WHICH INFO TO RETURN
    !********************************
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
!#ifdef VAR_MPI
!    if(LSMPIASYNCP)then
!       memrq = 1.5*memrq
!    endif
!#endif

  end function get_min_mem_req



  !> \brief Precondition singles 


  function precondition_singles_newarr(omega1,ppfock,qqfock) result(prec)

    implicit none
    type(tensor), intent(in) :: omega1,ppfock,qqfock
    type(tensor) :: prec
    integer, dimension(2) :: dims
    integer :: a,i
    if(omega1%mode/=2.or.ppfock%mode/=2.or.qqfock%mode/=2)then
      call lsquit("ERROR(precondition_singles_newarr):wrong number of modes&
      & for this operation",DECinfo%output)
    endif

    dims = omega1%dims
    call tensor_init(prec, dims,2)

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
       & Co,Cv,Co2,Cv2)

    implicit none
    type(array2), intent(inout) :: t1
    type(array2), intent(inout) :: Co,Cv,Co2,Cv2
    type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt

    ! Occupied X and Y matrices
    call getT1transformation_occ(t1,xocc,yocc,Co,Co2,Cv2)

    ! Virtual X and Y matrices
    call getT1transformation_virt(t1,xvirt,yvirt,Co,Cv,Cv2)


  end subroutine getT1transformation


  !> \brief T1 transformations for occupied X and Y matrices
  subroutine getT1transformation_occ(t1,xocc,yocc,Co,Co2,Cv2)

    implicit none
    type(array2), intent(inout) :: t1
    type(array2), intent(inout) :: Co,Co2,Cv2
    type(array2), intent(inout) :: xocc,yocc

    ! xocc
    call array2_zero(xocc)
    call array2_copy(xocc,Co)

    ! yocc
    call array2_zero(yocc)
    call array2_copy(yocc,Co2)
    call array2_matmul(Cv2,t1,yocc,'n','n',1.0E0_realk,1.0E0_realk)


  end subroutine getT1transformation_occ



  !> \brief T1 transformations for virtual X and Y matrices
  subroutine getT1transformation_virt(t1,xvirt,yvirt,Co,Cv,Cv2)

    implicit none
    type(array2), intent(inout) :: t1
    type(array2), intent(inout) :: Co,Cv,Cv2
    type(array2), intent(inout) :: xvirt,yvirt

    ! xvirt
    call array2_zero(xvirt)
    call array2_copy(xvirt,Cv)
    call array2_matmul(Co,t1,xvirt,'n','t',-1.0E0_realk,1.0E0_realk)

    ! yvirt
    call array2_zero(yvirt)
    call array2_copy(yvirt,Cv2)

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
    type(array2) :: Co,Co2,Cv
    integer :: nbasis,nocc,nvirt
    integer,dimension(2) :: bo,bv

    ! Dimensions
    nbasis = MyMolecule%nbasis
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nunocc
    bo(1) = nbasis
    bo(2) = nocc
    bv(1) = nbasis
    bv(2) = nvirt

    ! Initialize stuff in array2 format
    Co = array2_init(bo,MyMolecule%Co)
    Co2 = array2_init(bo,MyMolecule%Co)
    Cv = array2_init(bv,MyMolecule%Cv)

    ! Get T1-transformed Fock matrix in AO basis
    call Get_AOt1Fock(mylsitem,t1,fockt1,nocc,nvirt,nbasis,Co,Co2,Cv)

    ! Free stuff
    call array2_free(Co)
    call array2_free(Co2)
    call array2_free(Cv)


  end subroutine fullmolecular_get_AOt1Fock



  !> \brief Get T1 transformed Fock matrix in the AO basis using
  !> the fragment T1 amplitudes (dimension: occ AOS, virt AOS) and fragment MOs.
  subroutine fragment_get_AOt1Fock(MyFragment,fockt1)
    implicit none
    !> Fragment info
    type(decfrag), intent(inout) :: MyFragment
    !> T1-transformed Fock matrix in the AO basis (also initialized here)
    type(array2),intent(inout) :: fockt1
    type(array2) :: t1
    type(array2) :: Co,Co2,Cv
    integer :: nbasis,nocc,nvirt
    integer,dimension(2) :: bo,bv,vo

    ! Dimensions
    ! **********
    nbasis = MyFragment%nbasis
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
    Co = array2_init(bo,MyFragment%Co)
    Co2 = array2_init(bo,MyFragment%Co)   ! particle and hole coefficients are identical
    Cv = array2_init(bv,MyFragment%Cv)

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
    call Get_AOt1Fock(MyFragment%mylsitem,t1,fockt1,nocc,nvirt,nbasis,Co,Co2,Cv)

    ! Free stuff
    call array2_free(Co)
    call array2_free(Co2)
    call array2_free(Cv)


  end subroutine fragment_get_AOt1Fock


  subroutine Get_AOt1Fock_arraywrapper(mylsitem,t1,fockt1,nocc,nvirt,nbasis,Co,Co2,Cv)
    implicit none
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Singles amplitudes for fragment or full molecule
    !> Although this is intent(inout), it is unchanged at output
    type(tensor),intent(inout) :: t1
    !> T1-transformed Fock matrix in the AO basis (also initialized here)
    type(tensor),intent(inout) :: fockt1
    !> Number of occupied orbitals (fragment AOS or full molecule)
    integer,intent(in) :: nocc
    !> Number of virtual orbitals (fragment AOS or full molecule)
    integer,intent(in) :: nvirt
    !> Number of basis functions (atomic fragment extent or full molecule)
    integer,intent(in) :: nbasis
    !> Occupied MOs (particle)
    type(tensor),intent(inout) :: Co
    !> Occupied MOs (hole - currently particle=hole always)
    type(tensor),intent(inout) :: Co2
    !> Virtual MOs (hole)
    type(tensor),intent(inout) :: Cv

    type(array2) :: fockt1_a2
    type(array2) :: Co_a2
    type(array2) :: Co2_a2
    type(array2) :: Cv_a2
    type(array2) :: t1_a2

    fockt1_a2%dims = fockt1%dims
    Co_a2%dims     = Co%dims
    Co2_a2%dims    = Co2%dims
    Cv_a2%dims     = Cv%dims
    t1_a2%dims     = t1%dims

    call ass_D1to2( fockt1%elm1,fockt1_a2%val,fockt1%dims )
    call ass_D1to2( Co%elm1,    Co_a2%val,    Co%dims     )
    call ass_D1to2( Co2%elm1,   Co2_a2%val,   Co2%dims    )
    call ass_D1to2( Cv%elm1,    Cv_a2%val,    Cv%dims     )
    call ass_D1to2( t1%elm1,    t1_a2%val,    t1%dims     )

    call Get_AOt1Fock(mylsitem,t1_a2,fockt1_a2,nocc,nvirt,nbasis,Co_a2,Co2_a2,Cv_a2)


    nullify( fockt1_a2%val )
    nullify( Co_a2%val     )
    nullify( Co2_a2%val    )
    nullify( Cv_a2%val     )
    nullify( t1_a2%val     )

    !call print_norm(fockt1)
  end subroutine Get_AOt1Fock_arraywrapper

  !> \brief Get T1 transformed Fock matrix in the AO basis using
  !> either fragment or full molecular T1 amplitudes.
  subroutine Get_AOt1Fock_oa(mylsitem,t1,fockt1,nocc,nvirt,nbasis,Co,Co2,Cv)
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
    type(array2),intent(inout) :: Co
    !> Occupied MOs (hole - currently particle=hole always)
    type(array2),intent(inout) :: Co2
    !> Virtual MOs (hole)
    type(array2),intent(inout) :: Cv
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
    call getT1transformation_occ(t1,xocc,yocc,Co,Co2,Cv)

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
    type(tensor), intent(inout) :: fock_array
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
    call util_get_symm_part(D)
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
    call util_get_symm_part(D)
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


  !> \brief: calculate atomic and pair fragment contributions to CCSD
  !> correlation energy for full molecule calculation.
  !> Currently, only for occupied partitioning scheme.
  !> \author: Janus Juul Eriksen
  !> \date: February 2013
  subroutine ccsd_energy_full_occ(nocc,nvirt,nfrags,offset,ccsd_doubles,ccsd_singles,integral,occ_orbitals,&
                           & eccsdpt_matrix_cou,eccsdpt_matrix_exc)

    implicit none

    !> ccsd doubles amplitudes and VOVO integrals (ordered as (a,b,i,j))
    type(tensor), intent(inout) :: ccsd_doubles, integral
    !> ccsd singles amplitudes
    type(tensor), intent(inout) :: ccsd_singles
    !> dimensions
    integer, intent(in) :: nocc, nvirt, nfrags, offset
    !> occupied orbital information
    type(decorbital), dimension(nocc+offset), intent(inout) :: occ_orbitals
    !> etot
    real(realk), dimension(nfrags,nfrags), intent(inout) :: eccsdpt_matrix_cou, eccsdpt_matrix_exc
    !> integers
    integer :: i,j,a,b,atomI,atomJ
    !> energy reals
    real(realk) :: energy_tmp_1, energy_tmp_2

    ! *************************************************************
    ! ************** do energy for full molecule ******************
    ! *************************************************************

    ! ***********************
    !   do CCSD energy part
    ! ***********************

    ! ***note: we only run over nval (which might be equal to nocc_tot if frozencore = .false.)
    ! so we only assign orbitals for the space in which the core orbitals (the offset) are omited

    !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp_1,energy_tmp_2),&
    !$OMP REDUCTION(+:eccsdpt_matrix_cou),&
    !$OMP SHARED(ccsd_doubles,ccsd_singles,integral,nocc,nvirt,occ_orbitals,offset,DECinfo)
    do j=1,nocc
    atomJ = occ_orbitals(j+offset)%CentralAtom
       do i=1,nocc
       atomI = occ_orbitals(i+offset)%CentralAtom

          do b=1,nvirt
             do a=1,nvirt

                energy_tmp_1 = ccsd_doubles%elm4(a,b,i,j) * integral%elm4(a,b,i,j)
                if(DECinfo%use_singles)then
                   energy_tmp_2 = ccsd_singles%elm2(a,i) * ccsd_singles%elm2(b,j) * integral%elm4(a,b,i,j)
                else
                   energy_tmp_2 = 0.0E0_realk
                endif
                eccsdpt_matrix_cou(AtomI,AtomJ) = eccsdpt_matrix_cou(AtomI,AtomJ) &
                                        & + energy_tmp_1 + energy_tmp_2

             end do
          end do

       end do
    end do
    !$OMP END PARALLEL DO

    ! reorder from (a,b,i,j) to (a,b,j,i)
    call tensor_reorder(integral,[1,2,4,3])

    !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp_1,energy_tmp_2),&
    !$OMP REDUCTION(+:eccsdpt_matrix_exc),&
    !$OMP SHARED(ccsd_doubles,ccsd_singles,integral,nocc,nvirt,occ_orbitals,offset,DECinfo)
    do j=1,nocc
    atomJ = occ_orbitals(j+offset)%CentralAtom
       do i=1,nocc
       atomI = occ_orbitals(i+offset)%CentralAtom

          do b=1,nvirt
             do a=1,nvirt

                energy_tmp_1 = ccsd_doubles%elm4(a,b,i,j) * integral%elm4(a,b,i,j)
                if(DECinfo%use_singles)then
                   energy_tmp_2 = ccsd_singles%elm2(a,i) * ccsd_singles%elm2(b,j) * integral%elm4(a,b,i,j)
                else
                   energy_tmp_2 = 0.0E0_realk
                endif
                eccsdpt_matrix_exc(AtomI,AtomJ) = eccsdpt_matrix_exc(AtomI,AtomJ) &
                                        & + energy_tmp_1 + energy_tmp_2

             end do
          end do

       end do
    end do
    !$OMP END PARALLEL DO

    ! get total fourth--order energy contribution
    eccsdpt_matrix_cou = 2.0E0_realk * eccsdpt_matrix_cou - eccsdpt_matrix_exc

    ! for the pair fragment energy matrix,
    ! we only consider pairs IJ where J>I; thus, move contributions and set J<I contribs to zero.
    ! (must be consistent with printout in print_pair_fragment_energies)

    do AtomI=1,nfrags
       do AtomJ=AtomI+1,nfrags

          eccsdpt_matrix_cou(AtomI,AtomJ) = eccsdpt_matrix_cou(AtomI,AtomJ) &
                                              & + eccsdpt_matrix_cou(AtomJ,AtomI)
          eccsdpt_matrix_cou(AtomJ,AtomI) = 0.0_realk
       end do
    end do


    ! ******************************************************************
    ! ************** done w/ energy for full molecule ******************
    ! ******************************************************************

  end subroutine ccsd_energy_full_occ


  !> \brief: print out CCSD fragment and pair interaction energies for full molecule calculation
  !          Only for occupied partitioning scheme.
  !          This routine should print the information in the same way as kasper's routine,
  !          print_all_fragment_energies in dec_utils.F90
  !
  !> \author: Janus Juul Eriksen, modified by Pablo Baudin to print (T) contributions.
  !> \date: February 2013
  subroutine print_fragment_energies_full(nfrags,FragEnergies,ccenergies,dofrag,distancetable)

    implicit none

    !> number of atoms in molecule
    integer, intent(in) :: nfrags
    !> matrices containing Frag. energies and interatomic distances
    real(realk), intent(in) :: FragEnergies(nfrags,nfrags,4), distancetable(nfrags,nfrags)
    !> Total cc energies:
    real(realk), intent(in) :: ccenergies(4)
    !> vector handling how the orbitals are assigned?
    logical, intent(inout) :: dofrag(nfrags)

    !> local variables 
    character(len=30) :: CorrEnergyString
    integer :: iCorrLen, cc_sol, pT_full, pT_4, pT_5
    logical :: print_pair, print_4pT_5pT

    print_pair = count(dofrag)>1
    print_4pT_5pT = DECinfo%PL>0
    CorrEnergyString = 'correlation energy            '
    iCorrLen = 18
    cc_sol  = 1
    pT_full = 2
    pT_4    = 3
    pT_5    = 4

    ! Print Header:
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '             ==============================================='
    write(DECinfo%output,*) '             |   Print single and pair fragment energies   |'
    write(DECinfo%output,*) '             ==============================================='
    write(DECinfo%output,*)

    if(.not.DECinfo%CCDhack)then
       if( DECinfo%ccmodel == MODEL_RPA)then
          call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),dofrag,&
             & 'RPA occupied single energies','AF_RPA_OCC')
          if (print_pair) call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),&
             & dofrag,Distancetable, 'RPA occupied pair energies','PF_RPA_OCC')

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,A,A,A,g20.10)') 'RPA ', &
             & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol)
          write(DECinfo%output,*)

       else if( DECinfo%ccmodel == MODEL_MP2)then
          call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),dofrag,&
             & 'MP2 occupied single energies','AF_MP2_OCC')
          if (print_pair) call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),&
             & dofrag,Distancetable, 'MP2 occupied pair energies','PF_MP2_OCC')

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,A,A,A,g20.10)') 'MP2 ', &
             & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol)
          write(DECinfo%output,*)

       else if( DECinfo%ccmodel == MODEL_CC2 )then
          call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),dofrag,&
             & 'CC2 occupied single energies','AF_CC2_OCC')
          if (print_pair) call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),&
             & dofrag,Distancetable, 'CC2 occupied pair energies','PF_CC2_OCC')

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,A,A,A,g20.10)') 'CC2 ', &
             & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol)
          write(DECinfo%output,*)

       else if( DECinfo%ccmodel == MODEL_CCSD )then 
          call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),dofrag,&
             & 'CCSD occupied single energies','AF_CCSD_OCC')
          if (print_pair) call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),&
             & dofrag,Distancetable, 'CCSD occupied pair energies','PF_CCSD_OCC')

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,A,A,A,g20.10)') 'CCSD ', &
             & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol)
          write(DECinfo%output,*)

       else if( DECinfo%ccmodel == MODEL_CCSDpT )then
          call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),dofrag,&
             & 'CCSD occupied single energies','AF_CCSD_OCC')
          if (print_pair) call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),&
             & dofrag,Distancetable, 'CCSD occupied pair energies','PF_CCSD_OCC')

          call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,pT_full),dofrag,&
             & '(T) occupied single energies','AF_ParT_OCC_BOTH')
          if (print_4pT_5pT) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,pT_4),dofrag,&
                & '(T) occupied single energies (fourth order)','AF_ParT_OCC4')
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,pT_5),dofrag,&
                & '(T) occupied single energies (fifth order)','AF_ParT_OCC5')
          end if
          
          if (print_pair) call print_pair_fragment_energies(nfrags,FragEnergies(:,:,pT_full),&
             & dofrag,Distancetable, '(T) occupied pair energies','PF_ParT_OCC_BOTH')
          if (print_4pT_5pT.and.print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,pT_4),&
                & dofrag,Distancetable, '(T) occupied pair energies (fourth order)','PF_ParT_OCC4')
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,pT_5),&
                & dofrag,Distancetable, '(T) occupied pair energies (fifth order)','PF_ParT_OCC5')
          end if

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCSD ', &
             & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) correlation energy  : ', &
             & ccenergies(pT_full)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) 4th order energy    : ', &
             & ccenergies(pT_4)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) 5th order energy    : ', &
             & ccenergies(pT_5)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'Total CCSD(T) ', &
             & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol)+ccenergies(pT_full)
          write(DECinfo%output,*)

       else
          call lsquit("ERROR(print_fragment_energies_full) model not implemented",-1)
       endif
    else
       call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),dofrag,&
          & 'CCD occupied single energies','AF_CCD_OCC')
       if (print_pair) call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol),&
          & dofrag,Distancetable, 'CCD occupied pair energies','PF_CCD_OCC')

       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,A,A,A,g20.10)') 'CCD ', &
          & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol)
       write(DECinfo%output,*)

    endif

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)'============================================================================='

  end subroutine print_fragment_energies_full

#ifdef MOD_UNRELEASED
  !============================================================================!
  !                   MO-based CCSD residual subroutines                       !
  !============================================================================!
  

  !> Purpose: Calculate CCSD residual using T1-transformed equations.
  !           This algorithm is MO-based and read MO-integral in memory,
  !           (PDM if MPI). The MO-integral are assumed to be packed
  !           in batches using the routine get_t1_free_gmo.
  !
  !> Author:  Pablo Baudin
  !> Date:    November 2013
  subroutine get_mo_ccsd_residual(ccmodel,pgmo_diag,pgmo_up,t1,omega1,t2,omega2, &
             & govov,nbas,nocc,nvir,iter,MOinfo,MyLsItem,lampo,lampv, &
             & lamho,lamhv,deltafock,fock,ppfock,pqfock,qpfock,qqfock,local)

    implicit none

    !> CC model
    integer,intent(in) :: ccmodel
    !> MO pack integrals; amplitudes and residuals:
    integer, intent(in) :: nbas, nocc, nvir, iter
    type(tensor), intent(inout) :: pgmo_diag, pgmo_up
    type(tensor), intent(inout) :: govov
    type(tensor), intent(inout) :: t1
    type(tensor), intent(inout) :: omega1
    type(tensor), intent(inout) :: t2
    type(tensor), intent(inout) :: omega2

    !> Long-range correction to Fock matrix
    type(tensor), intent(in) :: deltafock
    !> AO Fock matrix:
    type(tensor), intent(inout) :: fock
    !> occupied-occupied block of the t1-fock matrix
    type(tensor), intent(inout) :: ppfock
    !> virtual-virtual block of the t1-fock matrix
    type(tensor), intent(inout) :: qqfock
    !> occupied-virtual block of the t1-fock matrix
    type(tensor), intent(inout) :: pqfock
    !> virtual-occupied block of the t1-fock matrix
    type(tensor), intent(inout) :: qpfock
    !> transformation matrices from AO to t1-MO:
    real(realk), intent(in), pointer :: lampo(:,:), lampv(:,:)
    real(realk), intent(in), pointer :: lamho(:,:), lamhv(:,:)

    !> LS item with information needed for integrals
    type(lsitem), intent(inout) :: MyLsItem
     logical, intent(in) :: local

    !> Batches info:
    type(MObatchInfo), intent(in) :: MOinfo
    integer :: P_sta, P_end, Q_sta, Q_end, ntot
    integer :: dimMO, dimP, dimQ, Nbat, PQ_batch

    !> transformation matrices from MO to t1-MO:
    real(realk), pointer :: xvir(:), yocc(:)
    !> Unpack molecular integral:
    real(realk), pointer :: gmo(:)
    !> Intermediate for B2 term:
    real(realk), pointer :: B2prep(:)
    !> 2 coulomb - exchange amplitudes:
    real(realk), pointer :: u2(:,:,:,:) 
    !> Intermediates used to calculate the E2, A1 and B1 terms:
    real(realk), pointer :: G_Pi(:),  H_aQ(:)
    !> T1-transforled MO integrals:
    real(realk), pointer :: gvoov(:), gvvoo(:)
    real(realk), pointer :: goooo(:), govoo(:), gvooo(:)

    !> MPI info:
    integer, pointer :: joblist(:)
    logical :: master
    integer(kind=ls_mpik) :: tile_master, myrank, nnod
    integer :: ccmodel_copy
 
    !> Working arrays:
    real(realk), pointer :: tmp0(:), tmp1(:), tmp2(:) 
    integer(kind=long) :: tmp_size, no2v2

    integer :: i, a, O, V, N, X

    !> debug:
    logical :: print_debug, local_moccsd
    real(realk) :: tcpu, twall, tcpu1, twall1
    integer :: idb, iub

    call time_start_phase(PHASE_WORK)

    ! Initialize stuff
    nullify(xvir)
    nullify(yocc)
    nullify(gmo)
    nullify(B2prep)
    nullify(u2)
    nullify(G_Pi)
    nullify(H_aQ)
    nullify(tmp0)
    nullify(tmp1)
    nullify(tmp2)
    nullify(gvoov)
    nullify(gvvoo)
    nullify(goooo)
    nullify(govoo)
    nullify(gvooo)
    nullify(joblist)

    ! find out how the MO int. are stored:
    local_moccsd = .false.
    if (pgmo_diag%atype=='RTAR') local_moccsd = .true.

    ! Set MPI related info
    ! ********************
    myrank        = int(0,kind=ls_mpik)
    nnod          = 1
    master        = .true.
#ifdef VAR_MPI
    myrank        = infpar%lg_mynum
    nnod          = infpar%lg_nodtot
    master        = (myrank == infpar%master)
    call mem_alloc(joblist,MOinfo%Nbatch)
#endif

    ! shortcuts:
    ntot = nocc + nvir
    O = nocc
    V = nvir
    N = ntot
    X = MOinfo%DimInd1(1)
    print_debug = (DECinfo%PL>3.or.DECinfo%cc_driver_debug.and.master)

    ! Allocate working memory:
    dimMO = MOinfo%DimInd1(1)
    tmp_size = max(O**4, V*O**3, V*V*O*O, X*X*N*N, X*O*O*V, X*O*V*V)
    tmp_size = int(i8*tmp_size, kind=long)
    call mem_alloc(tmp0, tmp_size)

    tmp_size = max(X*X*N*N, O*O*V*N, O*O*X*N, O**4)
    tmp_size = int(i8*tmp_size, kind=long)
    call mem_alloc(tmp1, tmp_size)

    tmp_size = max(X*O*V*N, O*O*V*V, X*X*N*N, X*O*O*N)
    tmp_size = int(i8*tmp_size, kind=long)
    call mem_alloc(tmp2, tmp_size)

    call mem_alloc(xvir, int(i8*nvir*ntot, kind=long))
    call mem_alloc(yocc, int(i8*nocc*ntot, kind=long))
    call mem_alloc(gmo, int(i8*dimMO*dimMO*ntot*ntot, kind=long)) 
    call mem_alloc(B2prep, int(i8*nocc*nocc*nocc*nocc, kind=long))
    call mem_alloc(u2, nvir,nvir,nocc,nocc)
    call mem_alloc(G_Pi, ntot*nocc)
    call mem_alloc(H_aQ, nvir*ntot)

    tmp_size = int(i8*nocc*nocc*nocc*nocc, kind=long)
    call mem_alloc(goooo, tmp_size)
    tmp_size = int(i8*nvir*nocc*nocc*nocc, kind=long)
    call mem_alloc(govoo, tmp_size)
    call mem_alloc(gvooo, tmp_size)
    tmp_size = int(i8*nocc*nvir*nocc*nvir, kind=long)
    call mem_alloc(gvoov, tmp_size)
    call mem_alloc(gvvoo, tmp_size)

    ! start timings:
    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

#ifdef VAR_MPI
    ! Change array type to be dense:
    if (.not.local.and.master) then
      call tensor_cp_tiled2dense(t2,.false.)
      if(iter==1) then
        call memory_allocate_tensor_dense(govov)
      else
        call tensor_cp_tiled2dense(govov,.false.)
      end if
    end if 
#endif

    !===========================================================================
    ! Calculate transformation matrix:
    !  xvir_ap = delta_ap - t_ai delta_ip
    xvir = 0.0E0_realk
    call daxpy(nvir*nocc, -1.0E0_realk, t1%elm1, 1, xvir, 1)
    do a = 1, nvir
      xvir(nvir*nocc + a + (a-1)*nvir) = 1.0E0_realk
    end do

    !  yocc_ip = delta_ip + t'_ai delta_ap
    yocc = 0.0E0_realk
    call mat_transpose(nvir,nocc,1.0E0_realk,t1%elm1,0.0E0_realk,yocc(1+nocc*nocc:))
    do i = 1, nocc
      yocc(i + (i-1)*nocc) = 1.0E0_realk
    end do

    ! Initialization
    !$OMP WORKSHARE
    u2 = 0.0E0_realk
    B2prep = 0.0E0_realk
    G_Pi   = 0.0E0_realk
    H_aQ   = 0.0E0_realk
    gvoov  = 0.0E0_realk
    gvvoo  = 0.0E0_realk
    goooo  = 0.0E0_realk
    govoo  = 0.0E0_realk
    gvooo  = 0.0E0_realk
    !$OMP END WORKSHARE
    Nbat = MOinfo%nbatch

    !===========================================================================
    ! Calculate  2*coulomb - exchange doubles amplitudes:
    ! ordered as u2[ab, ij] = 2*t2[ab, ij] - t2[ab, ji]
    call daxpy(nvir*nvir*nocc*nocc,2.0E0_realk,t2%elm1,1,u2,1)
    call array_reorder_3d(-1.0E0_realk,t2%elm1,nvir*nvir,nocc,nocc,[1,3,2],1.0E0_realk,u2)

    call LSTIMER('MO-CCSD init calc.',tcpu1,twall1,DECinfo%output)


#ifdef VAR_MPI
    !===========================================================================
    !                          MPI COMMUNICATIONS
    !
    ! ~ Distribution of workloads to the nodes
    ! ~ Wake up slaves and communicate data
    ! ~ Reduce MO batches on appropriate nodes if local_moccsd
    call get_mo_ccsd_joblist(MOinfo, joblist, pgmo_diag, pgmo_up)

    ! Wake up slaves and communicate important data
    call time_start_phase(PHASE_COMM)
    StartUpSlaves: if (master.and.nnod>1) then
      ccmodel_copy = ccmodel
      call ls_mpibcast(MOCCSDDATA,infpar%master,infpar%lg_comm)
      call mpi_communicate_moccsd_data(ccmodel_copy,pgmo_diag,pgmo_up,t1,t2,omega2, &
             & govov,nbas,nocc,nvir,iter,MOinfo,MyLsItem,local)
    end if StartUpSlaves
    call time_start_phase(PHASE_WORK)

    call LSTIMER('MO-CCSD MPI-comm.',tcpu1,twall1,DECinfo%output)

    ! If RTAR type of array is used then we reduce the partial MO batches
    ! to the MPI process that will treat it (depending on joblist)
    ! (For TDAR type of array the batches are already complete due to one-sided
    ! communication in get_t1_free_gmo)
    if (iter==1.and.local_moccsd) then 
      do PQ_batch=1,Nbat 
        tile_master = joblist(PQ_batch) - 1
        
        P_sta = MOinfo%StartInd1(PQ_batch)
        Q_sta = MOinfo%StartInd2(PQ_batch)

        if (P_sta==Q_sta) then
          idb = MOinfo%tileInd(PQ_batch,1)
          call time_start_phase(PHASE_COMM)
          call lsmpi_reduction(pgmo_diag%ti(idb)%t,pgmo_diag%ti(idb)%e,tile_master,infpar%lg_comm)
          call time_start_phase(PHASE_WORK)
        else
          iub = MOinfo%tileInd(PQ_batch,1)
          call time_start_phase(PHASE_COMM)
          call lsmpi_reduction(pgmo_up%ti(iub)%t,pgmo_up%ti(iub)%e,tile_master,infpar%lg_comm)
          call time_start_phase(PHASE_WORK)
        end if
      end do
    end if

   govov%access_type  = AT_ALL_ACCESS
   t2%access_type     = AT_ALL_ACCESS
   omega2%access_type = AT_ALL_ACCESS
   if (.not.local) then
     call memory_allocate_tensor_dense_pc(omega2)
     omega2%itype = TT_DENSE
     govov%itype  = TT_DENSE
   end if
#endif
    if (iter==1) call tensor_zero(govov)
    call tensor_zero(omega2)

    call LSTIMER('MO-CCSD INIT',tcpu1,twall1,DECinfo%output)



    !===========================================================================!
    !                        START LOOP OVER MO BATCHES                         !
    !===========================================================================!

    BatchPQ: do PQ_batch = 1, Nbat

      if (nnod>1) then 
        ! check if the current job is to be done by current node
        if (joblist(PQ_batch)/=(myrank+1)) cycle
      end if

      P_sta = MOinfo%StartInd1(PQ_batch)
      dimP  = MOinfo%dimInd1(PQ_batch)
      Q_sta = MOinfo%StartInd2(PQ_batch)
      dimQ  = MOinfo%dimInd2(PQ_batch)

      ! Get batch of MO integral
      if (P_sta==Q_sta) then
        idb = MOinfo%tileInd(PQ_batch,1)
        call unpack_gmo(gmo,pgmo_diag,idb,ntot,dimP,dimQ,.true.,tmp0)
      else
        iub = MOinfo%tileInd(PQ_batch,1)
        call unpack_gmo(gmo,pgmo_up,iub,ntot,dimP,dimQ,.false.,tmp0)
      end if

      ! Get intermediate for the calculation of residual
      call wrapper_get_intermediates(ccmodel,ntot,nocc,nvir,dimP,dimQ,P_sta,Q_sta,iter,gmo, &
                         & xvir,yocc,t2%elm1,u2,goooo,B2prep,omega2%elm1,G_Pi,H_aQ, &
                         & govov%elm1,gvoov,gvvoo,govoo,gvooo,tmp0,tmp1,tmp2)

    end do BatchPQ

    call LSTIMER('MO-CCSD main loop',tcpu1,twall1,DECinfo%output)



    !===========================================================================
    !           REDUCE DATA AND GET FINAL CONTRIBUTION TO RESIDUAL
    !
    ! Calculate norm of A2:
    if (print_debug.and.nnod==1) call print_norm(omega2,'debug: residual A2 norm:           ')

    if (ccmodel>MODEL_CC2) then
      ! get final B2 term and add it to residual:
      call dgemm('n','n',nvir*nvir,nocc*nocc,nocc*nocc,0.5E0_realk,t2%elm1,nvir*nvir, &
                & B2prep,nocc*nocc,0.5E0_realk,omega2%elm1,nvir*nvir)
       
      ! Calculate norm of A2 + B2 residual:
      if (print_debug.and.nnod==1) call print_norm(omega2,'debug: residual B2 norm:           ')
    end if

    call mpi_reduction_after_main_loop(ccmodel,ntot,nvir,nocc,iter,G_Pi,H_aQ,goooo, &
                & govoo,gvooo,gvoov,gvvoo,govov)

    call mem_dealloc(tmp1)
    call mem_dealloc(tmp2)

    if (ccmodel>MODEL_CC2) then
      call LSTIMER('MO-CCSD A2 B2 + comm',tcpu1,twall1,DECinfo%output)

      ! Get C2 and D2 terms
      call wrapper_get_C2_and_D2(tmp0,tmp1,tmp2,t2,u2,govov,gvoov,gvvoo,nocc,nvir,omega2)
    end if

    nullify(tmp1)
    nullify(tmp2)

#ifdef VAR_MPI
    call time_start_phase(PHASE_COMM)
    no2v2 = int(nvir*nvir*nocc*nocc, kind=long)
    call lsmpi_local_reduction(omega2%elm1,no2v2,infpar%master)
    call time_start_phase(PHASE_WORK)
#endif

    ! Slaves exit the routines:
    if(.not. master) then
      govov%access_type  = AT_MASTER_ACCESS
      t2%access_type     = AT_MASTER_ACCESS
      omega2%access_type = AT_MASTER_ACCESS
      ! Free arrays:
      call mem_dealloc(xvir)
      call mem_dealloc(yocc)
      call mem_dealloc(gmo)
      call mem_dealloc(B2prep)
      call mem_dealloc(u2)
      call mem_dealloc(G_Pi)
      call mem_dealloc(H_aQ)
      call mem_dealloc(tmp0)
       
      call mem_dealloc(gvoov)
      call mem_dealloc(gvvoo)
      call mem_dealloc(goooo)
      call mem_dealloc(govoo)
      call mem_dealloc(gvooo)

      call mem_dealloc(joblist)
 
      return
    endif

    ! Calculate norm of A2 + B2 + C2 + D2 residual:
    if (print_debug.and.ccmodel>MODEL_CC2) then 
      call print_norm(omega2,'debug: residual D2 norm:            ')
      call LSTIMER('MO-CCSD A2, B2, C2, D2',tcpu1,twall1,DECinfo%output)
    end if


    !===========================================================================
    !                          GET MO-FOCK MATRICES
    call get_MO_fock_matrices(ccmodel,nbas,nocc,nvir,lampo,lampv,lamho,lamhv,tmp0, &
                  & goooo,govoo,gvooo,gvoov,gvvoo,ppfock%elm1,pqfock%elm1, & 
                  & qpfock%elm1,qqfock%elm1,fock%elm1,deltafock%elm1,MyLsItem)

    if (print_debug) then
      call print_norm(ppfock,"MO-CCSD (ppfock):                ")
      call print_norm(pqfock,"MO-CCSD (pqfock):                ")
      call print_norm(qpfock,"MO-CCSD (qpfock):                ")
      call print_norm(qqfock,"MO-CCSD (qqfock):                ")
      call LSTIMER('MO-CCSD Fock mat',tcpu1,twall1,DECinfo%output)
    end if

 
    !===========================================================================
    !                       GET SINGLES CCSD RESIDUAL
    !
    !CCD can be achieved by not using singles residual updates here
    if(.not. DECinfo%CCDhack)then
      ! Get A1 term:
      ! Omega_ai = xvir_aP * G_Pi 
      call dgemm('n','n',nvir,nocc,ntot,1.0E0_realk,xvir,nvir,G_Pi,ntot, &
                & 0.0E0_realk,omega1%elm1,nvir)
      ! Calculate norm of A1:
      if (print_debug) call print_norm(omega1,'debug: residual A1 norm:         ')
       
      ! Calculate B1 term:
      ! Omega_ai += - H_aQ * yocc_iQ 
      call dgemm('n','t',nvir,nocc,ntot,-1.0E0_realk,H_aQ,nvir,yocc,nocc, &
                & 1.0E0_realk,omega1%elm1,nvir)
      ! Calculate norm of A1 + B1:
      if (print_debug) call print_norm(omega1,'debug: residual B1 norm:         ')
       
      ! Calculate C1 term:
      ! Reorder u2[ac, ik] -> u2[ai, kc]
      call array_reorder_4d(1.0E0_realk,u2,nvir,nvir,nocc,nocc,[1,3,4,2], & 
                & 0.0E0_realk,tmp0)
       
      ! Omega_ai += u2[ai, kc] * F_kc 
      call dgemv('n',nvir*nocc,nvir*nocc,1.0E0_realk,tmp0,nvir*nocc,pqfock%elm1,1, &
                & 1.0E0_realk,omega1%elm1,1)
      ! Calculate norm of A1 + B1 + C1:
      if (print_debug) call print_norm(omega1,'debug: residual C1 norm:         ')
       
      ! Calculate D1 term:
      ! Omega_ai += F_ai
      call daxpy(nocc*nvir,1.0E0_realk,qpfock%elm1,1,omega1%elm1,1)
       
      ! Calculate norm of full single residual:
      if (print_debug) then
        call print_norm(omega1,'debug: residual D1 norm:                 ')
        call LSTIMER('MO-CCSD singles',tcpu1,twall1,DECinfo%output)
      end if
    end if

    !===========================================================================
    !                   GET DOUBLES CCSD RESIDUAL AND FINALIZE
    !
    ! Get E2 term and introduce permutational symmetry
    call get_E2_and_permute(ccmodel,ntot,nocc,nvir,ppfock%elm1,qqfock%elm1,tmp0,&
         & t2%elm1,G_Pi,H_aQ,omega2%elm1)

    ! Calculate norm of full double residual:
    if (print_debug) then
      call print_norm(omega2,'debug: residual E2 norm:            ')
      call LSTIMER('MO-CCSD E2',tcpu1,twall1,DECinfo%output)
    end if

    ! Free arrays:
    call mem_dealloc(xvir)
    call mem_dealloc(yocc)
    call mem_dealloc(gmo)
    call mem_dealloc(B2prep)
    call mem_dealloc(u2)
    call mem_dealloc(G_Pi)
    call mem_dealloc(H_aQ)
    call mem_dealloc(tmp0)

    call mem_dealloc(gvoov)
    call mem_dealloc(gvvoo)
    call mem_dealloc(goooo)
    call mem_dealloc(govoo)
    call mem_dealloc(gvooo)

#ifdef VAR_MPI
    call mem_dealloc(joblist)

    ! Move dense part to tiles
    if(.not.local)then
      govov%access_type  = AT_MASTER_ACCESS
      t2%access_type     = AT_MASTER_ACCESS
      omega2%access_type = AT_MASTER_ACCESS
      if (iter==1) then
        call tensor_mv_dense2tiled(govov,.true.)
      else
        call memory_deallocate_tensor_dense_pc(govov)
        govov%itype      = TT_TILED_DIST
      end if
      call tensor_mv_dense2tiled(omega2,.true.)
      call tensor_mv_dense2tiled(t2,.true.)
    endif
#endif 

    call LSTIMER('MO-CCSD residual',tcpu,twall,DECinfo%output)

  end subroutine get_mo_ccsd_residual


  !> Purpose: Call subroutine to get intermediates and contributions to ccsd
  !           residual
  !
  !> Author:  Pablo Baudin
  !> Date:    November 2013
  subroutine wrapper_get_intermediates(ccmodel,ntot,nocc,nvir,dimP,dimQ,P_sta,Q_sta,iter,gmo, &
                                     & xvir,yocc,t2,u2,goooo,B2prep,omega2,G_Pi,H_aQ, &
                                     & govov,gvoov,gvvoo,govoo,gvooo,tmp0,tmp1,tmp2)

    implicit none

    !> CC model
    integer,intent(in) :: ccmodel
    !> dimensions for arrays:
    integer, intent(in) :: ntot, nocc, nvir
    integer, intent(in) :: dimP, dimQ, P_sta, Q_sta
    !> CC iterations counter:
    integer, intent(in) :: iter
    !> MO integral (pc|rd):
    real(realk), intent(in) :: gmo(dimP*dimQ*ntot*ntot)
    !> transformation matrices:
    real(realk), intent(in) :: xvir(nvir*ntot), yocc(nocc*ntot)
    !> doubles amplitudes:
    real(realk), intent(in) :: t2(nvir,nvir,nocc,nocc)
    !> 2*Coulomb - Exchange form of amplitudes:
    real(realk), intent(in) :: u2(nvir,nvir,nocc,nocc)
    !> Intermediates used to calculate the B2 term.
    !  B2prep = g_kilj + sigma_kilj
    real(realk), intent(inout) :: B2prep(nocc*nocc*nocc*nocc)
    !> doubles residual array:
    real(realk), intent(inout) :: omega2(nvir,nvir,nocc,nocc)
    !> Intermediates used to calculate the E2, A1 and B1 terms:
    real(realk), intent(inout) :: G_Pi(ntot*nocc),  H_aQ(nvir*ntot)
    !> T1-transformed MO integrals:
    real(realk), intent(inout) :: govov(:), gvoov(:), gvvoo(:)
    real(realk), intent(inout) :: goooo(:), govoo(:), gvooo(:)
    !> working arrays:
    real(realk), intent(inout) :: tmp0(:), tmp1(:), tmp2(:)


    call get_A2_and_B2prep_terms(ccmodel,ntot,nocc,nvir,dimP,dimQ,P_sta,Q_sta,gmo, &
                        & xvir,yocc,t2,goooo,B2prep,omega2,tmp0,tmp1,tmp2)
    call get_G_and_H_intermeditates(ntot,nocc,nvir,dimP,dimQ, &
                        &  P_sta,Q_sta,gmo,u2,G_Pi,H_aQ,tmp0,tmp1,tmp2)
    call get_MO_integrals(ntot,nocc,nvir,dimP,dimQ,P_sta,Q_sta,iter,gmo, &
                        & xvir,yocc,govov,gvoov,gvvoo,govoo,gvooo,tmp0,tmp1)

    ! If the PQ batch is an upper diagonal block, we repeat the oprerations
    ! with the transposed batch:
    if (P_sta<Q_sta) then
      ! transpose gmo to get batch QP:
      call array_reorder_3d(1.0E0_realk,gmo,dimP,dimQ,ntot*ntot, &
                          & [2,1,3],0.0E0_realk,tmp1)
      call dcopy(dimQ*dimP*ntot*ntot,tmp1,1,gmo,1)

      call get_A2_and_B2prep_terms(ccmodel,ntot,nocc,nvir,dimQ,dimP,Q_sta,P_sta,gmo, &
                          & xvir,yocc,t2,goooo,B2prep,omega2,tmp0,tmp1,tmp2)
      call get_G_and_H_intermeditates(ntot,nocc,nvir,dimQ,dimP, &
                          & Q_sta,P_sta,gmo,u2,G_Pi,H_aQ,tmp0,tmp1,tmp2)
      call get_MO_integrals(ntot,nocc,nvir,dimQ,dimP,Q_sta,P_sta,iter,gmo, &
                          & xvir,yocc,govov,gvoov,gvvoo,govoo,gvooo,tmp0,tmp1)
    end if

  end subroutine wrapper_get_intermediates  


  !> Purpose: Calculate the A2 part of the doubles amplitudes residual
  !           and some intermediates for the B2 part.
  !
  !> Author:  Pablo Baudin
  !> Date:    November 2013
  subroutine get_A2_and_B2prep_terms(ccmodel,ntot,nocc,nvir,dimP,dimQ,P_sta,Q_sta,gmo, &
                          & xvir,yocc,t2,goooo,B2prep,omega2,tmp0,tmp1,tmp2)

    implicit none

    !> CC model
    integer,intent(in) :: ccmodel
    !> dimensions for arrays:
    integer, intent(in) :: ntot, nocc, nvir
    integer, intent(in) :: dimP, dimQ, P_sta, Q_sta
    !> MO integral (pc|rd):
    real(realk), intent(in) :: gmo(dimP*dimQ*ntot*ntot)
    !> transformation matrices:
    real(realk), intent(in) :: xvir(nvir*ntot), yocc(nocc*ntot)
    !> doubles amplitudes:
    real(realk), intent(in) :: t2(nvir,nvir,nocc,nocc)
    !> full occ. T1 transformed integral:
    real(realk), intent(inout) :: goooo(nocc*nocc*nocc*nocc)
    !> Intermediates used to calculate the B2 term.
    !  B2prep = g_kilj + sigma_kilj
    real(realk), intent(inout) :: B2prep(nocc*nocc*nocc*nocc)
    !> doubles residual array:
    real(realk), intent(inout) :: omega2(nvir,nvir,nocc,nocc)
    !> working arrays:
    real(realk), intent(inout) :: tmp0(:), tmp1(:), tmp2(:)

    !> orbital indices for loops:
    integer :: i, j, k, ij, n_ij, r, d, dimI, I_sta, dimK, K_sta, dimC, C_sta
    !> variable indices for arrays:
    integer :: pos1, pos2, ncopy
    integer :: mv((nvir*nvir)/2), st
    real(realk), external :: ddot

    n_ij = nocc*(nocc+1)/2

    ! get set of virt orb in Q batch
    ! and set of occ orb in P_batch
    dimI = nocc - Q_sta + 1
    if (dimI<0) dimI = 0
    if (dimI>dimQ) dimI = dimQ
    dimC = dimQ - dimI

    if (Q_sta<=nocc) then
      C_sta = 1
    else
      C_sta = Q_sta - nocc
    end if

    I_sta = nocc - dimI + 1
    dimK = dimI
    if (P_sta/=Q_sta) then
      dimK = nocc - P_sta + 1
      if (dimK<0) dimK = 0
      if (dimK>dimP) dimK = dimP
    end if

    !===========================================================================
    ! Get g[aibj] A2.1 term and add to residual:

    ! transform s => j: g[PQrj]
    call dgemm('n','T',dimP*dimQ*ntot,nocc,ntot,1.0E0_realk,gmo,dimP*dimQ*ntot, &
              & yocc,nocc,0.0E0_realk,tmp1,dimP*dimQ*ntot)
  
    ! transform P => a: g[Qrja]
    pos1 = 1 + (P_sta-1)*nvir
    call dgemm('T','T',dimQ*ntot*nocc,nvir,dimP,1.0E0_realk,tmp1,dimP, &
              & xvir(pos1),nvir,0.0E0_realk,tmp2,dimQ*ntot*nocc)

    ! transform Q => i: g[rjai]
    pos1 = 1 + (Q_sta-1)*nocc
    call dgemm('T','T',ntot*nocc*nvir,nocc,dimQ,1.0E0_realk,tmp2,dimQ, &
              & yocc(pos1),nocc,0.0E0_realk,tmp1,ntot*nocc*nvir)

    ! transform r => b:  g[jaib]
    call dgemm('T','T',nocc*nvir*nocc,nvir,ntot,1.0E0_realk,tmp1,ntot, &
              & xvir,nvir,0.0E0_realk,tmp2,nocc*nvir*nocc)

    ! reorder and add to omega: g[abij]
    call array_reorder_4d(1.0E0_realk,tmp2,nocc,nvir,nocc,nvir,[2,4,3,1], &
              & 1.0E0_realk,omega2)
    !===========================================================================
     
    if (ccmodel==MODEL_CC2) return

    !===========================================================================
    ! add O4 int. to B2 term
    if (dimK>0) then
      ! reorder g[PQrs] to g[QsPr]:
      call array_reorder_4d(1.0E0_realk,gmo,dimP,dimQ,ntot,ntot,[2,4,1,3],0.0E0_realk,tmp0)
 
      ! Get g[QsKl from g[QsPr]:
      ncopy = dimQ*ntot*dimK
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ncopy,nocc,ntot,dimP,dimQ,tmp0,tmp2)&
      !$OMP PRIVATE(i,pos1,pos2)
      do i =1,nocc
        pos1 = 1 + (i-1)*dimQ*ntot*dimP
        pos2 = 1 + ncopy*(i-1)
        call dcopy(ncopy,tmp0(pos1),1,tmp2(pos2),1)
      end do
      !$OMP END PARALLEL DO
  
      ! transpose g[Q,sKl] to g[sKl, Q]
      call mat_transpose(dimQ, ntot*dimK*nocc, 1.0E0_realk, tmp2, 0.0E0_realk, tmp1)

      ! transform Q => i and get: g[sK li]
      pos1 = 1 + (Q_sta-1)*nocc
      call dgemm('n','t',ntot*dimK*nocc,nocc,dimQ,1.0E0_realk,tmp1,ntot*dimK*nocc, &
                & yocc(pos1),nocc,0.0E0_realk,tmp2,ntot*dimK*nocc)

      ! transform s => j and get: g[Kl ij]
      call dgemm('t','t',dimK*nocc*nocc,nocc,ntot,1.0E0_realk,tmp2,ntot, &
                & yocc,nocc,0.0E0_realk,tmp1,dimK*nocc*nocc)
     
      ! add to previous loops:
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(nocc,dimK,P_sta,tmp1,B2prep)&
      !$OMP PRIVATE(i,pos1,pos2)
      do i=1,nocc*nocc*nocc
        pos1 = 1 + (i-1)*dimK
        pos2 = P_sta + (i-1)*nocc
        call daxpy(dimK,1.0E0_realk,tmp1(pos1),1,B2prep(pos2),1)
      end do
      !$OMP END PARALLEL DO
      call array_reorder_4d(1.0E0_realk,tmp1,dimK,nocc,nocc,nocc, &
                            & [1,3,2,4],0.0E0_realk,tmp2)
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(nocc,dimK,P_sta,tmp2,goooo)&
      !$OMP PRIVATE(i,pos1,pos2)
      do i=1,nocc*nocc*nocc
        pos1 = 1 + (i-1)*dimK
        pos2 = P_sta + (i-1)*nocc
        call daxpy(dimK,1.0E0_realk,tmp2(pos1),1,goooo(pos2),1)
      end do
      !$OMP END PARALLEL DO
    end if
    !===========================================================================


    !===========================================================================
    if (dimC<=0) return
    ! else calculate sigma and A2.2 and B2.2 contribution:

    ! Extract t2[Cd, i<=j] from t2[cd, ij]:
    call dcopy(nvir*nvir*nocc*nocc,t2,1,tmp2,1)
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(nocc,nvir,C_sta,dimC,tmp1,tmp2)&
    !$OMP PRIVATE(i,j,d,pos1,pos2)
    do j = 1,nocc
      do i = 1,j
        do d = 1, nvir
          pos2 = C_sta + (d-1)*nvir + (i-1)*nvir*nvir + (j-1)*nvir*nvir*nocc
          pos1 = 1 + (d-1)*dimC + (j*(j-1)/2 + i-1)*dimC*nvir
          call dcopy(dimC,tmp2(pos2),1,tmp1(pos1),1)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! Calculate sigma[p r; i<=j]
    ! get g[PrQs]:
    if (dimK>0) then
      call mat_transpose(dimQ*ntot, dimP*ntot, 1.0E0_realk, tmp0, 0.0E0_realk, tmp2)
    else
      call array_reorder_4d(1.0E0_realk,gmo,dimP,dimQ,ntot,ntot,[1,3,2,4],0.0E0_realk,tmp2)
    end if
     
    ! get g[PrCd] from g[PrQs] 
    ncopy = dimP*ntot*dimC
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ncopy,ntot,nocc,nvir,dimI,dimC,dimP,dimQ,tmp0,tmp2)&
    !$OMP PRIVATE(i,j,pos1,pos2)
    do i =1,nvir
      pos1 = 1 + dimI*dimP*ntot + (nocc+i-1)*dimP*ntot*dimQ
      pos2 = 1 + ncopy*(i-1)
      call dcopy(ncopy,tmp2(pos1),1,tmp0(pos2),1)
    end do
    !$OMP END PARALLEL DO
     
    ! Get: sigma[pr, ij] = sum_cd g[pr, cd] * t2red[cd, ij]
    !n_ij= nocc*nocc
    call dgemm('n','n',dimP*ntot,n_ij,dimC*nvir,1.0E0_realk,tmp0,dimP*ntot, &
              & tmp1,dimC*nvir,0.0E0_realk,tmp2,dimP*ntot)
    !===========================================================================


    !===========================================================================
    ! add sigma to B2 term:
    !
    ! prep B2 term: sigma[K l, i j] <= sigma[P r i<=j]
    ! occupied indices do not need to be transformed.
    if (dimK>0) then
      !$OMP WORKSHARE 
      tmp1 = 0.0E0_realk
      !$OMP END WORKSHARE
      ! get occ indices from full: sigma[K l, i<=j] <= sigma[P r, i<=j]
      do i=1,n_ij
        do r=1,nocc
          pos1 = 1 + (r-1)*dimP + (i-1)*dimP*ntot
          pos2 = P_sta + (r-1)*nocc + (i-1)*nocc*nocc
          call dcopy(dimK,tmp2(pos1),1,tmp1(pos2),1) 
        end do
      end do
      ! Copy contracted array into proper place in the full array
      ! first spread i<=j indices to full ij indices
      !$OMP PARALLEL DEFAULT(NONE) SHARED(nocc,tmp1)&
      !$OMP PRIVATE(i,j,k,pos1,pos2)
      do j=nocc,2,-1
        do i=j,1,-1
          pos1=((i+j*(j-1)/2)-1)*nocc*nocc
          pos2=(i-1)*nocc*nocc+(j-1)*nocc*nocc*nocc
          !$OMP DO
          do k=1, nocc*nocc
            tmp1(k+pos2) = tmp1(k+pos1)
          end do
          !$OMP END DO
        enddo
      enddo
      !$OMP BARRIER
      !$OMP DO
      ! then square the array, copy full index ij to full index ji:
      do j=nocc,1,-1
        do i=j,1,-1
          pos1=1+(i-1)*nocc*nocc+(j-1)*nocc*nocc*nocc
          pos2=1+(j-1)*nocc*nocc+(i-1)*nocc*nocc*nocc
          if(i/=j) tmp1(pos2:pos2+nocc*nocc-1) = tmp1(pos1:pos1+nocc*nocc-1)
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      do j=1,nocc
        do i=j+1,nocc
          pos1=1+(i-1)*nocc*nocc+(j-1)*nocc*nocc*nocc
          call alg513(tmp1(pos1:nocc*nocc+pos1-1),nocc,nocc,nocc*nocc,mv,(nocc*nocc)/2,st)
        enddo
      enddo
      call daxpy(nocc*nocc*nocc*nocc,1.0E0_realk,tmp1,1,B2prep,1)
       
    end if
    !===========================================================================


    ! Transpose sigma matrix from sigma[p; r i<=j ] to sigma[r i<=j; p]
    call mat_transpose(dimP, ntot*n_ij, 1.0E0_realk, tmp2, 0.0E0_realk, tmp1)

    ! Get A2.2 term:
    ! Transform r -> b
    call dgemm('n','n',nvir,n_ij*dimP,ntot,1.0E0_realk,xvir,nvir, &
              & tmp1, ntot, 0.0E0_realk, tmp2, nvir)

    ! Transform p -> a; order is now: sigma[a b i j]
    call dgemm('n','t',nvir,nvir*n_ij,dimP,1.0E0_realk,xvir(1+(P_sta-1)*nvir),nvir, &
              & tmp2, nvir*n_ij, 0.0E0_realk, tmp1, nvir)


    ! Sum up sigma PQ batches contributions to A2.2 part of CCSD residual:
    !$OMP PARALLEL DEFAULT(NONE) SHARED(nocc,nvir,tmp1)&
    !$OMP PRIVATE(i,j,k,pos1,pos2)
    do j=nocc,2,-1
      do i=j,1,-1
        pos1=((i+j*(j-1)/2)-1)*nvir*nvir
        pos2=(i-1)*nvir*nvir+(j-1)*nocc*nvir*nvir
        !$OMP DO
        do k=1, nvir*nvir
          tmp1(k+pos2) = tmp1(k+pos1)
        end do
        !$OMP END DO
      enddo
    enddo
    !$OMP BARRIER
    !$OMP DO
    do j=nocc,1,-1
      do i=j,1,-1
        pos1=1+(i-1)*nvir*nvir+(j-1)*nocc*nvir*nvir
        pos2=1+(j-1)*nvir*nvir+(i-1)*nocc*nvir*nvir
        if(i/=j) tmp1(pos2:pos2+nvir*nvir-1) = tmp1(pos1:pos1+nvir*nvir-1)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    do j=1,nocc
      do i=j+1,nocc
        pos1=1+(i-1)*nvir*nvir+(j-1)*nocc*nvir*nvir
        call alg513(tmp1(pos1:nvir*nvir+pos1-1),nvir,nvir,nvir*nvir,mv,(nvir*nvir)/2,st)
      enddo
    enddo
    call daxpy(nocc*nocc*nvir*nvir,1.0E0_realk,tmp1,1,omega2,1)

  end subroutine get_A2_and_B2prep_terms


  !> Purpose: Calculate intermediates G and H used for the construction of
  !           the E2, A1, and B1 term of the CCSD residual.
  !
  !              G_Pi = sum_ckd u2[cd, ki] g[pd, kc]
  !              H_aQ = sum_ckl u2[ac, kl] g[kq, lc]
  !
  !> Author:  Pablo Baudin
  !> Date:    November 2013
  subroutine get_G_and_H_intermeditates(ntot,nocc,nvir,dimP,dimQ,P_sta,Q_sta, &
                            & gmo,u2,G_Pi,H_aQ,tmp0,tmp1,tmp2)

    implicit none

    !> dimensions for arrays:
    integer, intent(in) :: ntot, nocc, nvir
    integer, intent(in) :: dimP, dimQ, P_sta, Q_sta
    !> MO integral (pc|rd):
    real(realk), intent(in) :: gmo(dimP*dimQ*ntot*ntot)
    !> 2*Coulomb - Exchange form of amplitudes:
    real(realk), intent(in) :: u2(nvir,nvir,nocc,nocc)
    !> Intermediates used to calculate the E2, A1 and B1 terms:
    real(realk), intent(inout) :: G_Pi(ntot*nocc),  H_aQ(nvir*ntot)
    !> working arrays:
    real(realk), intent(inout) :: tmp0(:), tmp1(:), tmp2(:)

    !> orbital indices for loops:
    integer :: i, k, d, c, l, q, dimI, I_sta, dimK, K_sta, dimD, D_sta
    !> variable indices for arrays:
    integer :: pos1, pos2, ncopy

    ! get set of occ/virt orb in Q batch (i/d)
    ! and set of occ orb in P batch (k) 
    dimI = nocc - Q_sta + 1
    if (dimI<0) dimI = 0
    if (dimI>dimQ) dimI = dimQ
    dimD = dimQ - dimI

    if (Q_sta<=nocc) then
      D_sta = 1
    else
      D_sta = Q_sta - nocc
    end if

    dimK = dimI
    if (P_sta/=Q_sta) then
      dimK = nocc - P_sta + 1
      if (dimK<0) dimK = 0
      if (dimK>dimP) dimK = dimP
    end if
    K_sta = nocc - dimK + 1
    
    !========================================================
    ! Get G_Pi
    if (dimD>0) then

      ! Get u2[cD, ki] from u2[cd, ki]; D in Q batch
      ncopy = nvir*dimD*nocc*nocc
      pos1 = D_sta+dimD-1
      call dcopy(ncopy,u2(:,D_sta:pos1,:,:),1,tmp0,1)
       
      ! Get g[PD, kc] from g[PQ, rs]; D in Q batch
      ncopy = dimP*dimD
      pos2 = 1
      do c=1,nvir
        do k=1,nocc
          pos1 = 1 + dimP*dimI + dimP*dimQ*ntot*(nocc+c-1) + (k-1)*dimP*dimQ
          call dcopy(ncopy,gmo(pos1),1,tmp2(pos2),1)
          pos2 = pos2 + ncopy
        end do
      end do
       
      call array_reorder_4d(1.0E0_realk,tmp2,dimP,dimD,nocc,nvir, &
                            & [1,4,2,3],0.0E0_realk,tmp1)

      ! Get G_Pi:
      call dgemm('n','n',dimP,nocc,nvir*dimD*nocc,1.0E0_realk,tmp1,dimP, &
                & tmp0,nvir*dimD*nocc,1.0E0_realk,G_Pi(P_sta),ntot)
    end if

    !========================================================
    ! Get H_aQ
    if (dimK>0) then 

      ! Get u2[ac Kl] from u2[ac, kl]; K in P batch
      ncopy = nvir*nvir*dimK*nocc
      pos1 = P_sta+dimK-1
      call dcopy(ncopy,u2(:,:,P_sta:pos1,:),1,tmp0,1)

      ! get g[cKl, Q] from g[PQrs]
      do c=1,nvir
        do l=1,nocc
          do q=1,dimQ
            pos1 = 1 + (q-1)*dimP + (l-1)*dimP*dimQ + (nocc+c-1)*dimP*dimQ*ntot
            pos2 = c + (l-1)*nvir*dimK + (q-1)*nvir*dimK*nocc
            call dcopy(dimK,gmo(pos1),1,tmp1(pos2),nvir)
          end do 
        end do 
      end do 
      
      ! Get H_aQ:
      call dgemm('n','n',nvir,dimQ,nvir*dimK*nocc,1.0E0_realk,tmp0,nvir, &
                & tmp1,nvir*dimK*nocc,1.0E0_realk,H_aQ(1+(Q_sta-1)*nvir),nvir)
    end if

  end subroutine get_G_and_H_intermeditates

 
  !> Purpose: Get MO integrals needed for the C2, B2 and E2 terms of the
  !           ccsd residual
  !
  !> Author: Pablo Baudin
  !> Date:   November 2013
  subroutine get_MO_integrals(ntot,nocc,nvir,dimP,dimQ,P_sta,Q_sta,iter,gmo, &
                             & xvir,yocc,govov,gvoov,gvvoo,govoo,gvooo,tmp0,tmp1)
    implicit none

    !> dimensions for arrays:
    integer, intent(in) :: ntot, nocc, nvir
    integer, intent(in) :: dimP, dimQ, P_sta, Q_sta
    !> CC iterations counter:
    integer, intent(in) :: iter
    !> MO integral (pc|rd):
    real(realk), intent(in) :: gmo(dimP*dimQ*ntot*ntot)
    !> transformation matrices:
    real(realk), intent(in) :: xvir(nvir*ntot), yocc(nocc*ntot)
    !> T1-transformed MO integrals:
    real(realk), intent(inout) :: govov(:), gvoov(:), gvvoo(:)
    real(realk), intent(inout) :: govoo(:), gvooo(:)
    !> working arrays:
    real(realk), intent(inout) :: tmp0(:), tmp1(:)

    !> orbital indices for loops:
    integer :: i, a, j, b, s, dimI, I_sta, dimK, K_sta
    integer :: dimA, A_sta, dimB, B_sta
    !> variable indices for arrays:
    integer :: pos1, pos2, ncopy

    ! get set of occ/virt orb in Q batch (k/a)
    ! and set of occ orb in P batch (i) 
    dimK = nocc - Q_sta + 1
    if (dimK<0) dimK = 0
    if (dimK>dimQ) dimK = dimQ
    dimA = dimQ - dimK

    if (Q_sta<=nocc) then
      A_sta = 1
    else
      A_sta = Q_sta - nocc
    end if

    dimI = dimK
    if (P_sta/=Q_sta) then
      dimI = nocc - P_sta + 1
      if (dimI<0) dimI = 0
      if (dimI>dimP) dimI = dimP
    end if
    I_sta = nocc - dimI + 1

    !====================================================================
    ! Get govov integrals:  g[iajb]
    if ((iter==1).and.(dimI>0)) then 
      do b=1,nvir
        do j=1,nocc
          do a=1,dimA
            pos1 = 1 + dimP*dimK + (a-1)*dimP + (j-1)*dimP*dimQ &
                 & + (nocc+b-1)*dimP*dimQ*ntot

            pos2 = P_sta + (A_sta+a-2)*nocc + (j-1)*nocc*nvir + (b-1)*nocc*nvir*nocc

            call dcopy(dimI,gmo(pos1),1,govov(pos2),1)
          end do
        end do
      end do
    end if
 
    !====================================================================
    ! Get gvooo integrals: g[aijk]
    ! 1) get g[PQjs]
    ncopy = dimP*dimQ*nocc
    pos2 = 1 
    do b=1,ntot
      pos1 = 1 + (b-1)*dimP*dimQ*ntot
      call dcopy(ncopy,gmo(pos1),1,tmp0(pos2),1)
      pos2 = pos2 + ncopy
    end do
    
    ! 2) transform s to k g[PQjs]
    call dgemm('n','t',dimP*dimQ*nocc,nocc,ntot,1.0E0_realk,tmp0,dimP*dimQ*nocc, &
              & yocc,nocc,0.0E0_realk,tmp1,dimP*dimQ*nocc)

    ! 3) transpose and get g[Qjk,P]
    call mat_transpose(dimP,dimQ*nocc*nocc,1.0E0_realk,tmp1,0.0E0_realk,tmp0)

    ! 4) transform Q to i:
    pos1 = 1 + (Q_sta-1)*nocc
    call dgemm('n','n',nocc,nocc*nocc*dimP,dimQ,1.0E0_realk,yocc(pos1),nocc, &
              & tmp0,dimQ,0.0E0_realk,tmp1,nocc)

    ! 4) transform P to a:
    pos1 = 1 + (P_sta-1)*nvir
    call dgemm('n','t',nvir,nocc*nocc*nocc,dimP,1.0E0_realk,xvir(pos1),nvir, &
              & tmp1,nocc*nocc*nocc,1.0E0_realk,gvooo,nvir)


    !====================================================================
    ! Get gvoov integrals: g[aijb]
    ! 1) get g[PQjb]
    ncopy = dimP*dimQ*nocc
    pos2 = 1 
    do b=1,nvir
      pos1 = 1 + (nocc+b-1)*dimP*dimQ*ntot
      call dcopy(ncopy,gmo(pos1),1,tmp0(pos2),1)
      pos2 = pos2 + ncopy
    end do
    
    ! 2) transpose and get g[QjbP]:
    call mat_transpose(dimP,dimQ*nocc*nvir,1.0E0_realk,tmp0,0.0E0_realk,tmp1)
  
    ! 3) transform Q to i:
    pos1 = 1 + (Q_sta-1)*nocc
    call dgemm('n','n',nocc,nocc*nvir*dimP,dimQ,1.0E0_realk,yocc(pos1),nocc, &
              & tmp1,dimQ,0.0E0_realk,tmp0,nocc)

    ! 4) transform P to a:
    pos1 = 1 + (P_sta-1)*nvir
    call dgemm('n','t',nvir,nocc*nocc*nvir,dimP,1.0E0_realk,xvir(pos1),nvir, &
              & tmp0,nocc*nocc*nvir,1.0E0_realk,gvoov,nvir)


    !====================================================================
    dimB = dimA
    !dimI = dimK
    B_sta = A_sta
    ! Get gvvoo integrals: g[abij]
    if (dimB>0) then
      ! 1) get g[PBis]
      ncopy = dimP*dimB
      pos2 = 1 
      do s=1,ntot
        do i=1,nocc
          pos1 = 1 + dimP*dimK + (i-1)*dimP*dimQ + (s-1)*dimP*dimQ*ntot
          call dcopy(ncopy,gmo(pos1),1,tmp0(pos2),1)
          pos2 = pos2 + ncopy
        end do
      end do
     
      ! 2) transform s to j:
      call dgemm('n','t',dimP*dimB*nocc,nocc,ntot,1.0E0_realk,tmp0,dimP*dimB*nocc, &
                & yocc,nocc,0.0E0_realk,tmp1,dimP*dimB*nocc)

      ! 3) get g[KBij] from g[PBij]
      if (dimI>0) then
        do i=1,nocc*nocc
          do s=1,dimB
            pos1 = 1 + (s-1)*dimP + (i-1)*dimP*dimB
            pos2 = P_sta + (B_sta+s-2)*nocc + (i-1)*nocc*nvir
            call daxpy(dimI,1.0E0_realk,tmp1(pos1),1,govoo(pos2),1)
          end do
        end do
      end if

      ! 4) transform P to a:
      pos1 = 1 + (P_sta-1)*nvir
      call dgemm('n','n',nvir,dimB*nocc*nocc,dimP,1.0E0_realk,xvir(pos1),nvir, &
                & tmp1,dimP,0.0E0_realk,tmp0,nvir)

      ! 5) copy tmp0 to gvvoo array:
      ncopy = nvir*dimB
      do i=1,nocc*nocc
        pos1 = 1 + (i-1)*ncopy
        pos2 = 1 + (B_sta-1)*nvir + (i-1)*nvir*nvir
        call daxpy(ncopy,1.0E0_realk,tmp0(pos1),1,gvvoo(pos2),1)
      end do

    end if

  end subroutine get_MO_integrals
  

  !> Purpose: MPI reduction for intermediate arrays and integrals after the 
  !           main loop.
  !
  !> Author:  Pablo Baudin
  !> Date:    April 2014
  subroutine  mpi_reduction_after_main_loop(ccmodel,nt,nv,no,iter,G_Pi,H_aQ,goooo, &
                & govoo,gvooo,gvoov,gvvoo,govov)

    implicit none
  
    integer, intent(in) :: ccmodel, nt, nv, no, iter
    real(realk), intent(inout) :: G_Pi(:), H_aQ(:)
    real(realk), intent(inout) :: goooo(:), govoo(:), gvooo(:)
    real(realk), intent(inout) :: gvoov(:), gvvoo(:)
    type(tensor), intent(inout) :: govov
    
    integer(kind=long) :: no2v2

#ifdef VAR_MPI
    call time_start_phase(PHASE_COMM)
    no2v2 = int(nv*nv*no*no, kind=long)
    call lsmpi_local_reduction(goooo,no**4,infpar%master)
    call lsmpi_local_reduction(govoo,nv*no**3,infpar%master)
    call lsmpi_local_reduction(gvooo,nv*no**3,infpar%master)
    call lsmpi_local_reduction(G_Pi,nt*no,infpar%master)
    call lsmpi_local_reduction(H_aQ,nv*nt,infpar%master)

    if (ccmodel>MODEL_CC2) then
      ! ALL REDUCE FOR C2 AND D2 TERMS WITH MPI:
      call lsmpi_allreduce(gvoov,no2v2,infpar%lg_comm)
      call lsmpi_allreduce(gvvoo,no2v2,infpar%lg_comm)
      if (iter==1) call lsmpi_allreduce(govov%elm1,no2v2,infpar%lg_comm)
    else if(ccmodel==MODEL_CC2) then
      call lsmpi_local_reduction(gvoov,no2v2,infpar%master)
      call lsmpi_local_reduction(gvvoo,no2v2,infpar%master)
      if (iter==1) call lsmpi_local_reduction(govov%elm1,no2v2,infpar%master)
    end if
    call time_start_phase(PHASE_WORK)
#endif

  end subroutine  mpi_reduction_after_main_loop


  subroutine wrapper_get_C2_and_D2(tmp0,tmp1,tmp2,t2,u2,govov,gvoov,gvvoo,no,nv,omega2)

    implicit none
    
    real(realk), intent(inout) :: tmp0(:)
    real(realk), pointer :: tmp1(:), tmp2(:)
    real(realk), intent(in) :: u2(:,:,:,:)
    type(tensor) :: t2
    type(tensor), intent(inout) :: omega2
    type(tensor), intent(inout) :: govov
    real(realk), intent(inout) :: gvoov(:), gvvoo(:)
    integer, intent(in) :: no, nv

    type(tensor) :: u2a, govova, gvvooa, gvoova
    real(realk) :: time_CND_work,time_CND_comm

    time_CND_work = 0.0E0_realk
    time_CND_comm = 0.0E0_realk

    call tensor_init(u2a, [nv,nv,no,no],4)
    call array_reorder_3d(1.0E0_realk,u2,nv,nv,no*no,[2,1,3],0.0E0_realk,u2a%elm1)

    ! gvoov [aijb] must be ordered as [ajbi]:
    call tensor_init(gvoova,[nv,no,nv,no],4)
    call array_reorder_4d(1.0E0_realk,gvoov,nv,no,no,nv,[1,3,4,2],0.0E0_realk,gvoova%elm1)

    ! gvvoo [abij] must be ordered as [aijb]:
    call tensor_init(gvvooa,[nv,no,no,nv],4)
    call array_reorder_4d(1.0E0_realk,gvvoo,nv,nv,no,no,[1,3,4,2],0.0E0_realk,gvvooa%elm1)

    call get_cnd_terms_mo_3n4(tmp0,tmp1,tmp2,t2,u2a,govov,gvoova,gvvooa, &
       & no,nv,omega2,4,.false.,0_long,time_CND_work,time_CND_comm)

    call tensor_free(u2a)
    call tensor_free(gvoova)
    call tensor_free(gvvooa)

  end subroutine wrapper_get_C2_and_D2


  !> Purpose: Calculate fock matrices in MO basis using T1-transformed MO
  !           integrals and T1 transformation matrices.
  !
  !> Author:  Pablo Baudin
  !> Date:    Movember 2013
  subroutine get_MO_fock_matrices(ccmodel,nbas,nocc,nvir,lampo,lampv,lamho,lamhv,tmp0, &
             & goooo,govoo,gvooo,gvoov,gvvoo,Foo,Fov,Fvo,Fvv,fock,deltafock,MyLsItem)
 
    implicit none

    !> CC model:
    integer, intent(in) :: ccmodel
    integer, intent(in) :: nbas, nocc, nvir
    !> Transformation matrices:
    real(realk), intent(in) :: lampo(nbas,nocc), lampv(nbas,nvir)
    real(realk), intent(in) :: lamho(nbas,nocc), lamhv(nbas,nvir)
    !> Working array:
    real(realk) :: tmp0(:)
    !> T1-transformed MO integrals:
    real(realk), intent(in) :: goooo(:), govoo(:), gvooo(:)
    real(realk), intent(in) :: gvoov(:), gvvoo(:)
    !> AO fock matrix:
    real(realk), intent(in) :: fock(nbas*nbas)
    !> T1-transformed MO inactive Fock matrices:
    real(realk), intent(inout) :: Foo(nocc*nocc), Fov(nocc*nvir)
    real(realk), intent(inout) :: Fvo(nvir*nocc), Fvv(nvir*nvir)
    !> Long-range correction to Fock matrix
    real(realk), intent(in) :: deltafock(nbas,nbas)
    !> LS item with information needed for integrals
    type(lsitem), intent(inout) :: MyLsItem
  
    !> Inactive AO Fock matrix:
    type(matrix) :: iFock
    integer :: pos1, i
    real(realk), external :: ddot

    ! allocate and get 1-electron AO fock matrix:
    call mat_init(iFock,nbas,nbas)
    call mat_zero(iFock)
    call ii_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsItem%setting,&
         & iFock%elms,nbas,nbas,AORdefault,AORdefault)

    ! KK: Add long-range Fock correction
    call daxpy(nbas*nbas,1.0E0_realk,deltafock,1,iFock%elms,1)


    !Transform 1-electron inactive Fock matrix into the different MO subspaces
    ! -> Fop
    call dgemm('t','n',nocc,nbas,nbas,1.0E0_realk,lampo,nbas,iFock%elms,nbas, &
              & 0.0E0_realk,tmp0,nocc)

    if (ccmodel>MODEL_CC2) then
      ! -> Foo
      call dgemm('n','n',nocc,nocc,nbas,1.0E0_realk,tmp0,nocc,lamho,nbas, &
                & 0.0E0_realk,Foo,nocc)
    end if

    ! -> Fov
    call dgemm('n','n',nocc,nvir,nbas,1.0E0_realk,tmp0,nocc,lamhv,nbas, &
              & 0.0E0_realk,Fov,nocc)
    ! -> Fvp
    call dgemm('t','n',nvir,nbas,nbas,1.0E0_realk,lampv,nbas,iFock%elms,nbas, &
              & 0.0E0_realk,tmp0,nvir)
    ! -> Fvo
    call dgemm('n','n',nvir,nocc,nbas,1.0E0_realk,tmp0,nvir,lamho,nbas, &
              & 0.0E0_realk,Fvo,nvir)

    if (ccmodel>MODEL_CC2) then
      ! -> Fvv
      call dgemm('n','n',nvir,nvir,nbas,1.0E0_realk,tmp0,nvir,lamhv,nbas, &
                & 0.0E0_realk,Fvv,nvir)
    end if

    ! Free the 1-electron AO fock matrix
    call mat_free(iFock)

    !===================================================================
    ! Get two-electron contribution to MO Fock matrix:

    if (ccmodel>MODEL_CC2) then
      ! Foo:
      call array_reorder_3d(1.0E0_realk,goooo,nocc,nocc*nocc,nocc,[1,3,2], &
                & 0.0E0_realk,tmp0)
      do i=1, nocc
        pos1 = 1 + (i-1)*nocc*nocc*(nocc+1)
        call daxpy(nocc*nocc,2.0E0_realk,goooo(pos1),1,Foo,1)
        call daxpy(nocc*nocc,-1.0E0_realk,tmp0(pos1),1,Foo,1)
      end do
    end if

    ! Fov:
    call array_reorder_4d(1.0E0_realk,govoo,nocc,nvir,nocc,nocc,[3,2,4,1], &
              & 0.0E0_realk,tmp0)
    do i=1, nocc
      pos1 = 1 + (i-1)*nocc*nvir*(nocc+1)
      call daxpy(nocc*nvir,2.0E0_realk,govoo(pos1),1,Fov,1)
      call daxpy(nocc*nvir,-1.0E0_realk,tmp0(pos1),1,Fov,1)
    end do

    ! Fvo:
    call array_reorder_3d(1.0E0_realk,gvooo,nvir,nocc*nocc,nocc,[1,3,2], &
              & 0.0E0_realk,tmp0)
    do i=1, nocc
      pos1 = 1 + (i-1)*nvir*nocc*(nocc+1)
      call daxpy(nvir*nocc,2.0E0_realk,gvooo(pos1),1,Fvo,1)
      call daxpy(nvir*nocc,-1.0E0_realk,tmp0(pos1),1,Fvo,1)
    end do

    if (ccmodel>MODEL_CC2) then
      ! Get Fvv:
      call array_reorder_3d(1.0E0_realk,gvoov,nvir,nocc*nocc,nvir,[1,3,2], &
                & 0.0E0_realk,tmp0)
      do i=1, nocc
        pos1 = 1 + (i-1)*nvir*nvir*(nocc+1)
        call daxpy(nvir*nvir,2.0E0_realk,gvvoo(pos1),1,Fvv,1)
        call daxpy(nvir*nvir,-1.0E0_realk,tmp0(pos1),1,Fvv,1)
      end do
    else if (ccmodel==MODEL_CC2) then
      ! get Block diag. MO-Fock matrices for CC2 model:
      ! -> Foo
      call dgemm('t','n',nocc,nbas,nbas,1.0E0_realk,lampo,nbas,fock,nbas, &
                & 0.0E0_realk,tmp0,nocc)
      call dgemm('n','n',nocc,nocc,nbas,1.0E0_realk,tmp0,nocc,lamho,nbas, &
                & 0.0E0_realk,Foo,nocc)
      ! -> Fvv
      call dgemm('t','n',nvir,nbas,nbas,1.0E0_realk,lampv,nbas,fock,nbas, &
                & 0.0E0_realk,tmp0,nvir)
      call dgemm('n','n',nvir,nvir,nbas,1.0E0_realk,tmp0,nvir,lamhv,nbas, &
                & 0.0E0_realk,Fvv,nvir)
    end if
 

  end subroutine get_MO_Fock_matrices


  !> Purpose: Calculate E2 term of the CCSD residual using G and H
  !           intermediates and the fock matrix.
  !           E2 is added to get the full double residual on which we
  !           apply a permutational operator.
  !
  !> Author:  Pablo Baudin
  !> Date:    Novemeber 2013
  subroutine get_E2_and_permute(ccmodel,ntot,nocc,nvir,ppfock,qqfock,tmp0,t2,G_Pi,H_aQ,omega2)

    implicit none

    integer, intent(in) :: ccmodel, ntot, nocc, nvir
    real(realk), intent(in) :: ppfock(:), qqfock(:)
    real(realk), intent(inout) :: tmp0(:)
    real(realk), intent(in) :: t2(nvir,nvir,nocc,nocc)
    real(realk), intent(inout) :: G_Pi(:), H_aQ(:)
    real(realk), intent(inout) :: omega2(nvir,nvir,nocc,nocc)

    integer :: i
    real(realk), pointer :: tmp1(:)
    real(realk), external :: ddot

    call mem_alloc(tmp1, int(i8*nocc*nocc*nvir*nvir, kind=long))

    !===================================================================
    ! Calculate contribution from E2.1 to residual:
    ! transpose t2: -> t2[ij, ac]
    call mat_transpose(nvir*nvir, nocc*nocc, 1.0E0_realk, t2, 0.0E0_realk, tmp0)

    ! for CCSD H_aQ := - F_bc + H_bc
    if (ccmodel>MODEL_CC2) then
      call daxpy(nvir*nvir,-1.0E0_realk,qqfock,1,H_aQ(1+nvir*nocc),1)
       
      ! Get E2.1[ijab] = - t2[ij, ac] * tmp_bc
      call dgemm('n','t',nocc*nocc*nvir,nvir,nvir,-1.0E0_realk,tmp0,nocc*nocc*nvir, &
                & H_aQ(1+nvir*nocc),nvir,0.0E0_realk,tmp1,nocc*nocc*nvir)

      ! Transpose and add to omega2:
      call mat_transpose(nocc*nocc, nvir*nvir, 1.0E0_realk, tmp1, 1.0E0_realk, omega2)
    else if (ccmodel==MODEL_CC2) then
      ! for CC2:
      ! Get E2.1[ijab] = t2[ij, ac] * F__bc
      call dgemm('n','t',nocc*nocc*nvir,nvir,nvir,1.0E0_realk,tmp0,nocc*nocc*nvir, &
                & qqfock,nvir,0.0E0_realk,tmp1,nocc*nocc*nvir)

      ! Transpose and add to omega2:
      call mat_transpose(nocc*nocc, nvir*nvir, 1.0E0_realk, tmp1, 0.5E0_realk, omega2)
    end if


    !===================================================================
    ! Calculate contribution from E2.2 to residual:
    ! Get G_kj from G_pi:
    if (ccmodel>MODEL_CC2) then
      do i=1,nocc
        call dcopy(nocc,G_Pi(1+(i-1)*ntot),1,tmp0(1+(i-1)*nocc),1)
      end do

      ! Sum F_kj and G_kj:
      call daxpy(nocc*nocc,1.0E0_realk,ppfock,1,tmp0,1)

      ! Omega2 += - t[abik] * tmp0_kj
      call dgemm('n','n',nvir*nvir*nocc,nocc,nocc,-1.0E0_realk,t2,nvir*nvir*nocc, &
                & tmp0,nocc,1.0E0_realk,omega2,nvir*nvir*nocc)
    else if (ccmodel==MODEL_CC2) then
      ! Omega2 += - t[abik] * F_kj
      call dgemm('n','n',nvir*nvir*nocc,nocc,nocc,-1.0E0_realk,t2,nvir*nvir*nocc, &
                & ppfock,nocc,1.0E0_realk,omega2,nvir*nvir*nocc)
    end if

    !===================================================================
    ! Introduce permutational symmetry:
    call dcopy(nvir*nvir*nocc*nocc,omega2,1,tmp0,1)
      
    call array_reorder_4d(1.0E0_realk,tmp0,nvir,nvir,nocc,nocc,[2,1,4,3], &
                         & 1.0E0_realk,omega2)

    call mem_dealloc(tmp1)

  end subroutine get_E2_and_permute
#endif

  subroutine ccsd_debug_print(ccmodel,print_nr,master,local,&
        &scheme,print_debug,o2v2,w1,omega2,govov,gvvooa,gvoova)
     implicit none
     integer,intent(in) :: ccmodel
     integer,intent(in) :: print_nr
     integer,intent(in) :: scheme
     integer(kind=8)    :: o2v2
     logical,intent(in) :: print_debug,master,local
     type(mpi_realk) :: w1
     type(tensor) :: omega2,govov,gvvooa,gvoova
     character(tensor_MSG_LEN) :: msg

#ifdef VAR_LSDEBUG
     select case(print_nr)
     case(1)
        if(print_debug)then
           !DEBUG PRINT NORM OMEGA
#ifdef VAR_MPI
           if(.not.local)then
              if(alloc_in_dummy)then
                 if(omega2%lock_set(1))call lsmpi_win_flush(omega2%wi(1))
              else
                 call tensor_unlock_wins(omega2,check=.true.)
              endif
           endif
#endif
           call print_norm(omega2," NORM(omega2):",print_on_rank = 0)

           !DEBUG PRINT NORM GOVOV
#ifdef VAR_MPI
           if(.not.local)then
              if(alloc_in_dummy)then
                 if(govov%lock_set(1))call lsmpi_win_flush(govov%wi(1))
              else
                 call tensor_unlock_wins(govov,check=.true.)
              endif
           endif
#endif
            call print_norm(govov," NORM(govov):",print_on_rank = 0)

           if (ccmodel>MODEL_CC2) then
              !DEBUG PRINT NORM GVVOO
#ifdef VAR_MPI
              if(.not.local.and..not.scheme==4)then
                 if(alloc_in_dummy)then
                    if(gvvooa%lock_set(1))call lsmpi_win_flush(gvvooa%wi(1))
                 else
                    call tensor_unlock_wins(gvvooa,check=.true.)
                 endif
              endif
#endif
              call print_norm(gvvooa," NORM(gvvoo):",print_on_rank = 0)
               
              !DEBUG PRINT NORM GVOOV
#ifdef VAR_MPI
              if(.not.local.and..not.scheme==4)then
                 if(alloc_in_dummy)then
                    if(gvoova%lock_set(1))call lsmpi_win_flush(gvoova%wi(1))
                 else
                    call tensor_unlock_wins(gvoova,check=.true.)
                 endif
              endif
#endif
              call print_norm(gvoova," NORM(gvoov):",print_on_rank = 0)
           endif
        endif
     case(2)
        if(print_debug.and.ccmodel>MODEL_CC2)then
#ifdef VAR_MPI
           if(.not.local)then
              if(alloc_in_dummy)then
                 if(omega2%lock_set(1))call lsmpi_win_flush(omega2%wi(1))
              else
                 call tensor_unlock_wins(omega2,check=.true.)
              endif
           endif
#endif
           call print_norm(omega2," NORM(omega2 after B2.2):",print_on_rank = 0)
        endif
     case(3)
        if(print_debug.and.ccmodel>MODEL_CC2)then
#ifdef VAR_MPI
           if(.not.local)then
              if(alloc_in_dummy)then
                 if(omega2%lock_set(1))call lsmpi_win_flush(omega2%wi(1))
              else
                 call tensor_unlock_wins(omega2,check=.true.)
              endif
           endif
#endif
           call print_norm(omega2," NORM(omega2 after CND):",print_on_rank=0)
        endif
     case default
        print *,"WARNING(ccsd_debug_print):unknown debug print selected"
     end select
#endif
  end subroutine ccsd_debug_print


  function get_wsize_for_ccsd_int_direct(wnr,no,os,nv,vs,nb,bs,&
        &nba,nbg,nbuffs,s,setting,intspec) result(wsize)
     implicit none
     integer, intent(in) :: wnr,no,os,nv,vs,nb,bs,nba,nbg,nbuffs,s
     type(lssetting),intent(inout) :: setting
     Character,intent(in) :: intspec(5)
     integer(kind=long) :: wsize
     integer(kind=long) :: maxsize64,nor,nvr, MAX_INTEGRAL_BUF
     integer :: nbuffs,me

     nor = (i8*(no*(no+1))/2)
     nvr = (i8*(nv*(nv+1))/2)

     maxsize64 = 0
     select case(wnr)
     case(0)
        maxsize64 = int((i8*nb*nb)*nba*nbg,kind=8)
        if(s==2) maxsize64 = max(maxsize64,int((2*vs*vs*os)*os,kind=8))
        if(s==0)then

           maxsize64 = max(max(nbg*nv,nba*no),max(nbg*no,nba*nv))
           maxsize64 = max(maxsize64,nbuffs*max(max(i8*bs,i8*vs),i8*os)**4)

        endif
     case(1)
        maxsize64 = max(int((i8*nb*nb)*nba*nbg,kind=8),int((i8*nv*nv*no)*nba,kind=8))
        maxsize64 = max(maxsize64,int((i8*no*no*nv)*nbg,kind=8))
        if(s==4.or.s==3) maxsize64 = max(maxsize64,int((i8*no*no*nv)*nba,kind=8))
        if(s==2) maxsize64 = max(maxsize64,int((2*vs*vs*os)*os,kind=8))
        if(s==0)then
           maxsize64 = nbuffs * int(max(max(os,vs),bs),kind=8)**4
           maxsize64 = max(maxsize64,nbuffs*int((i8*vs**2)*os**2,kind=8))
           MAX_INTEGRAL_BUF = 0
           call simulate_intloop_and_get_worksize(MAX_INTEGRAL_BUF,nb,nbg,nba,bs,intspec,setting)
           maxsize64 = max(maxsize64,MAX_INTEGRAL_BUF)
        endif
     case(2)
        maxsize64 = max(int((i8*nb)*nb*nba*nbg,kind=8),(i8*no*no)*nv*nv)
        maxsize64 = max(maxsize64,int(nor*no*no,kind=8))
        if(s==2)then
           maxsize64 = max(maxsize64,int((no*no*no)*no,kind=8))
           maxsize64 = max(maxsize64,int((no*no*no)*nbg,kind=8))
        endif
        if(s==0)then
           maxsize64 = 0
        endif
     case(3)
        maxsize64 = max(int((i8*nv)*no*nba*nbg,kind=8),int((i8*no*no)*nba*nbg,kind=8))
        maxsize64 = max(maxsize64,int((i8*no*no*nv)*nba,kind=8))
        maxsize64 = max(maxsize64,int((2_long*nor)*nba*nbg,kind=8)) 
        maxsize64 = max(maxsize64,int((i8*nor)*nv*nba,kind=8)) 
        maxsize64 = max(maxsize64,int((i8*nor)*nv*nbg,kind=8)) 
        maxsize64 = max(maxsize64,int((i8*no)*nor*nba,kind=8)) 
        maxsize64 = max(maxsize64,int((i8*no)*nor*nbg,kind=8)) 
        if(s==2)then
            maxsize64 = max(maxsize64,int((2_long*vs*vs*os)*os,kind=8))
            maxsize64 = max(maxsize64,int((i8*no*no)*no*no,kind=8))
        endif
        if(s==0)then
           maxsize64 = 0
        endif
     case default
        call lsquit("ERROR(get_wsize_for_ccsd_int_direct):unknown identifier",-1)
     end select

     wsize = maxsize64

     !Sanity checks for matrix sizes which need to be filled
     if(wsize>MAXINT)then
        me = 0
#ifdef VAR_MPI
        me = infpar%lg_mynum
#endif
        print *,me,"ERROR(get_wsize_for_ccsd_int_direct)",wnr,no,os,nv,vs,nb,nba,nbg,s,wsize,MAXINT
        call lsquit("ERROR(CCSD):matrix sizes too large, please recompile with 64bit integers",-1)
     endif

  end function get_wsize_for_ccsd_int_direct

end module ccsd_module


#ifdef VAR_MPI
!> \brief slave function for data preparation
!> \author Patrick Ettenhuber
!> \date March 2012
subroutine ccsd_data_preparation()
  use, intrinsic :: iso_c_binding, only:c_ptr
  use precision
  use lstiming
  use lsparameters, only: CCSDDATA
  use dec_typedef_module
  use typedeftype,only:lsitem,tensor
  use infpar_module
  use lsmpi_type, only:ls_mpibcast,ls_mpibcast,LSMPIBROADCAST,MPI_COMM_NULL,&
  &ls_mpiInitBuffer,ls_mpi_buffer,ls_mpiFinalizeBuffer
  use lsmpi_op, only:mpicopy_lsitem
  use daltoninfo, only:ls_free
  use memory_handling, only: mem_alloc, mem_dealloc
  use tensor_interface_module, only: tensor_ainit,tensor_free,&
      &memory_allocate_tensor_dense,memory_deallocate_tensor_dense,&
      &memory_deallocate_window,&
      &lspdm_use_comm_proc,get_tensor_from_parr,TT_DENSE,AT_ALL_ACCESS,AT_MASTER_ACCESS
  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
  use decmpi_module, only: mpi_communicate_ccsd_calcdata
  use array2_simple_operations,only: array2_free,array2_init
  use array4_simple_operations,only: array4_free,array4_init
  use ccsd_module, only:get_ccsd_residual_integral_driven,yet_another_ccsd_residual,&
     &RN_RESIDUAL_INT_DRIVEN,RN_YET_ANOTHER_RES
  implicit none
  !type(mp2_batch_construction) :: bat
  type(array2) :: deltafock,ppfock,qqfock,pqfock,qpfock,omega1,fock
  type(array2) :: xocc,xvirt,yocc,yvirt
  type(tensor)  :: om2,t2,govov
  type(lsitem) :: MyLsItem
  logical :: local
  integer :: nbas,nocc,nvirt,scheme,ccmodel,res_nr
  integer(kind=long) :: nelms
  integer      :: iter,k,n4,i
  real(realk),pointer  :: xodata(:),xvdata(:),yodata(:),yvdata(:),&
                         & df(:),f(:),ppf(:),qqf(:),pqf(:),qpf(:),om1(:),t2data(:),&
                         & t2d(:,:,:,:),xod(:,:),xvd(:,:),yod(:,:),yvd(:,:)
  type(c_ptr)           :: xoc,xvc,yoc,yvc,&
                         & dfc,fc,ppfc,qqfc,pqfc,qpfc,om1c,t2c,&
                         & t2dc,xodc,xvdc,yodc,yvdc
  integer(kind=ls_mpik) :: xow,xvw,yow,yvw,&
                         & dfw,fw,ppfw,qqfw,pqfw,qpfw,om1w,t2w,&
                         & t2dw,xodw,xvdw,yodw,yvdw
  logical :: parent
  integer :: addr01(infpar%pc_nodtot)
  integer :: addr02(infpar%pc_nodtot)
  integer :: addr03(infpar%pc_nodtot)
  integer(kind=ls_mpik) :: lg_me,lg_nnod
  integer(kind=8) :: ne


  lg_me   = infpar%lg_mynum
  lg_nnod = infpar%lg_nodtot
  parent  = (infpar%parent_comm == MPI_COMM_NULL)

  call time_start_phase(PHASE_COMM)
  !note that for the slave all allocatable arguments are just dummy indices
  !the allocation and broadcasting happens in here
  call mpi_communicate_ccsd_calcdata(ccmodel,om2,t2,govov,xodata,xvdata,yodata,yvdata,&
     &MyLsItem,nbas,nvirt,nocc,iter,local,res_nr)

  if(local)then
      call tensor_ainit(t2,    [nvirt,nvirt,nocc,nocc], 4, local=local, atype='LDAR' )
      call tensor_ainit(govov, [nocc,nvirt,nocc,nvirt], 4, local=local, atype='LDAR' )
      call tensor_ainit(om2,   [nvirt,nvirt,nocc,nocc], 4, local=local, atype='LDAR' )
  endif
  
  ! Quantities, that need to be defined and setset
  ! ********************************************************
   

  !split messages in 2GB parts, compare to counterpart in
  !mpi_communicate_ccsd_calcdate

  nelms = nbas*nocc
  call mem_alloc( yodata, nelms )
  call mem_alloc( xodata, nelms )
  call ls_mpibcast( xodata, nelms, infpar%master, infpar%lg_comm)
  call ls_mpibcast( yodata, nelms, infpar%master, infpar%lg_comm)

  nelms = nbas*nvirt
  call mem_alloc( xvdata, nelms )
  call mem_alloc( yvdata, nelms )
  call ls_mpibcast( xvdata, nelms, infpar%master, infpar%lg_comm )
  call ls_mpibcast( yvdata, nelms, infpar%master, infpar%lg_comm )
  
  call time_start_phase(PHASE_WORK)

  ! Quantities, that need to be defined but not set
  ! ********************************************************
  nelms=nbas*nbas
  call mem_alloc( df, nelms )
  call mem_alloc(  f, nelms )

  nelms=nocc*nocc
  call mem_alloc( ppf, nelms )

  nelms=nvirt*nocc
  call mem_alloc( pqf, nelms )
  call mem_alloc( qpf, nelms )
  call mem_alloc( om1, nelms )

  nelms=nvirt*nvirt
  call mem_alloc( qqf, nelms )

  ! Calculate contribution to integrals/amplitudes for slave
  ! ********************************************************
  select case (res_nr)
  case( RN_RESIDUAL_INT_DRIVEN )
     call get_ccsd_residual_integral_driven(ccmodel,df,om2,t2,f,govov,nocc,nvirt,&
        ppf,qqf,pqf,qpf,xodata,xvdata,yodata,yvdata,nbas,MyLsItem,om1,iter,local)
  case( RN_YET_ANOTHER_RES )
     call yet_another_ccsd_residual(ccmodel,df,om2,t2,f,govov,nocc,nvirt,&
        ppf,qqf,pqf,qpf,xodata,xvdata,yodata,yvdata,nbas,MyLsItem,om1,iter,local)
  end select

  ! FREE EVERYTHING
  ! ***************
  if(local)then
     call tensor_free(om2)
     call tensor_free(t2)
     call tensor_free(govov)
  else
     !FIXME : IS THIS DEALLOC CORRECT AND IF YES DO IT SOMEWHERE ELSE
     call memory_deallocate_tensor_dense(om2)
     call memory_deallocate_tensor_dense(t2)
  endif

  call mem_dealloc(      f )
  call mem_dealloc(     df )
  call mem_dealloc(    ppf )
  call mem_dealloc(    pqf )
  call mem_dealloc(    qpf )
  call mem_dealloc(    qqf )
  call mem_dealloc(    om1 )
  call mem_dealloc( xodata )
  call mem_dealloc( yodata )
  call mem_dealloc( yvdata )
  call mem_dealloc( xvdata )

  call ls_free(MyLsItem)

end subroutine ccsd_data_preparation

subroutine calculate_E2_and_permute_slave()
  use precision
  use lstiming
  use dec_typedef_module
  use typedeftype,only:lsitem,tensor
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
  integer :: no,nv,nb,s,ccmodel
  type(tensor) :: t2,omega2
  real(realk),pointer :: w1(:)
  logical :: lo
  integer(kind=8) :: o2v2
  real(realk) :: time_E2_work,time_E2_comm

  time_E2_work = 0.0E0_realk
  time_E2_comm = 0.0E0_realk


  call time_start_phase(PHASE_COMM)

  call share_E2_with_slaves(ccmodel,ppf,qqf,t2,xo,yv,Gbi,Had,no,nv,nb,omega2,s,lo)

  call time_start_phase(PHASE_WORK)

  o2v2 = int((i8*no)*no*nv*nv,kind=8)
  call mem_alloc(ppf,no*no)
  call mem_alloc(qqf,nv*nv)

  call time_start_phase(PHASE_COMM)

  call ls_mpibcast(ppf,no*no,infpar%master,infpar%lg_comm)
  call ls_mpibcast(qqf,nv*nv,infpar%master,infpar%lg_comm)

  call time_start_phase(PHASE_WORK)

  call mem_alloc(w1,o2v2)
  call calculate_E2_and_permute(ccmodel,ppf,qqf,w1,t2,xo,yv,Gbi,Had,no,nv,nb,&
     &omega2,o2v2,s,.false.,lo,time_E2_work,time_E2_comm)

  call mem_dealloc(ppf)
  call mem_dealloc(qqf)
  call mem_dealloc(xo)
  call mem_dealloc(yv)
  call mem_dealloc(Gbi)
  call mem_dealloc(Had)
  call mem_dealloc(w1)
end subroutine calculate_E2_and_permute_slave


#ifdef MOD_UNRELEASED
!> Purpose: Intermediate routine for the slaves, they get data
!           from the local master and then call the routine to 
!           calculate MO-CCSD residual.
!
!> Author:  Pablo Baudin
!> Date:    January 2014
subroutine moccsd_data_slave()

  use daltoninfo
  use dec_typedef_module
  use ccsd_module
  use tensor_type_def_module
  use infpar_module
  use lsmpi_type
  use tensor_interface_module
  use typedeftype, only: lsitem
  use decmpi_module, only: mpi_communicate_moccsd_data

  implicit none

  !> CC model
  integer :: ccmodel
  !> MO pack integrals; amplitudes and residuals:
  integer :: nbas, nocc, nvir, iter
  type(tensor) :: pgmo_diag, pgmo_up
  type(tensor) :: govov
  type(tensor) :: t1
  type(tensor) :: omega1
  type(tensor) :: t2
  type(tensor) :: omega2

  !> Long-range correction to Fock matrix
  type(tensor) :: deltafock
  !> AO fock matrix
  type(tensor) :: fock
  !> occupied-occupied block of the t1-fock matrix
  type(tensor) :: ppfock
  !> virtual-virtual block of the t1-fock matrix
  type(tensor) :: qqfock
  !> occupied-virtual block of the t1-fock matrix
  type(tensor) :: pqfock
  !> virtual-occupied block of the t1-fock matrix
  type(tensor) :: qpfock
  !> transformation matrices from AO to t1-MO:
  real(realk), pointer :: lampo(:,:), lampv(:,:)
  real(realk), pointer :: lamho(:,:), lamhv(:,:)

  !> LS item with information needed for integrals
  type(lsitem) :: MyLsItem

  !> Batches info:
  type(MObatchInfo) :: MOinfo
  
  integer :: k
  integer(kind=long) :: nelms
  logical :: local


  call mpi_communicate_moccsd_data(ccmodel,pgmo_diag,pgmo_up,t1,t2,omega2, &
         & govov,nbas,nocc,nvir,iter,MOinfo,MyLsItem,local)
  
  !==============================================================================
  ! Initialize arrays:
  if (local) then 
    call tensor_ainit(t1    , [nvir,nocc],           2, local=local, atype='LDAR' )
    call tensor_ainit(t2    , [nvir,nvir,nocc,nocc], 4, local=local, atype='LDAR' )
    call tensor_ainit(omega2, [nvir,nvir,nocc,nocc], 4, local=local, atype='LDAR' )
    call tensor_ainit(govov , [nocc,nvir,nocc,nvir], 4, local=local, atype='LDAR' )
  else
    call memory_allocate_tensor_dense(t2)
    call memory_allocate_tensor_dense(govov)
  end if

  !split messages in 2GB parts, compare to counterpart in
  !ccsd_data_preparation

  nelms = int(i8*nvir*nvir*nocc*nocc,kind=8)
  call ls_mpibcast(t2%elm1,nelms,infpar%master,infpar%lg_comm)
  if (iter/=1) then
    call ls_mpibcast(govov%elm1,nelms,infpar%master,infpar%lg_comm)
  endif

  !==============================================================================
  ! the slave call the routine to get MO-CCSD residual:
  call get_mo_ccsd_residual(ccmodel,pgmo_diag,pgmo_up,t1,omega1,t2,omega2, &
         & govov,nbas,nocc,nvir,iter,MOinfo,MyLsItem,lampo,lampv, &
         & lamho,lamhv,deltafock,fock,ppfock,pqfock,qpfock,qqfock,local)

  ! deallocate slave stuff:
  call ls_free(MyLsItem)
  call mem_dealloc(MOinfo%dimInd1)
  call mem_dealloc(MOinfo%dimInd2)
  call mem_dealloc(MOinfo%StartInd1)
  call mem_dealloc(MOinfo%StartInd2)
  call mem_dealloc(MOinfo%dimTot)
  call mem_dealloc(MOinfo%tileInd)

  if (local) then 
    call tensor_free(t1)
    call tensor_free(t2)
    call tensor_free(omega2)
    call tensor_free(govov)
  else
    call memory_deallocate_tensor_dense(omega2)
    call memory_deallocate_tensor_dense(t2)
    call memory_deallocate_tensor_dense(govov)
  end if
end subroutine moccsd_data_slave
#endif
#endif
