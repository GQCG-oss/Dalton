!> @file
!> Residual and energy for RPA model
!> \author Johannes Rekkedal and Thomas Bondo
module rpa_module

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

    ! DEC DEPENDENCIES (within deccc directory)   
    ! *****************************************
#ifdef VAR_LSMPI
  use decmpi_module!, only: mpi_communicate_ccsd_calcdata,distribute_mpi_jobs
#endif
    use dec_fragment_utils
    use array_memory_manager!, only: DENSE,TILED,TILED_DIST,SCALAPACK,&
!         & NO_PDM,MASTER_INIT,REPLICATED,ALL_INIT,ass_1to3,ass_1to2,&
!         & ass_1to4,ass_2to1,&
!         & ass_4to1,ARR_MSG_LEN
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
    use dec_pdm_module!, only: precondition_doubles_parallel
    use array_operations!, only: array_init, array_change_atype_to_rep,&
!         & array_change_atype_to_d,print_norm
    use ccintegrals!, only: get_gmo_simple,getL,dec_fock_transformation
    use ccsd_module


    public :: RPA_residual,RPA_energy,SOSEX_contribution

    private


contains

  !\brief Calculate RPA residual for current doubles amplitudes
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_residual(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    type(array4), intent(in) :: gmo
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    type(array4) :: tmp
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k
    real(realk) :: starttime,stoptime


    call cpu_time(starttime)

    !get the MP2 part of the residual
    !call getDoublesResidualMP2_simple(Omega2,t2,gmo,pfock,qfock, &
    !           & nocc,nvirt)

    call RPA_fock_part(omega2,t2,gmo,pfock,qfock,nocc,nvirt)
    call RPA_residual_add(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    !for debugging
    !call array4_add_to(omega2,2.0E0_realk,gmo)

    call cpu_time(stoptime)

    !call lsquit('RPA_residual: Needs implementation',-1)

  end subroutine RPA_residual



  !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_fock_part(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

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


    !For debugging
    !call array4_add_to(omega2,2.0E0_realk,gmo)

    return
  end subroutine RPA_fock_part

  !\brief Calculate additional linear and quadratic terms of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_residual_add(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    type(array4), intent(in) :: gmo
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    type(array4) :: Sckdl,Dckbj
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k,dim1
    real(realk) :: starttime,stoptime

    Sckdl = array4_init([nvirt,nocc,nvirt,nocc])
    Sckdl = array4_duplicate(t2)
    Dckbj = array4_init([nvirt,nocc,nvirt,nocc])

    do a=1,nvirt
     do i=1,nocc
        Sckdl%val(a,i,a,i)=Sckdl%val(a,i,a,i)+1._realk
     enddo
    enddo


    dim1=nocc*nvirt
    call dgemm('n','n',dim1,dim1,dim1, &
         1.0E0_realk,gmo%val,dim1,Sckdl%val,dim1,0.0E0_realk,Dckbj%val,dim1)

    call dgemm('n','n',dim1,dim1,dim1, &
         2.0E0_realk,Sckdl%val,dim1,Dckbj%val,dim1,1.0E0_realk,omega2%val,dim1)
    

    call array4_free(Dckbj)
    call array4_free(Sckdl)

    return

  end subroutine RPA_residual_add








  !\brief Calculate RPA energy
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  function RPA_energy(t2,gmo) result(energy)
    implicit none
    type(array4) :: J,X,gmo
    type(array4), intent(in) :: t2
    real(realk) :: energy

    !Test for understanding the structure, 
    !gmo would be g_aibj
    J = array4_duplicate(gmo)
    !call array4_scale(J,2.E0_realk)
    energy=t2*J

    !call lsquit('rpa_energy: Needs implementation',-1)

  end function RPA_energy

  !\brief Calculate RPA energy
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  function SOSEX_contribution(t2,gmo) result(energy)
    implicit none
    type(array4) :: J,gmo
    type(array4), intent(in) :: t2
    real(realk) :: energy

    !Test for understanding the structure, 
    !gmo would be g_aibj
    call array4_reorder(gmo,[1,4,3,2])
    J = array4_duplicate(gmo)
    call array4_scale(J,-0.5E0_realk)
    energy=t2*J
    call array4_reorder(gmo,[1,4,3,2])

    !call lsquit('rpa_energy: Needs implementation',-1)

  end function SOSEX_contribution

end module rpa_module
