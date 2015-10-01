!> @file
!> Residual and energy for RPA model
!> \author Johannes Rekkedal and Thomas Bondo
module rpa_module

  use precision
  use lstiming!, only: lstimer
  use typedeftype!, only: lsitem,lssetting
  use matrix_module!, only:matrix
  use matrix_operations!, only: mat_init, mat_zero, mat_free
  use screen_mod!, only: DECscreenITEM
  use memory_handling!, only: mem_dealloc, mem_alloc
  use dec_typedef_module
  use reorder_frontend_module
  use BUILDAOBATCH!,only:build_batchesofaos,determine_maxbatchorbitalsize,&
!       & determine_MaxOrbitals
  use screen_mod!,only: free_decscreen, DECSCREENITEM
  use integralinterfaceDEC
  use integralinterfaceMod!, only: ii_get_h1, ii_get_h1_mixed_full,&
  use tensor_parameters_module
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
  use lsmpi_module
#endif

    ! DEC DEPENDENCIES (within deccc directory)   
    ! *****************************************
#ifdef VAR_MPI
  use decmpi_module!, only: mpi_communicate_ccsd_calcdata,distribute_mpi_jobs
#endif
  use ccsd_module
    use dec_fragment_utils
    use tensor_interface_module
    use array2_simple_operations!, only: array2_init, array2_add,&
!         & array2_transpose, array2_free, array2_add_to
    use array4_simple_operations!, only: array4_init, operator(*),&
!         & array_reorder_4d, mat_transpose, array4_contract1,&
!         & array4_reorder, array4_free, array4_contract2_middle,&
!         & array4_read, array4_contract2, mat_transpose_p,mat_transpose_pl,&
!         & array4_alloc, array4_add_to, array4_scale, array4_contract3,&
!         & array4_read_file_type2, array4_write_file_type2,&
!         & array4_open_file, array4_read_file, array4_close_file,&
!         & array4_write_file
    use ccintegrals!, only: get_gmo_simple,getL,dec_fock_transformation
    use ccsd_module


    public :: RPA_residual,RPA_energy,SOSEX_contribution,RPA_multiplier,&
              &rpa_residualdeb,RPA_residual_par,&
              & get_rpa_energy_arrnew,get_sosex_cont_arrnew!,&
              !& rpa_fock_para2!,RPA_fock_para,RPA_residual_par_add,
    

    private


contains

  !\brief Calculate RPA residual for current doubles amplitudes
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_residual(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    !type(array4), intent(inout) :: omega2,t2
    !type(array4), intent(in) :: gmo
    !type(array2), intent(inout) :: pfock,qfock
    !type(array4) :: tmp
    type(tensor), intent(inout) :: omega2,t2
    type(tensor), intent(in) :: gmo
    type(tensor), intent(inout) :: pfock,qfock
    type(tensor) :: tmp
    integer, intent(in) :: nocc,nvirt
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k
    real(realk) :: starttime,stoptime
    character(tensor_MSG_LEN) :: msg


    call cpu_time(starttime)


    !call tensor_zero(omega2)
    call RPA_fock_part(omega2,t2,gmo,pfock,qfock,nocc,nvirt)
    !msg = 'Norm of fockpart'
    !call print_norm(omega2,msg)
    call RPA_residual_add(omega2,t2,gmo,nocc,nvirt)

    call cpu_time(stoptime)

  end subroutine RPA_residual



  !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_fock_part(omega2,t2,gmo,pfock,qfock,no,nv)

    implicit none
    !type(array4), intent(inout) :: omega2,t2
    !type(array4), intent(in) :: gmo
    !type(array2), intent(inout) :: pfock,qfock
    !type(array4) :: tmp
    type(tensor), intent(inout) :: omega2,t2
    type(tensor), intent(in) :: gmo
    type(tensor), intent(inout) :: pfock,qfock
    type(tensor) :: tmp
    real(realk),pointer :: w_o2v2(:)
    integer, intent(in) :: no,nv
    integer, dimension(4) :: tmp_dims
    integer(KIND=LONG) :: o2v2
    integer :: a,b,c,i,j,k,nv2,no2,o2v,v2o

    no2=no*no
    nv2=nv*nv
    o2v=no*no*nv
    v2o=nv*nv*no
    o2v2=i8*no2*nv2
    call mem_alloc(w_o2v2,no2*nv2)

    !calculate first part of doubles E term and its permutation
    ! (-1) t [a b i k] * F [k j] =+ Omega [a b i j]
    call dgemm('n','n',v2o,no,no,-1.0E0_realk,t2%elm1,v2o,pfock%elm1,no,0.0E0_realk,w_o2v2,v2o)
    !call dgemm('n','n',v2o,no,no,-1.0E0_realk,t2%elm4,v2o,pfock%elm4,no,1.0E0_realk,omega2%elm4,v2o)

    !calculate second part of doubles E term
    ! F [a c] * t [c b i j] =+ Omega [a b i j]
    call dgemm('n','n',nv,o2v,nv,1.0E0_realk,qfock%elm1,nv,t2%elm1,nv,1.0E0_realk,w_o2v2,nv)

    call tensor_convert(w_o2v2,omega2)
    call array_reorder_4d(1.0E0_realk,omega2%elm1,nv,nv,no,no,[2,1,4,3],1.0E0_realk,w_o2v2)
    call tensor_convert(w_o2v2,omega2)
    call mem_dealloc(w_o2v2)

    !call tensor_convert(t2%elm4,t_par)
    !call tensor_gather(1.0E0_realk,omega2,0.0E0_realk,w_o2v2,o2v2,oo=[2,1,4,3])
    !call tensor_scatter(1.0E0_realk,w_o2v2,1.0E0_realk,omega2,o2v2)

    !For debugging
    !call array4_add_to(omega2,2.0E0_realk,gmo)

    return
  end subroutine RPA_fock_part

  !\brief Calculate additional linear and quadratic terms of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_residual_add(omega2,t2,gmo,nocc,nvirt)

    implicit none
    !type(array4), intent(inout) :: omega2,t2
    !type(array4), intent(in) :: gmo
    type(tensor), intent(inout) :: omega2,t2
    type(tensor), intent(in) :: gmo
    integer, intent(in) :: nocc,nvirt
    type(tensor) :: Sckdl
    real(realk),pointer :: Dckbj(:),omegatmp(:)
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k,dim1
    real(realk) :: starttime,stoptime
    character(tensor_MSG_LEN) :: msg

    dim1 = nocc*nvirt
    call tensor_minit(Sckdl, [nocc,nvirt,nocc,nvirt],4,atype='TDAR')

    call array_reorder_4d(1.0E0_realk,t2%elm4,nvirt,nvirt,nocc,nocc,&
      & [3,1,4,2],0.0E0_realk,Sckdl%elm4)

    do a=1,nvirt
     do i=1,nocc
        Sckdl%elm4(i,a,i,a)=Sckdl%elm4(i,a,i,a)+1._realk
     enddo
    enddo

    call mem_alloc(Dckbj,dim1*dim1)

    call dgemm('n','n',dim1,dim1,dim1, &
         & 1.0E0_realk,Sckdl%elm4,dim1,gmo%elm1,dim1,0.0E0_realk,Dckbj,dim1)

    call mem_alloc(omegatmp,dim1*dim1)

    call dgemm('n','n',dim1,dim1,dim1, &
         & 2.0E0_realk,Dckbj,dim1,Sckdl%elm4,dim1,0.0E0_realk,omegatmp,dim1)
    
    call mem_dealloc(Dckbj)
    call tensor_free(Sckdl)


    call array_reorder_4d(1.0E0_realk,omegatmp,nocc,nvirt,nocc,nvirt,&
    & [4,2,3,1],1.0E0_realk,omega2%elm1)

    call mem_dealloc(omegatmp)


  end subroutine RPA_residual_add

  !\brief Calculate RPA multipliers
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date December 2013
  subroutine RPA_multiplier(omega2,t2_final,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2,gmo
    type(array4), intent(in) :: t2_final
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    real(realk) :: starttime,stoptime


    call cpu_time(starttime)

    !get the MP2 part of the residual
    !call getDoublesResidualMP2_simple(Omega2,t2,gmo,pfock,qfock, &
    !           & nocc,nvirt)

    call RPA_fock_multi_part(omega2,t2,gmo,pfock,qfock,nocc,nvirt)
    call RPA_multi_add(omega2,t2,t2_final,gmo,pfock,qfock,nocc,nvirt)

    !for debugging
    !call array4_add_to(omega2,2.0E0_realk,gmo)

    !adds the part not containing multipliers
    !remember to understand how to know whether we should sum
    !over OCC/VIRT particles in P or [P]
    call array4_add_to(omega2,1.0E0_realk,gmo)
    call array4_reorder(gmo,[1,4,3,2])
    call array4_add_to(omega2,-0.5E0_realk,gmo)
    call array4_reorder(gmo,[1,4,3,2])

    call cpu_time(stoptime)

    !call lsquit('RPA_residual: Needs implementation',-1)

  end subroutine RPA_multiplier

  !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date December 2013
  subroutine RPA_fock_multi_part(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    type(array4), intent(in) :: gmo
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    type(array4) :: tmp
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k


    ! In 1 and 2 i ,j should be in P and a b in [P]
    ! In 3 and 4 i ,j should be in [P] while a and b in P
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
  end subroutine RPA_fock_multi_part

  !\brief Calculate RPA residual for current doubles amplitudes
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_residualdeb(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    real(realk), intent(in) :: gmo(nvirt*nocc*nvirt*nocc)
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    type(array4) :: tmp
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k
    real(realk) :: starttime,stoptime
    character(tensor_MSG_LEN) :: msg


    call cpu_time(starttime)

    !get the MP2 part of the residual
    !call getDoublesResidualMP2_simple(Omega2,t2,gmo,pfock,qfock, &
    !           & nocc,nvirt)
    !write(*,*) 'I am now in residualdeb, and everything is ok'

    call RPA_fock_partdeb(omega2,t2,pfock,qfock,nocc,nvirt)

    call RPA_residual_adddeb(omega2,t2,gmo,nocc,nvirt)

    !call print_norm
  !  msg = 'Norm of omega2'
  !  call print_norm(omega2%val,i8*nvirt*nvirt*nocc*nocc,msg)

    !for debugging
    !call array4_add_to(omega2,2.0E0_realk,gmo)

    call cpu_time(stoptime)

    !call lsquit('RPA_residual: Needs implementation',-1)

  end subroutine RPA_residualdeb

  !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_fock_partdeb(omega2,t2,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    type(array4) :: tmp
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k
    character(tensor_MSG_LEN) :: msg


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
  end subroutine RPA_fock_partdeb

  !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
!  subroutine RPA_fock_part2(omega2,t2,pfock,qfock,nocc,nvirt)
!
!    implicit none
!    type(tensor), intent(inout) :: omega2,t2
!    type(tensor), intent(inout) :: pfock,qfock
!    integer, intent(in) :: nocc,nvirt
!    type(tensor) :: tmp
!    integer, dimension(4) :: tmp_dims
!    integer :: a,b,c,i,j,k
!    character(tensor_MSG_LEN) :: msg
!
!
!    ! 1
!    !call array2_transpose(qfock)
!    
!    tmp = tensor_minit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')
!    call tensor_contract_outer_indices_rl(-1.0E0_realk,t2,pfock,0.0E0_realk,tmp)
!    call tensor_add(omega2,1.0E0_realk,tmp)
!    ! 1
!    call tensor_contract_outer_indices_rl(1.0E0_realk,qfock,t2,0.0E0_realk,tmp)
!    call tensor_add(omega2,1.0E0_realk,tmp)
!
!    call tensor_contract_outer_indices_lr(-1.0E0_realk,pfock,t2,0.0E0_realk,tmp)
!    call tensor_add(omega2,1.0E0_realk,tmp)
!    ! 2
!    call tensor_contract_outer_indices_lr(1.0E0_realk,t2,qfock,0.0E0_realk,tmp)
!    call tensor_add(omega2,1.0E0_realk,tmp)
!
!    call tensor_free(tmp)
!
!    !call array4_contract1(t2,qfock,omega2,.false.)
!    !call array2_transpose(qfock)
!
!    !! 3
!    !call array4_reorder(t2,[4,3,2,1])
!    !tmp = array4_init([nocc,nvirt,nocc,nvirt])
!    !call array4_contract1(t2,pfock,tmp,.true.)
!    !call array4_reorder(t2,[4,3,2,1])
!    !call array4_reorder(tmp,[4,3,2,1])
!    !call array4_add_to(omega2,-1.0E0_realk,tmp)
!    !call array4_free(tmp)
!
!    !! 4
!    !call array4_reorder(t2,[2,1,3,4])
!    !tmp = array4_init([nocc,nvirt,nvirt,nocc])
!    !call array4_contract1(t2,pfock,tmp,.true.)
!    !call array4_reorder(t2,[2,1,3,4])
!    !call array4_reorder(tmp,[2,1,3,4])
!    !call array4_add_to(omega2,-1.0E0_realk,tmp)
!    !call array4_free(tmp)
!
!
!    !For debugging
!    !call array4_add_to(omega2,2.0E0_realk,gmo)
!
!    return
!  end subroutine RPA_fock_part2


  !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_fock_para(omega2,t2,iajb,oof,vvf,no,nv,local)

    implicit none
    type(tensor), intent(inout) :: t2,omega2,iajb
    logical,intent(in) :: local
    !type(array4), intent(inout) :: omega2
    integer, intent(inout) :: no,nv
    !type(array2), intent(inout) :: pfock,qfock
    type(tensor), intent(inout) :: vvf,oof
    !real(realk),pointer, intent(inout) :: pfock(:),qfock(:)
    !real(realk) intent(inout) :: pfock(no,no),qfock(nv,nv)
    !type(tensor) :: tmpt2,omegaw2
    type(tensor) :: E1,E2, Pijab_om2
    real(realk),pointer :: tmp(:,:),w1(:),w2(:),w3(:),w4(:),omegw(:),w5(:)
    real(realk), pointer :: w_o2v2(:)
    integer, dimension(4) :: tmp_dims
    integer(kind=long) :: o2v2
    integer :: ord(4)
    integer :: no2,nv2,o2v,v2o,os,vs
    integer(kind=8) :: w3size
    integer :: faip,faiv,tlp,tlv,dim1,fai1,fai2,tl1,tl2,fri,tri
    integer :: i 
    real(realk) :: tw,tc
    integer(kind=ls_mpik) me,mode,nod, nnod
    character(tensor_MSG_LEN) :: msg
    logical :: master,trafo1,trafo2,trafoi

   me   = 0
   nnod = 1

   !GET SLAVES
#ifdef VAR_MPI
   me   = infpar%lg_mynum
   nnod = infpar%lg_nodtot
   mode = MPI_MODE_NOCHECK
   if(.not.local.and.me == infpar%master)call rpa_fock_communicate_data(t2,omega2,iajb,oof,vvf,no,nv)
   omega2%access_type = AT_ALL_ACCESS
   iajb%access_type   = AT_ALL_ACCESS
   t2%access_type     = AT_ALL_ACCESS
   oof%access_type    = AT_ALL_ACCESS
   vvf%access_type    = AT_ALL_ACCESS
   if(.not.local) call tensor_lock_local_wins(omega2,'e',mode)
    write(*,*) 'Total nodes in fock_para',nnod,'from',me
#endif

   vs = t2%tdim(1)
   os = t2%tdim(3)
   no = iajb%dims(1)
   nv = iajb%dims(2)

   call tensor_ainit(E1,[nv,nv],2,tdims = [vs,vs],local=local, atype="TDAR")
   call tensor_ainit(E2,[no,no],2,tdims = [os,os],local=local, atype="TDAR")

   call tensor_cp_data(vvf,E1)
   call tensor_cp_data(oof,E2)

   ord = [1,4,2,3]
   !call tensor_contract( 1.0E0_realk,t2,vvf,[2],[2],1,0.0E0_realk,omega2,ord,force_sync=.true.)
   call tensor_contract( 1.0E0_realk,t2,E1,[2],[2],1,0.0E0_realk,omega2,ord,force_sync=.true.)

   ord = [1,2,3,4]
   !call tensor_contract(-1.0E0_realk,t2,oof,[4],[1],1,1.0E0_realk,omega2,ord,force_sync=.true.)
   call tensor_contract(-1.0E0_realk,t2,E2,[4],[1],1,1.0E0_realk,omega2,ord,force_sync=.true.)


   call tensor_ainit(Pijab_om2,omega2%dims,4,local=local,tdims=int(omega2%tdim,kind=tensor_int),&
      &atype="TDAR",fo=int(omega2%offset,kind=tensor_int))

   call tensor_free(E1)
   call tensor_free(E2)

#ifdef VAR_MPI
   if(.not.local) call tensor_lock_local_wins(Pijab_om2,'e',mode)
   if(.not.local) call tensor_unlock_local_wins(omega2)
#endif

   !INTRODUCE PERMUTATION
   ord = [2,1,4,3]
   call tensor_add(Pijab_om2,1.0E0_realk,omega2, a = 0.0E0_realk, order = ord )

#ifdef VAR_MPI
   if(.not.local) call tensor_lock_local_wins(omega2,'e',mode)
   if(.not.local) call tensor_unlock_local_wins(Pijab_om2)
#endif

   call tensor_add(omega2,1.0E0_realk,Pijab_om2)

   !ADD INTEGRAL CONTRIB
   !ord = [2,4,1,3]
   !call tensor_add(omega2,1.0E0_realk,iajb, order = ord )

#ifdef VAR_MPI
   if(.not.local) call tensor_unlock_local_wins(omega2)
   omega2%access_type = AT_MASTER_ACCESS
   iajb%access_type   = AT_MASTER_ACCESS
   t2%access_type     = AT_MASTER_ACCESS
   oof%access_type    = AT_MASTER_ACCESS
   vvf%access_type    = AT_MASTER_ACCESS
#endif

   call tensor_free(Pijab_om2)
!
!    return
  end subroutine RPA_fock_para


 !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
!  subroutine RPA_fock_para2(omega2,t2,iajb,pfock,qfock,no,nv)
!
!    implicit none
!    type(tensor), intent(inout) :: t2,omega2,iajb
!    !type(array4), intent(inout) :: omega2
!    integer, intent(inout) :: no,nv
!    !type(array2), intent(inout) :: pfock,qfock
!    type(tensor), intent(inout) :: pfock,qfock
!    !real(realk), intent(inout) :: pfock(no,no),qfock(nv,nv)
!    type(tensor) :: tmpt2,omegaw2
!    real(realk),pointer :: tmp(:,:),w1(:),w3(:),w4(:),omegw(:),w5(:)
!    real(realk), pointer :: w_o2v2(:),w2(:)
!    integer, dimension(4) :: tmp_dims
!    type(tensor) :: t_par
!    integer(kind=long) :: o2v2
!    integer :: no2,nv2,o2v,v2o
!    integer(kind=8) :: w3size,b0
!    integer :: dim1,fai1,fai2,tl1,tl2,fri,tri
!    integer :: i 
!    real(realk) :: tw,tc
!    integer(kind=ls_mpik) me,mode,nod, nnod
!    character(tensor_MSG_LEN) :: msg
!    logical :: master,lock_outside,lock_safe,local,trafo1,trafo2,trafo
!
!    master=.true.
!#ifdef VAR_MPI
!    master        = .false.
!    master        = (infpar%lg_mynum == infpar%master)
!    nnod          = infpar%lg_nodtot
!    me            = infpar%lg_mynum
!    mode          = int(MPI_MODE_NOCHECK,kind=ls_mpik)
!#endif
!
!    
!    dim1=no*nv
!
!    no2  = no**2
!    nv2  = nv**2
!    o2v2 = (i8*no2)*(i8*nv2)
!    o2v  = no2*nv 
!    v2o  = nv2*no 
!    
!
!#ifdef VAR_MPI
!     StartUpSlaves: if(master .and. infpar%lg_nodtot>1) then
!       call ls_mpibcast(RPAGETFOCK,infpar%master,infpar%lg_comm)
!       call rpa_fock_communicate_data(t2,omega2,iajb,pfock,qfock,no,nv)
!    endif StartUpSlaves
!#endif
!
!
!#ifdef VAR_MPI
!
!    !Setting transformation variables for each rank
!    !**********************************************
!    call mo_work_dist(v2o,fai1,tl1,trafo1)
!    call mo_work_dist(o2v,fai2,tl2,trafo2)
!
!
!    call mo_work_dist(nv*nv*no,fri,tri,trafo)
!    !call mem_alloc(w2,nv2*no,no)
!    !call mem_alloc(w_o2v2,nv2*no,no)
!    call mem_alloc(w_o2v2,tl1*no)
!    call mem_alloc(w2,tl1*no)
!
!
!    !t_par = tensor_ainit([nv,nv,no,no],4,atype='TDAR',local=.false.)
!    !call tensor_convert(t2%elm4,t_par)
!#ifdef VAR_MPI
!    call tensor_lock_wins(t2,'s',mode)
!    call tensor_two_dim_1batch(t2,[1,2,3,4],'g',w2,3,fai1,tl1,.true.,debug=.false.)
!
!    call tensor_unlock_wins(t2)
!#endif
!    !call tensor_gather(1.0E0_realk,t2,0.0E0_realk,w2,o2v2)
!
!    ! (-1) t [a b i k] * F [k j] =+ Omega [a b i j]
!
!    call dgemm('n','n',tl1,no,no,-1.0E0_realk,w2,tl1,pfock%elm1,no,0.0E0_realk,w_o2v2,tl1)
!
!
!#ifdef VAR_MPI
!    call tensor_lock_wins(omega2,'s',mode)
!    call tensor_two_dim_1batch(omega2,[1,2,3,4],'a',w_o2v2,3,fai1,tl1,.false.,debug=.false.)
!    !call tensor_two_dim_2batch(omega2,[1,2,3,4],'a',w_o2v2,3,fai1,tl1,.true.)
!    call tensor_unlock_wins(omega2)
!#endif
!
!    call lsmpi_barrier(infpar%lg_comm)
!
!    write(*,*) 'Done with occupied'
!      
!#endif
!
!    call mem_dealloc(w2)
!    call mem_dealloc(w_o2v2)
!    !DO ALL THINGS DEPENDING ON 2
!    call mem_alloc(w_o2v2,tl2*nv)
!    !call mem_alloc(w2,tl2*nv)
!    call mem_alloc(w2,tl2*nv)
!    !call tensor_convert(t2,w2)
!    
!#ifdef VAR_MPI
!    call tensor_lock_wins(t2,'s',mode)
!    call tensor_two_dim_2batch(t2,[1,2,3,4],'g',w2,3,fai2,tl2,.true.)
!
!    call tensor_unlock_wins(t2)
!#endif
!    ! F[a c] * t [c b i j] =+ Omega [a b i j]
!    write(*,*) 'Starting with dgemm virtual'
!
!      call dgemm('n','n',nv,tl2,nv,1.0E0_realk,qfock%elm1,nv,w2,nv,0.0E0_realk,w_o2v2,nv)
!      write(*,*) 'Done with dgemm virtual'
!
!
!#ifdef VAR_MPI
!    call tensor_lock_wins(omega2,'s',mode)
!    call tensor_two_dim_2batch(omega2,[1,2,3,4],'a',w_o2v2,3,fai2,tl2,.true.)
!    call lsmpi_barrier(infpar%lg_comm)
!    call tensor_unlock_wins(omega2)
!    !call tensor_two_dim_1batch(omega2,[2,1,4,3],'a',w_o2v2,3,fai2,tl2,.false.,debug=.false.)
!    !call tensor_two_dim_2batch(omega2,[2,1,4,3],'a',w_o2v2,3,fai2,tl2,.true.)
!    !call lsmpi_barrier(infpar%lg_comm)
!    write(*,*) 'Done with virtual'
!
!
!#endif
!
!   call mem_dealloc(w_o2v2)
!   call mem_dealloc(w2)
!
!#ifdef VAR_MPI
!   if(master) then
!     call mem_alloc(w_o2v2,no2*nv2)
!     write(*,*) 'lock omega2'
!     call tensor_lock_wins(omega2,'s',mode)
!     write(*,*) 'gather omega2'
!     call tensor_gather(1.0E0_realk,omega2,0.0E0_realk,w_o2v2,o2v2,oo=[2,1,4,3])
!     write(*,*) 'unlock omega2'
!     call tensor_unlock_wins(omega2,.true.)
!     call tensor_lock_wins(omega2,'s',mode)
!     write(*,*) 'scatter to omega2'
!     call tensor_scatter(1.0E0_realk,w_o2v2,1.0E0_realk,omega2,o2v2)
!     call tensor_unlock_wins(omega2,.true.)
!     call mem_dealloc(w_o2v2)
!   endif
!#endif
!
!   return
!  end subroutine RPA_fock_para2




  !\brief Calculate additional linear and quadratic terms of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_residual_adddeb(omega2,t2,gmo,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    real(realk), intent(in) :: gmo(nvirt*nocc*nvirt*nocc)
    integer, intent(in) :: nocc,nvirt
    type(array4) :: Sckdl,Dckbj
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k,dim1
    real(realk) :: starttime,stoptime
    character(tensor_MSG_LEN) :: msg

    !Sckdl = array4_init([nvirt,nocc,nvirt,nocc])
    Sckdl = array4_duplicate(t2)
    call array4_reorder(Sckdl,[2,1,4,3])
    Dckbj = array4_init([nocc,nvirt,nvirt,nocc])


    do a=1,nvirt
     do i=1,nocc
        Sckdl%val(i,a,i,a)=Sckdl%val(i,a,i,a)+1._realk
     enddo
    enddo

    !gmo_CKLD We need it to be g_CKDL

    dim1=nocc*nvirt
    call dgemm('n','n',dim1,dim1,dim1, &
       &  1.0E0_realk,gmo,dim1,Sckdl%val,dim1,0.0E0_realk,Dckbj%val,dim1)


    call array4_reorder(omega2,[2,1,4,3])
    call dgemm('n','n',dim1,dim1,dim1, &
        & 2.0E0_realk,Sckdl%val,dim1,Dckbj%val,dim1,1.0E0_realk,omega2%val,dim1)

    call array4_reorder(omega2,[2,1,4,3])
   

    call array4_free(Dckbj)
    call array4_free(Sckdl)


  end subroutine RPA_residual_adddeb

!\brief Calculate RPA residual for current doubles amplitudes
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_residualpar(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    !type(tensor), intent(inout) :: omega2,t2
    real(realk),pointer, intent(inout) :: gmo(:)
    !type(tensor), intent(inout) :: gmo
    type(array2), intent(inout) :: pfock,qfock
    !type(tensor) :: foo,fvv
    integer, intent(in) :: nocc,nvirt
    !real(realk),pointer :: foo(:,:),fvv(:,:)
    !real(realk), intent(inout) :: pfock(nocc,nocc),qfock(nvirt,nvirt)
    !type(array4) :: tmp
    real(realk),pointer :: w2(:)!,foo(:,:),fvv(:,:)
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k
    real(realk) :: starttime,stoptime
    type(tensor) :: om_par,t_par,t_par1
    integer :: nnod,mynum,fai,tl
    integer :: nvir,noc
    logical :: master
    character(tensor_MSG_LEN) :: msg
    type(tensor) :: Sckdl
    

    nvir=nvirt
    noc=nocc

    !call mem_alloc(foo,noc,noc)
    !call mem_alloc(fvv,nvir,nvir)

    !foo = tensor_minit([nocc,nocc],2,atype='TDAR')
    !fvv = tensor_minit([nvirt,nvirt],2,atype='TDAR',local=.false.)
    !call tensor_convert(pfock%val,foo)
    !call tensor_convert(qfock%val,fvv)
    !call mem_alloc(foo,nocc,nocc)
    !call mem_alloc(fvv,nvirt,nvirt)
    !foo(:,:)=pfock%val(:,:)
    !fvv(:,:)=qfock%val(:,:)


    call cpu_time(starttime)

    call tensor_minit(t_par, [nvirt,nvirt,nocc,nocc],4,atype='TDAR')

    call array4_reorder(t2,[1,3,2,4])
    call tensor_convert(t2%val,t_par)
    call array4_reorder(t2,[1,3,2,4])


    call RPA_fock_partdeb(omega2,t2,pfock,qfock,nocc,nvirt)
    !call RPA_fock_para(omega2,t_par,foo,fvv,noc,nvir)


    !call mem_dealloc(foo)
    !call mem_dealloc(fvv)

    !call array4_reorder(omega2,[1,3,2,4])

    !call tensor_gather(1.0E0_realk,om_par,0.0E0_realk,omega2%val,i8*nvirt*nocc*nvirt*nocc)

    !call array4_reorder(omega2,[1,3,2,4])


    call tensor_minit(Sckdl, [nvirt,nvirt,nocc,nocc],4,atype='TDAR')
    !call copy_array(t2,Sckdl)
    call array4_reorder(t2,[1,3,2,4])
    call tensor_convert(t2%val,Sckdl)
    !Sckdl = array4_duplicate(t2)
    call array4_reorder(t2,[1,3,2,4])

    do a=1,nvirt
    do i=1,nocc
    Sckdl%elm4(a,a,i,i)=Sckdl%elm4(a,a,i,i)+1._realk
    !remember to change to a,a,i,i,when t2 = array
    enddo
    enddo

    !call RPA_residual_addpar(omega2,Sckdl,gmo,noc,nvir)

    ! MPI: here you should start the slaves!!

    !call lsmpi_barrier(infpar%lg_comm)
    !call lsmpi_barrier(infpar%lg_comm)
    !print*, infpar%lg_mynum, 'par residual'
    !call lsmpi_barrier(infpar%lg_comm)
    




    call cpu_time(stoptime)

    !call lsquit('RPA_residual: Needs implementation',-1)

  end subroutine RPA_residualpar


!\brief Calculate RPA residual for current doubles amplitudes
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_residual_par(omega2,t2,iajb,oof,vvf,nocc,nvirt,local)

    implicit none
    type(tensor), intent(inout) :: omega2
    type(tensor), intent(inout) :: t2
    integer, intent(in) :: nocc,nvirt
    type(tensor), intent(inout) :: iajb
    type(tensor),intent(inout) :: oof,vvf
    logical,intent(in) :: local
    real(realk),pointer :: w2(:,:,:,:)!,foo(:,:),fvv(:,:)
    type(tensor) :: E1,E2, Pijab_om2
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k
    real(realk) :: starttime,stoptime,tcpu1,twall1,tcpu2,twall2
    type(tensor) :: gtmp
    type(tensor) :: tg
    type(tensor) :: tmptens
    integer :: nnod,mynum,fai,tl
    integer :: nv,no,dim1
    integer :: ord(4)
    integer :: no2,nv2,o2v,v2o,os,vs
    integer :: fdim1(4), sdim1(4), fdim2(4), sdim2(4)
    logical :: master
    character(tensor_MSG_LEN) :: msg
    type(tensor) :: Sckdl
    integer(kind=ls_mpik) me,mode,nod
    logical :: trafo1,trafo2,trafoi

    
    nv=nvirt
    no=nocc
    dim1=nocc*nvirt

!    call mem_alloc(foo,nocc*nocc)
!    call mem_alloc(fvv,nvirt*nvirt)



    !t_par = tensor_minit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')


    !call RPA_fock_para2(omega2,t2,pfock,qfock,noc,nvir)
    !call LSTIMER('START',tcpu1,twall1,DECinfo%output)
    !call RPA_fock_para(omega2,t2,iajb,oof,vvf,noc,nvir,local)
    !call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    !write(DECinfo%output,'(a,g20.6,a)') 'Total CPU  time used in RPA fock       :',tcpu2-tcpu1,' s'
    !write(DECinfo%output,'(a,g20.6,a)') 'Total Wall time used in RPA fock        :',twall2-twall1,' s'


   me   = 0
   nnod = 1

   !Waking up slaves
#ifdef VAR_MPI
   me   = infpar%lg_mynum
   nnod = infpar%lg_nodtot
   mode = MPI_MODE_NOCHECK
   if(.not.local.and.me == infpar%master)call rpa_residual_communicate_data(t2,omega2,iajb,oof,vvf,no,nv)
   omega2%access_type = AT_ALL_ACCESS
   iajb%access_type   = AT_ALL_ACCESS
   t2%access_type     = AT_ALL_ACCESS
   oof%access_type    = AT_ALL_ACCESS
   vvf%access_type    = AT_ALL_ACCESS
   if(.not.local) call tensor_lock_local_wins(omega2,'e',mode)
    !write(*,*) 'Total nodes in fock_para',nnod,'from',me
   os     = iajb%tdim(1)
   vs     = iajb%tdim(2)
   sdim1 = [vs,os,os,vs]
   sdim2 = [vs,os,vs,os]
#endif

   vs = t2%tdim(1)
   os = t2%tdim(3)
   no = iajb%dims(1)
   nv = iajb%dims(2)

   call tensor_ainit(E1,[nv,nv],2,tdims = [vs,vs],local=local, atype="TDAR")
   call tensor_ainit(E2,[no,no],2,tdims = [os,os],local=local, atype="TDAR")

   call tensor_cp_data(vvf,E1)
   call tensor_cp_data(oof,E2)
   !copies to E1 and E2
   !oof is the occupied-occupied part of the fock matrix
   !vvf is the virtual-virtual part of the fock matrix

   ord = [1,4,2,3]
   !call tensor_contract( 1.0E0_realk,t2,vvf,[2],[2],1,0.0E0_realk,omega2,ord,force_sync=.true.)
   !sum_c t^ac_ij*f_cb
   call tensor_contract( 1.0E0_realk,t2,E1,[2],[2],1,0.0E0_realk,omega2,ord,force_sync=.true.)

   ord = [1,2,3,4]
   !call tensor_contract(-1.0E0_realk,t2,oof,[4],[1],1,1.0E0_realk,omega2,ord,force_sync=.true.)
   !sum_k t^ab_(ik)*f_kj
   call tensor_contract(-1.0E0_realk,t2,E2,[4],[1],1,1.0E0_realk,omega2,ord,force_sync=.true.)


   call tensor_ainit(Pijab_om2,omega2%dims,4,local=local,tdims=int(omega2%tdim,kind=tensor_int),&
      &atype="TDAR",fo=int(omega2%offset,kind=tensor_int))

   call tensor_free(E1)
   call tensor_free(E2)

#ifdef VAR_MPI
   if(.not.local) call tensor_lock_local_wins(Pijab_om2,'e',mode)
   if(.not.local) call tensor_unlock_local_wins(omega2)
#endif

   !INTRODUCE PERMUTATION
   ord = [2,1,4,3]
   call tensor_add(Pijab_om2,1.0E0_realk,omega2, a = 0.0E0_realk, order = ord )

#ifdef VAR_MPI
   if(.not.local) call tensor_lock_local_wins(omega2,'e',mode)
   if(.not.local) call tensor_unlock_local_wins(Pijab_om2)
#endif

   call tensor_add(omega2,1.0E0_realk,Pijab_om2)


   call tensor_free(Pijab_om2)

   !ADD INTEGRAL iajb 
   ord = [2,4,1,3]
   call tensor_add(omega2,2.0E0_realk,iajb, order = ord )
!  write(msg,*) 'Norm of gmo'
!   call print_norm(iajb,msg)


   !ADD 2t*iajb
#ifdef VAR_MPI
   call tensor_ainit(tg,[nvirt,nocc,nocc,nvirt],4,tdims=sdim1,&
     & local=local, atype="TDAR")
   if(.not.local) call tensor_lock_local_wins(tg,'e',mode)
#else
   call tensor_ainit(tg,[nvirt,nocc,nocc,nvirt],4,local=local, atype="TDAR")
#endif


   !contraction t2*g_iajb
   ord = [1,2,3,4]
   call tensor_contract( 2.0E0_realk,t2,iajb,[2,4],[2,1],2,0.0E0_realk,tg,ord,force_sync=.true.)

#ifdef VAR_MPI
   if(.not.local) call tensor_unlock_local_wins(tg)
#endif

   !add contraction t2*g_iajb to residual
   ord = [1,4,2,3]
   call tensor_add(omega2,1.0E0_realk,tg, order = ord )

   !ADD 2t*iajb*t
   ord = [1,3,2,4]
   call tensor_contract( 1.0E0_realk,tg,t2,[3,4],[3,1],2,1.0E0_realk,omega2,ord,force_sync=.true.)



   call tensor_free(tg)

   !ADD contraction g_iajb t2 to residual
   ord = [2,3,1,4]
   call tensor_contract( 2.0E0_realk,iajb,t2,[3,4],[3,1],2,1.0E0_realk,omega2,ord,force_sync=.true.)

   
  ! write(msg,*) 'Norm of omega2'
  ! call print_norm(omega2,msg)


#ifdef VAR_MPI
   if(.not.local) call tensor_unlock_local_wins(omega2)
   omega2%access_type = AT_MASTER_ACCESS
   iajb%access_type   = AT_MASTER_ACCESS
   t2%access_type     = AT_MASTER_ACCESS
   oof%access_type    = AT_MASTER_ACCESS
   vvf%access_type    = AT_MASTER_ACCESS
#endif









    !call cpu_time(starttime)
    !call LSTIMER('START',tcpu1,twall1,DECinfo%output)
    !call tensor_cp_tiled2dense(iajb,.true.)
    !call RPA_residual_par_add(omega2,t2,iajb,no,nv,local)
    !call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    !write(DECinfo%output,'(a,g20.6,a)') 'Total CPU  time used in RPA          :',tcpu2-tcpu1,' s'
    !write(DECinfo%output,'(a,g20.6,a)') 'Total Wall time used in RPA          :',twall2-twall1,' s'


    !call cpu_time(stoptime)


  end subroutine RPA_residual_par


  !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
!  subroutine RPA_fock_partpar(omega2,t2,pfock,qfock,nocc,nvirt)
!
!    implicit none
!    type(array4), intent(inout) :: omega2,t2
!    integer, intent(in) :: nocc,nvirt
!    type(array2), intent(inout) :: qfock,pfock
!    !real(realk), intent(inout) :: qfock(nvirt,nvirt),pfock(nocc,nocc)
!    !real(realk), intent(inout),dimension(nvirt,nvirt) :: qfock
!    !real(realk), intent(inout),dimension(nocc,nocc) :: pfock
!    type(array4) :: tmp
!    integer, dimension(4) :: tmp_dims
!    integer :: a,b,c,i,j,k,noc,nvir
!    nvir =nvirt
!    noc = nocc
!
!
!!#ifdef VAR_MPI
!!    StartUpSlaves: if(master .and. infpar%lg_nodtot>1) then
!!      call ls_mpibcast(RPAGETRESIDUAL,infpar%master,infpar%lg_comm)
!!      call rpa_fock_communicate_data(t2,omega2,pfock,qfock,nvir,noc)
!!    endif StartUpSlaves
!!#endif
!!
!!    call mo_work_dist(nvirt*nocc,fai,tl)
!!    call mem_alloc(w2,tl*nocc*nvirt*nocc)
!!
!!    t_par = tensor_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')
!!
!    call array4_reorder(t2,[1,3,2,4])
!!    call tensor_convert(t2%val,t_par)
!
!!    call tensor_two_dim_1batch(t_par,[1,3,2,4],'g',w2,2,fai,tl,.false.,debug=.true.)
!
!    ! 1
!    call array2_transpose(qfock)
!    call array4_reorder(t2,[3,4,1,2])
!    tmp = array4_init([nvirt,nocc,nvirt,nocc])
!    call array4_contract1(t2,qfock,tmp,.true.)
!    call array4_reorder(tmp,[3,4,1,2])
!    call array4_add_to(omega2,1.0E0_realk,tmp)
!    call array4_free(tmp)
!    call array4_reorder(t2,[3,4,1,2])
!
!    ! 2
!    call array4_contract1(t2,qfock,omega2,.false.)
!    call array2_transpose(qfock)
!
!    ! 3
!    call array4_reorder(t2,[4,3,2,1])
!    tmp = array4_init([nocc,nvirt,nocc,nvirt])
!    call array4_contract1(t2,pfock,tmp,.true.)
!    call array4_reorder(t2,[4,3,2,1])
!    call array4_reorder(tmp,[4,3,2,1])
!    call array4_add_to(omega2,-1.0E0_realk,tmp)
!    call array4_free(tmp)
!
!    ! 4
!    call array4_reorder(t2,[2,1,3,4])
!    tmp = array4_init([nocc,nvirt,nvirt,nocc])
!    call array4_contract1(t2,pfock,tmp,.true.)
!    call array4_reorder(t2,[2,1,3,4])
!    call array4_reorder(tmp,[2,1,3,4])
!    call array4_add_to(omega2,-1.0E0_realk,tmp)
!    call array4_free(tmp)
!
!
!    !For debugging
!    !call array4_add_to(omega2,2.0E0_realk,gmo)
!
!    return
!  end subroutine RPA_fock_partpar


  !\brief Calculate additional linear and quadratic terms of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
!  subroutine RPA_residual_par_add(omega2,t2,iajb,nocc,nvirt,local)
!
!    implicit none
!    !type(array4), intent(inout) :: omega2!,u2
!    type(tensor), intent(inout) :: t2 ,omega2
!    !real(realk), intent(inout),pointer :: gmo(:)
!    type(tensor), intent(inout):: iajb
!    logical,intent(in) :: local
!    integer,intent(inout) :: nocc,nvirt
!    type(tensor) :: tg
!    type(tensor) :: tmptens
!    integer, dimension(4) :: tmp_dims
!    integer :: a,b,c,i,j,k,dim1,no,nv
!    integer :: fai,tl,mynum,nnod,me
!    integer :: ord(4)
!    integer(kind=ls_mpik) :: mode
!    real(realk) :: starttime,stoptime
!    real(realk),pointer :: w2(:),w3(:),omegw(:),w4(:)
!    character(tensor_MSG_LEN) :: msg
!    logical :: master,trafo
!    integer :: fdim1(4), sdim1(4), fdim2(4), sdim2(4)
!    integer :: os, vs
!
!   me   = 0
!   nnod = 1
!
!   !GET SLAVES
!#ifdef VAR_MPI
!   me   = infpar%lg_mynum
!   nnod = infpar%lg_nodtot
!   mode = MPI_MODE_NOCHECK
!   if(.not.local.and.me == infpar%master)call rpa_res_communicate_data(iajb,t2,omega2,nvirt,nocc)
!   omega2%access_type = AT_ALL_ACCESS
!   iajb%access_type   = AT_ALL_ACCESS
!   t2%access_type     = AT_ALL_ACCESS
!   if(.not.local) call tensor_lock_local_wins(omega2,'e',mode)
!   os     = iajb%tdim(1)
!   vs     = iajb%tdim(2)
!   sdim1 = [vs,os,os,vs]
!   sdim2 = [vs,os,vs,os]
!    write(*,*) 'Total nodes in rest of rpa residual',nnod,'from',me
!#endif
!
!   !ADD INTEGRAL iajb 
!   ord = [2,4,1,3]
!   call tensor_add(omega2,2.0E0_realk,iajb, order = ord )
!!   write(msg,*) 'Norm of gmo'
!!   call print_norm(iajb,msg)
!
!
!   !ADD 2t*iajb
!#ifdef VAR_MPI
!   call tensor_ainit(tg,[nvirt,nocc,nocc,nvirt],4,tdims=sdim1,&
!     & local=local, atype="TDAR")
!   if(.not.local) call tensor_lock_local_wins(tg,'e',mode)
!#else
!   call tensor_ainit(tg,[nvirt,nocc,nocc,nvirt],4,local=local, atype="TDAR")
!#endif
!
!
!   ord = [1,2,3,4]
!   call tensor_contract( 2.0E0_realk,t2,iajb,[2,4],[2,1],2,0.0E0_realk,tg,ord,force_sync=.true.)
!
!#ifdef VAR_MPI
!   if(.not.local) call tensor_unlock_local_wins(tg)
!#endif
!
!   ord = [1,4,2,3]
!   call tensor_add(omega2,1.0E0_realk,tg, order = ord )
!
!   !ADD 2t*iajb*t
!   ord = [1,3,2,4]
!   call tensor_contract( 1.0E0_realk,tg,t2,[3,4],[3,1],2,1.0E0_realk,omega2,ord,force_sync=.true.)
!
!
!
!   call tensor_free(tg)
!
!   !ADD 2iajb*t
!   ord = [2,3,1,4]
!   call tensor_contract( 2.0E0_realk,iajb,t2,[3,4],[3,1],2,1.0E0_realk,omega2,ord,force_sync=.true.)
!
!   
!  ! write(msg,*) 'Norm of omega2'
!  ! call print_norm(omega2,msg)
!
!#ifdef VAR_MPI
!   if(.not.local) call tensor_unlock_local_wins(omega2)
!   omega2%access_type = AT_MASTER_ACCESS
!   iajb%access_type   = AT_MASTER_ACCESS
!   t2%access_type     = AT_MASTER_ACCESS
!#endif
!
!
!  end subroutine RPA_residual_par_add


  !\brief Calculate additional linear and quadratic terms of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
!  subroutine RPA_residual_addpar(omega2,u2,gmo,nocc,nvirt)
!
!    implicit none
!    type(array4), intent(inout) :: omega2!,u2
!    type(tensor), intent(inout) :: u2 !,omega2
!    real(realk), intent(inout),pointer :: gmo(:)
!    !type(tensor), intent(inout):: gmo
!    integer,intent(inout) :: nocc,nvirt
!    !type(array4) :: Sckdl,Dckbj
!    type(tensor) :: Sckdl,Dckbj
!    type(tensor) :: t_par,omegaw1
!    integer, dimension(4) :: tmp_dims
!    integer :: a,b,c,i,j,k,dim1
!    integer :: fai,tl,mynum,nnod
!    real(realk) :: starttime,stoptime
!    real(realk),pointer :: w2(:),w3(:),omegw(:),w4(:)
!    character(tensor_MSG_LEN) :: msg
!    logical :: master
!
!    master=.true.
!#ifdef VAR_MPI
!    master        = .false.
!    mynum         = infpar%lg_mynum
!    master        = (infpar%lg_mynum == 0)
!    nnod          = infpar%lg_nodtot
!#endif
!
!
!    dim1=nocc*nvirt
!#ifdef VAR_MPI
!    StartUpSlaves: if(master .and. infpar%lg_nodtot>1) then
!      call ls_mpibcast(RPAGETRESIDUAL,infpar%master,infpar%lg_comm)
!      call rpa_res_communicate_data(gmo,u2,omega2,nvirt,nocc)
!    endif StartUpSlaves
!#endif
!   
!
!    omegaw1 = tensor_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR',local=.false.)
!    !write(*,*) 'omegaw1 tdim=',omegaw1%tdim
!    !write(*,*) 'omegaw1 dims=',omegaw1%dims
!    !write(*,*) 'omegaw1 tile=',omegaw1%ntiles
!
!!#ifdef VAR_MPI
!!    if(master) then
!!      call array4_reorder(omega2,[1,3,2,4])
!!      call tensor_convert(omega2%val,omegaw1)
!!      !write(msg,*) 'Norm of omegaw1',infpar%lg_mynum
!!      write(*,*) 'omegaw1 tdim=',omegaw1%tdim
!!      !call print_norm(omegaw1,msg)
!!    endif
!!#else
!!#endif
!
!    call array4_reorder(omega2,[1,3,2,4])
!    call tensor_convert(omega2%val,omegaw1)
!!#endif
!
!
!    call mo_work_dist(nvirt*nocc,fai,tl)
!    !call mem_alloc(w2,tl*nocc*nvirt)
!!#ifdef VAR_MPI
!!    write(*,*) 'size of tile',tl,fai,infpar%lg_mynum,nvirt*nocc!*nvirt*nocc
!!#endif
!    call mem_alloc(w2,tl*nocc*nvirt)
!    call mem_alloc(w3,tl*nocc*nvirt)
!    call mem_alloc(w4,nocc*nvirt*nocc*nvirt)
!    call mem_alloc(omegw,tl*nocc*nvirt)
!    !call mem_alloc(omegw,tl*tl)
!
!    t_par = tensor_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR',local=.false.)
!    !call array4_reorder(u2,[1,3,2,4])
!    !call copy_array(u2,t_par)
!    call tensor_convert(u2%elm4,t_par)
!
!    call tensor_two_dim_1batch(t_par,[4,2,3,1],'g',w2,2,fai,tl,.false.,debug=.true.)
!
!
!    omegw=0.0_realk
!
!    !When fock part is parallelized instead of zero 1.0_realk
!    call dgemm('n','n',tl,dim1,dim1, &
!         1.0E0_realk,w2,tl,gmo,dim1,0.0E0_realk,w3,tl)
!!#ifdef VAR_MPI
!!    write(msg,*) 'Norm of w3',infpar%lg_mynum
!!#else
!!    write(msg,*) 'Norm of w3'
!!#endif
!!    call print_norm(w3,i8*tl*nvirt*nocc,msg)
!
!!#ifdef VAR_MPI
!    !call tensor_gather(1.0E0_realk,u2,0.0E0_realk,w4,i8*dim1*dim1,oo=[1,3,2,4])
!    !call tensor_gather(1.0E0_realk,t_par,0.0E0_realk,w4,i8*dim1*dim1,oo=[1,3,2,4])
!    call tensor_gather(1.0E0_realk,t_par,0.0E0_realk,w4,i8*dim1*dim1,oo=[4,2,3,1])
!!    write(msg,*) 'Norm of w4',infpar%lg_mynum
!!#else
!!    write(msg,*) 'Norm of w4'
!
!!#endif
!    !call sleep(1)
!    !call lsmpi_barrier(infpar%lg_comm)
!!    call print_norm(w4,i8*dim1*dim1,msg)
!    !call sleep(1)
!    !call lsmpi_barrier(infpar%lg_comm)
!    !stop
!    
!    call dgemm('n','n',tl,dim1,dim1, &
!         2.0E0_realk,w3,tl,w4,dim1,0.0E0_realk,omegw,tl)
!
!!#ifdef VAR_MPI
!!    write(msg,*) 'Norm of omegw',infpar%lg_mynum
!    !call sleep(1)
!    !call lsmpi_barrier(infpar%lg_comm)
!!    call print_norm(omegw,i8*tl*nvirt*nocc,msg)
!    !call sleep(1)
!    !call lsmpi_barrier(infpar%lg_comm)
!    !stop
!!#endif
!    
!
!    !call tensor_two_dim_1batch(omegaw1,[1,3,2,4],'a',omegw,2,fai,tl,.false.,debug=.true.)
!    call tensor_two_dim_1batch(omegaw1,[4,2,3,1],'a',omegw,2,fai,tl,.false.,debug=.true.)
!
!#ifdef VAR_MPI
!
!    call lsmpi_barrier(infpar%lg_comm)
!
!#endif
!
!!#ifdef VAR_MPI
!    if(master) then
!
!      call tensor_gather(1.0E0_realk,omegaw1,0.0E0_realk,omega2%val,i8*nvirt*nocc*nvirt*nocc)
!      
!      call array4_reorder(omega2,[1,3,2,4])
!
!     ! write(msg,*) 'Norm of omega2'
!     ! call print_norm(omega2%val,i8*dim1*dim1,msg)
!!      write(*,*) 'checkpoint 3',infpar%lg_mynum
!      !call array4_reorder(omega2,[2,1,4,3])
!    endif
!
!
!    call tensor_free(t_par)
!    call tensor_free(omegaw1)
!    call mem_dealloc(w2)
!    call mem_dealloc(w3)
!    call mem_dealloc(w4)
!    call mem_dealloc(omegw)
!
!
!  end subroutine RPA_residual_addpar


  !\brief Calculate additional linear and quadratic terms of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_multi_add(omega2,t2,t2_final,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2,t2
    type(array4), intent(in) :: gmo,t2_final
    type(array2), intent(inout) :: pfock,qfock
    integer, intent(in) :: nocc,nvirt
    type(array4) :: Sckdl
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k,dim1
    real(realk) :: starttime,stoptime

    Sckdl = array4_init([nvirt,nocc,nvirt,nocc])



    dim1=nocc*nvirt
    call dgemm('n','n',dim1,dim1,dim1, &
        & 1.0E0_realk,gmo%val,dim1,t2%val,dim1,0.0E0_realk,omega2%val,dim1)

    call dgemm('n','n',dim1,dim1,dim1, &
        & 1.0E0_realk,t2%val,dim1,gmo%val,dim1,1.0E0_realk,omega2%val,dim1)
    
       !gmo*t2*m2
    call dgemm('n','n',dim1,dim1,dim1, &
        & 1.0E0_realk,gmo%val,dim1,t2_final%val,dim1,0.0E0_realk,Sckdl%val,dim1)

    call dgemm('n','n',dim1,dim1,dim1, &
        & 1.0E0_realk,Sckdl%val,dim1,t2%val,dim1,1.0E0_realk,omega2%val,dim1)

       !m2*t2*gmo
    call dgemm('n','n',dim1,dim1,dim1, &
        & 1.0E0_realk,t2%val,dim1,t2_final%val,dim1,0.0E0_realk,Sckdl%val,dim1)

    call dgemm('n','n',dim1,dim1,dim1, &
        & 1.0E0_realk,Sckdl%val,dim1,gmo%val,dim1,1.0E0_realk,omega2%val,dim1)

    call array4_free(Sckdl)

    return

  end subroutine RPA_multi_add
 

  function get_rpa_energy_arrnew(t2,gmo,nocc,nvirt) result(ecorr)

    implicit none
    type(tensor), intent(in) :: t2
    type(tensor), intent(inout) :: gmo
    integer, intent(in) :: nocc,nvirt
    real(realk) :: ecorr,ecorr_d
    integer :: a,i,b,j

    ecorr = 0.0E0_realk
    ecorr_d = 0.0E0_realk

    if(t2%itype==TT_DENSE.and.gmo%itype==TT_DENSE)then
      !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) PRIVATE(i,a,j,b) SHARED(nocc,&
      !$OMP nvirt,t2,gmo) REDUCTION(+:ecorr_d)
      do j=1,nocc
      do b=1,nvirt
      do i=1,nocc
      do a=1,nvirt
      ecorr_d = ecorr_d + t2%elm4(a,b,i,j)* &
       & (1.0E0_realk*gmo%elm4(i,a,j,b))
      end do
      end do
      end do
      end do
      !$OMP END PARALLEL DO

      ecorr = ecorr_d

    else if(t2%itype==TT_TILED_DIST.and.gmo%itype==TT_TILED_DIST)then

      ecorr=get_rpa_energy_parallel(t2,gmo)

    endif


  end function get_rpa_energy_arrnew


  function get_sosex_cont_arrnew(t2,gmo,nocc,nvirt) result(ecorr)

    implicit none
    type(tensor), intent(in) :: t2
    type(tensor), intent(inout) :: gmo
    integer, intent(in) :: nocc,nvirt
    real(realk) :: ecorr,ecorr_d
    integer :: a,i,b,j

    ecorr = 0.0E0_realk
    ecorr_d = 0.0E0_realk

    if(t2%itype==TT_DENSE.and.gmo%itype==TT_DENSE)then
      !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) PRIVATE(i,a,j,b) SHARED(nocc,&
      !$OMP nvirt,t2,gmo) REDUCTION(+:ecorr_d)
      do j=1,nocc
      do b=1,nvirt
      do i=1,nocc
      do a=1,nvirt
      ecorr_d = ecorr_d + t2%elm4(a,b,i,j)* &
       & (-0.5E0_realk*gmo%elm4(i,b,j,a))
      end do
      end do
      end do
      end do
      !$OMP END PARALLEL DO

      ecorr = ecorr_d

    else if(t2%itype==TT_TILED_DIST.and.gmo%itype==TT_TILED_DIST)then

      ecorr=get_sosex_cont_parallel(t2,gmo)

    endif


  end function get_sosex_cont_arrnew


!  function get_rpa_energy_parallel(t2,gmo) result(Ec)
!    implicit none
!    !> two electron integrals in the mo-basis
!    type(tensor), intent(inout) :: gmo
!    !> doubles amplitudes
!    type(tensor), intent(in) :: t2
!    !> on return Ec contains the correlation energy
!    real(realk) :: E1,E2,Ec
!    real(realk),pointer :: t(:,:,:,:)
!    integer :: lt,i,j,a,b,o(t2%mode),da,db,di,dj
!
!#ifdef VAR_MPI
!    !Get the slaves to this routine
!    if(infpar%lg_mynum==infpar%master)then
!      call pdm_tensor_sync(infpar%lg_comm,JOB_GET_MP2_ENERGY,t2,gmo)
!    endif
!    call memory_allocate_tensor_dense(gmo)
!    call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)
!
!    E2=0.0E0_realk
!    Ec=0.0E0_realk
!    do lt=1,t2%nlti
!      call ass_D1to4(t2%ti(lt)%t,t,t2%ti(lt)%d)
!      !get offset for global indices
!      call get_midx(t2%ti(lt)%gt,o,t2%ntpm,t2%mode)
!      do j=1,t2%mode
!        o(j)=(o(j)-1)*t2%tdim(j)
!      enddo
!
!      da = t2%ti(lt)%d(1)
!      db = t2%ti(lt)%d(2)
!      di = t2%ti(lt)%d(3)
!      dj = t2%ti(lt)%d(4)
!      !count over local indices
!      do j=1,dj
!        do i=1,di
!          do b=1,db
!            do a=1,da
!     
!              E2 = E2 + t(a,b,i,j)*&
!              & (1.0E0_realk*  gmo%elm4(i+o(3),a+o(1),j+o(4),b+o(2)))
!   
!            enddo 
!          enddo
!        enddo
!      enddo
!      nullify(t)
!    enddo
!
!    call tensor_deallocate_dense(gmo)
!    
!    call lsmpi_local_reduction(E2,infpar%master)
!
!    Ec = E2
!#else
!    Ec = 0.0E0_realk
!#endif
!  end function get_rpa_energy_parallel
!
!
!  function get_sosex_cont_parallel(t2,gmo) result(Ec)
!    implicit none
!    !> two electron integrals in the mo-basis
!    type(tensor), intent(inout) :: gmo
!    !> doubles amplitudes
!    type(tensor), intent(in) :: t2
!    !> on return Ec contains the correlation energy
!    real(realk) :: E2,Ec
!    real(realk),pointer :: t(:,:,:,:)
!    integer :: lt,i,j,a,b,o(t2%mode),da,db,di,dj
!
!#ifdef VAR_MPI
!    !Get the slaves to this routine
!    if(infpar%lg_mynum==infpar%master)then
!      call pdm_tensor_sync(infpar%lg_comm,JOB_GET_MP2_ENERGY,t2,gmo)
!    endif
!    call memory_allocate_tensor_dense(gmo)
!    call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)
!
!    E2=0.0E0_realk
!    Ec=0.0E0_realk
!    do lt=1,t2%nlti
!      call ass_D1to4(t2%ti(lt)%t,t,t2%ti(lt)%d)
!      !get offset for global indices
!      call get_midx(t2%ti(lt)%gt,o,t2%ntpm,t2%mode)
!      do j=1,t2%mode
!        o(j)=(o(j)-1)*t2%tdim(j)
!      enddo
!
!      da = t2%ti(lt)%d(1)
!      db = t2%ti(lt)%d(2)
!      di = t2%ti(lt)%d(3)
!      dj = t2%ti(lt)%d(4)
!      do j=1,dj
!        do i=1,di
!          do b=1,db
!            do a=1,da
!     
!              E2 = E2 + t(a,b,i,j)*&
!              & (-0.5E0_realk* gmo%elm4(i+o(3),b+o(2),j+o(4),a+o(1)) )
!   
!            enddo 
!          enddo
!        enddo
!      enddo
!      nullify(t)
!    enddo
!
!    call tensor_deallocate_dense(gmo)
!    
!    call lsmpi_local_reduction(E2,infpar%master)
!
!    Ec = E2
!#else
!    Ec = 0.0E0_realk
!#endif
!  end function get_sosex_cont_parallel



  !\brief Calculate RPA energy
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  function RPA_energy(t2,gmo) result(energy)
    implicit none
    type(array4) :: J,X,gmo
    type(array4), intent(in) :: t2
    real(realk) :: energy

    !write(*,*) 'In rpa_energy'
    !Test for understanding the structure, 
    J = array4_duplicate(gmo)
    !call array4_scale(J,2.E0_realk)
    energy=t2*J

    call array4_free(J)

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
    call array4_free(J)

    !call lsquit('rpa_energy: Needs implementation',-1)

  end function SOSEX_contribution

end module rpa_module

#ifdef VAR_MPI
!subroutine rpa_res_slave()
!  use dec_typedef_module
!  use typedeftype,only:lsitem
!  use tensor_interface_module
!  use decmpi_module,only:rpa_res_communicate_data
!  use infpar_module
!  use rpa_module
!  implicit none
!  !> number of orbitals:
!  type(tensor) :: omega2,t2
!  type(tensor) :: gmo
!  logical :: local
!  !type(array4) :: omega2!,t2
!  !real(realk),pointer :: gmo(:)!,t2(:,:,:,:)
!  type(array2)  :: pfock,qfock
!  integer :: nbas, nocc, nvirt
!  !> how to pack integrals:
!  
!  !print*, infpar%lg_mynum,'rpa_res_slave'
!  call rpa_res_communicate_data(gmo,t2,omega2,nvirt,nocc)
!  local = .false.
!  call RPA_residual_par_add(omega2,t2,gmo,nocc,nvirt,local)
!
!end subroutine rpa_res_slave

subroutine rpa_fock_slave()
  use dec_typedef_module
  use typedeftype,only:lsitem
  use tensor_interface_module
  use decmpi_module,only:rpa_fock_communicate_data
  use infpar_module
  use rpa_module
  implicit none
  !> number of orbitals:
  type(tensor) :: t2
  type(tensor) :: omega2,iajb!,t2
  type(tensor)  :: oof,vvf
  real(realk),pointer :: gmo(:)
  logical :: local
  !real(realk),pointer  :: pfock(:),qfock(:)
  integer :: nbas, nocc, nvirt
  !> how to pack integrals:
  !print*, infpar%lg_mynum,'rpa_fock_slave'
  call rpa_fock_communicate_data(t2,omega2,iajb,oof,vvf,nocc,nvirt)
  local = .false.
  !call RPA_fock_para(omega2,t2,iajb,oof,vvf,nocc,nvirt,local)
  call RPA_residual_par(omega2,t2,iajb,oof,vvf,nocc,nvirt,local)
  !call RPA_fock_para2(omega2,t2,pfock,qfock,nocc,nvirt)

end subroutine rpa_fock_slave
#endif
