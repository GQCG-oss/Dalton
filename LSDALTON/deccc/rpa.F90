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
  use integralinterfaceDEC
  use integralinterfaceMod!, only: ii_get_h1, ii_get_h1_mixed_full,&
  use ccsd_module
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
  !use integralparameters!, only: AORdefault
#endif

    ! DEC DEPENDENCIES (within deccc directory)   
    ! *****************************************
#ifdef VAR_MPI
  use decmpi_module!, only: mpi_communicate_ccsd_calcdata,distribute_mpi_jobs
#endif
    use dec_fragment_utils
    use tensor_interface_module
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
    use ccintegrals!, only: get_gmo_simple,getL,dec_fock_transformation
    use ccsd_module


    public :: RPA_residual,RPA_energy,SOSEX_contribution,RPA_multiplier,&
              &rpa_residualdeb,RPA_residual_par_add,RPA_residual_par,&
              & RPA_fock_para,get_rpa_energy_arrnew,get_sosex_cont_arrnew,&
              & rpa_fock_para2
    

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

    !Sckdl = array4_init([nvirt,nocc,nvirt,nocc])
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
    character(ARR_MSG_LEN) :: msg


    call cpu_time(starttime)

    !get the MP2 part of the residual
    !call getDoublesResidualMP2_simple(Omega2,t2,gmo,pfock,qfock, &
    !           & nocc,nvirt)
    !write(*,*) 'I am now in residualdeb, and everything is ok'
    msg = 'Norm of t2'
    call print_norm(t2%val,i8*nvirt*nvirt*nocc*nocc,msg)

    call RPA_fock_partdeb(omega2,t2,pfock,qfock,nocc,nvirt)
   ! msg = 'Norm of fockpart'
   ! call print_norm(omega2%val,i8*nvirt*nvirt*nocc*nocc,msg)

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
    character(ARR_MSG_LEN) :: msg


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
!    type(array), intent(inout) :: omega2,t2
!    type(array), intent(inout) :: pfock,qfock
!    integer, intent(in) :: nocc,nvirt
!    type(array) :: tmp
!    integer, dimension(4) :: tmp_dims
!    integer :: a,b,c,i,j,k
!    character(ARR_MSG_LEN) :: msg
!
!
!    ! 1
!    !call array2_transpose(qfock)
!    
!    tmp = array_minit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')
!    call array_contract_outer_indices_rl(-1.0E0_realk,t2,pfock,0.0E0_realk,tmp)
!    call array_add(omega2,1.0E0_realk,tmp)
!    ! 1
!    call array_contract_outer_indices_rl(1.0E0_realk,qfock,t2,0.0E0_realk,tmp)
!    call array_add(omega2,1.0E0_realk,tmp)
!
!    call array_contract_outer_indices_lr(-1.0E0_realk,pfock,t2,0.0E0_realk,tmp)
!    call array_add(omega2,1.0E0_realk,tmp)
!    ! 2
!    call array_contract_outer_indices_lr(1.0E0_realk,t2,qfock,0.0E0_realk,tmp)
!    call array_add(omega2,1.0E0_realk,tmp)
!
!    call array_free(tmp)
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
  subroutine RPA_fock_para(omega2,t2,pfock,qfock,no,nv)

    implicit none
    type(array), intent(inout) :: t2,omega2
    !type(array4), intent(inout) :: omega2
    integer, intent(inout) :: no,nv
    !type(array2), intent(inout) :: pfock,qfock
    type(array), intent(inout) :: pfock,qfock
    !real(realk),pointer, intent(inout) :: pfock(:),qfock(:)
    !real(realk) intent(inout) :: pfock(no,no),qfock(nv,nv)
    type(array) :: tmpt2,omegaw2
    real(realk),pointer :: tmp(:,:),w1(:),w2(:),w3(:),w4(:),omegw(:),w5(:)
    real(realk), pointer :: w_o2v2(:)
    integer, dimension(4) :: tmp_dims
    integer(kind=long) :: o2v2
    integer :: no2,nv2,o2v,v2o
    integer(kind=8) :: w3size
    integer :: faip,faiv,tlp,tlv,dim1,fai1,fai2,tl1,tl2,fri,tri
    integer :: i 
    real(realk) :: tw,tc
    integer(kind=ls_mpik) me,mode,nod, nnod
    character(ARR_MSG_LEN) :: msg
    logical :: master,lock_outside,lock_safe,local

    master=.true.
    local = .true.
#ifdef VAR_MPI
    master        = .false.
    master        = (infpar%lg_mynum == 0)
    nnod          = infpar%lg_nodtot
    me            = infpar%lg_mynum
    mode          = int(MPI_MODE_NOCHECK,kind=ls_mpik)
#endif

    
    dim1=no*nv

    no2  = no**2
    nv2  = nv**2
    o2v2 = (i8*no2)*(i8*nv2)
    o2v  = no2*nv 
    v2o  = nv2*no 
    call mem_alloc(w_o2v2,t2%nelms)
    
    if(local) then


      write(*,*) 'before first dgemm'
      !calculate first part of doubles E term and its permutation
      ! (-1) t [a b i k] * F [k j] =+ Omega [a b i j]
      call dgemm('n','n',v2o,no,no,-1.0E0_realk,t2%elm1,v2o,pfock%elm1,no,0.0E0_realk,w_o2v2,v2o)
      !call dgemm('n','n',v2o,no,no,-1.0E0_realk,t2%elm4,v2o,pfock%elm4,no,1.0E0_realk,omega2%elm4,v2o)
      write(*,*) 'after first dgemm'

      !calculate second part of doubles E term
      ! F [a c] * t [c b i j] =+ Omega [a b i j]
      call dgemm('n','n',nv,o2v,nv,1.0E0_realk,qfock%elm1,nv,t2%elm1,nv,1.0E0_realk,w_o2v2,nv)
      write(*,*) 'after second dgemm'

      call mem_alloc(w3,no2*nv2)

      call array_scatter(1.0E0_realk,w_o2v2,0.0E0_realk,omega2,o2v2)
      call array_gather(1.0E0_realk,omega2,0.0E0_realk,w_o2v2,o2v2,oo=[2,1,4,3])
      call array_scatter(1.0E0_realk,w_o2v2,1.0E0_realk,omega2,o2v2)


!      !INTRODUCE PERMUTATION
!#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
!      call assign_in_subblocks(w1,'=',omega2%elm1,o2v2)
!#else
!      !$OMP WORKSHARE
!      w_o2v2(1_long:o2v2) = omega2%elm1(1_long:o2v2)
!      !$OMP END WORKSHARE
!#endif
!
!      call array_reorder_4d(1.0E0_realk,w_o2v2,nv,nv,no,no,[2,1,4,3],1.0E0_realk,omega2%elm1)

   else

#ifdef VAR_MPI
     StartUpSlaves: if(master .and. infpar%lg_nodtot>1) then
       call ls_mpibcast(RPAGETFOCK,infpar%master,infpar%lg_comm)
       call rpa_fock_communicate_data(t2,omega2,pfock,qfock,nv,no)
    endif StartUpSlaves
#endif


#ifdef VAR_MPI

    omega2%access_type = ALL_ACCESS
    t2%access_type     = ALL_ACCESS
    qfock%access_type    = ALL_ACCESS
    pfock%access_type    = ALL_ACCESS
    nnod               = infpar%lg_nodtot
    me                 = infpar%lg_mynum
    mode               = int(MPI_MODE_NOCHECK,kind=ls_mpik)
    lock_safe          = .false.
    !FIXME: the code has to work with lock_outside=.true.
    lock_outside       = .false.

    !Setting transformation variables for each rank
    !**********************************************
    call mo_work_dist(v2o,fai1,tl1)
    call mo_work_dist(o2v,fai2,tl2)


    w3size = max(tl1*no,tl2*nv)
    if(nnod>1)w3size = max(w3size,2*omega2%tsize)
    call mem_alloc(w3,w3size)





    ! (-1) t [a b i k] * F [k j] =+ Omega [a b i j]
    !if(me==0) call array_convert(t2,w_o2v2,t2%nelms)
    if(.not.lock_outside)then
      call time_start_phase(PHASE_COMM, at = tw)
      call array_gather(1.0E0_realk,t2,0.0E0_realk,w_o2v2,o2v2)
      do nod=1,nnod-1
      call mo_work_dist(nv*nv*no,fri,tri,nod)
      if(me==0)then
        do i=1,no
        call dcopy(tri,w_o2v2(fri+(i-1)*no*nv*nv),1,w3(1+(i-1)*tri),1)
        enddo
        !write(*,*) 'printing w3', w3
      endif
      if(me==0.or.me==nod)then
        write(*,*) 'send recv', me
        call ls_mpisendrecv(w3(1:no*tri),int((i8*no)*tri,kind=long),infpar%lg_comm,infpar%master,nod)
        write(*,*) 'after recv',me
      endif
      enddo
      if(me==0)then
        do i=1,no
        call dcopy(tl1,w_o2v2(fai1+(i-1)*no*nv*nv),1,w3(1+(i-1)*tl1),1)
        enddo
      endif
      w_o2v2=0.0E0_realk
      call time_start_phase(PHASE_WORK, at = tc)
    else
      call arr_unlock_wins(t2)
    endif

    print *,me,"contraction done",omega2%tdim
    call lsmpi_barrier(infpar%lg_comm)

    if(.not.lock_outside)then
      call dgemm('n','n',tl1,no,no,-1.0E0_realk,w3,tl1,pfock%elm1,no,0.0E0_realk,w_o2v2(fai1),v2o)
      call time_start_phase(PHASE_COMM, at = tw)
      call lsmpi_local_reduction(w_o2v2,o2v2,infpar%master)
      call array_scatteradd_densetotiled(omega2,1.0E0_realk,w_o2v2,o2v2,infpar%master)
      call time_start_phase(PHASE_WORK, at = tc)
    endif


    !DO ALL THINGS DEPENDING ON 2


    ! F[a c] * t [c b i j] =+ Omega [a b i j]
    if(.not.lock_outside)then
      call time_start_phase(PHASE_COMM, at = tw)
      call array_gather(1.0E0_realk,t2,0.0E0_realk,w_o2v2,o2v2)
      do nod=1,nnod-1
      call mo_work_dist(nv*no*no,fri,tri,nod)
      if(me==0)then
        do i=1,tri
        call dcopy(nv,w_o2v2(1+(fri+i-2)*nv),1,w3(1+(i-1)*nv),1)
        enddo
      endif
      if(me==0.or.me==nod)then
        call ls_mpisendrecv(w3(1:nv*tri),int((i8*nv)*tri,kind=long),infpar%lg_comm,infpar%master,nod)
        call time_start_phase(PHASE_WORK, at = tc)
      endif
      enddo
      if(me==0)then
        do i=1,tl2
        call dcopy(nv,w_o2v2(1+(fai2+i-2)*nv),1,w3(1+(i-1)*nv),1)
        enddo
      endif
      w_o2v2=0.0E0_realk
      else
        call arr_unlock_wins(t2)
      endif


      if(.not.lock_outside)then
        call dgemm('n','n',nv,tl2,nv,1.0E0_realk,qfock%elm1,nv,w3,nv,0.0E0_realk,w_o2v2(1+(fai2-1)*nv),nv)
        call time_start_phase(PHASE_COMM, at = tw)
        call lsmpi_local_reduction(w_o2v2,o2v2,infpar%master)
        call array_scatteradd_densetotiled(omega2,1.0E0_realk,w_o2v2,o2v2,infpar%master)
        call time_start_phase(PHASE_WORK, at = tc)
      endif



      !INTRODUCE PERMUTATION
      omega2%access_type = MASTER_ACCESS
      t2%access_type     = MASTER_ACCESS
      pfock%access_type    = MASTER_ACCESS
      qfock%access_type    = MASTER_ACCESS

      if(.not.lock_outside)then
         call time_start_phase(PHASE_COMM, at = tw)
         call array_gather(1.0E0_realk,omega2,0.0E0_realk,w_o2v2,o2v2,wrk=w3,iwrk=w3size)
         call array_gather(1.0E0_realk,omega2,1.0E0_realk,w_o2v2,o2v2,oo=[2,1,4,3],wrk=w3,iwrk=w3size)
         call array_scatter_densetotiled(omega2,w_o2v2,o2v2,infpar%master)
         call time_start_phase(PHASE_WORK, at = tc)
      else
         if(me==0)then
            call arr_lock_wins(omega2,'s',mode)
            call array_gather(1.0E0_realk,omega2,0.0E0_realk,w_o2v2,o2v2,oo=[2,1,4,3])
            call arr_unlock_wins(omega2,.true.)
            call arr_lock_wins(omega2,'s',mode)
            call array_scatter(1.0E0_realk,w_o2v2,1.0E0_realk,omega2,o2v2)
            call arr_unlock_wins(omega2,.true.)
         endif
      endif

      call mem_dealloc(w3)
      lock_outside     = lock_safe

#endif

   endif

   call mem_dealloc(w_o2v2)





    return
  end subroutine RPA_fock_para


 !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_fock_para2(omega2,t2,pfock,qfock,no,nv)

    implicit none
    type(array), intent(inout) :: t2,omega2
    !type(array4), intent(inout) :: omega2
    integer, intent(inout) :: no,nv
    !type(array2), intent(inout) :: pfock,qfock
    type(array), intent(inout) :: pfock,qfock
    !real(realk), intent(inout) :: pfock(no,no),qfock(nv,nv)
    type(array) :: tmpt2,omegaw2
    real(realk),pointer :: tmp(:,:),w1(:),w2(:),w3(:),w4(:),omegw(:),w5(:)
    real(realk), pointer :: w_o2v2(:)
    integer, dimension(4) :: tmp_dims
    type(array) :: t_par
    integer(kind=long) :: o2v2
    integer :: no2,nv2,o2v,v2o
    integer(kind=8) :: w3size,b0
    integer :: faip,faiv,tlp,tlv,dim1,fai1,fai2,tl1,tl2,fri,tri
    integer :: i 
    real(realk) :: tw,tc
    integer(kind=ls_mpik) me,mode,nod, nnod
    character(ARR_MSG_LEN) :: msg
    logical :: master,lock_outside,lock_safe,local

    master=.true.
#ifdef VAR_MPI
    master        = .false.
    master        = (infpar%lg_mynum == infpar%master)
    nnod          = infpar%lg_nodtot
    me            = infpar%lg_mynum
    mode          = int(MPI_MODE_NOCHECK,kind=ls_mpik)
#endif

    
    dim1=no*nv
    b0= i8*0

    no2  = no**2
    nv2  = nv**2
    o2v2 = (i8*no2)*(i8*nv2)
    o2v  = no2*nv 
    v2o  = nv2*no 
    

#ifdef VAR_MPI
     StartUpSlaves: if(master .and. infpar%lg_nodtot>1) then
       call ls_mpibcast(RPAGETFOCK,infpar%master,infpar%lg_comm)
       call rpa_fock_communicate_data(t2,omega2,pfock,qfock,nv,no)
    endif StartUpSlaves
#endif


!#ifdef VAR_MPI

    !Setting transformation variables for each rank
    !**********************************************
    call mo_work_dist(v2o,fai1,tl1)
    call mo_work_dist(o2v,fai2,tl2)


    w3size = max(tl1*no,tl2*nv)
    if(nnod>1)w3size = max(w3size,2*omega2%tsize)
    call mem_alloc(w3,w3size)
    call mo_work_dist(nv*nv*no,fri,tri)
    call mem_alloc(w2,tl1*no)
    call mem_alloc(w_o2v2,tl1*no)


    t_par = array_ainit([nv,nv,no,no],4,atype='TDAR',local=.false.)
    call array_convert(t2%elm4,t_par)
    !call array_two_dim_1batch(t2,[1,2,3,4],'g',t_par%elm1,1,0,no2*nv2,.false.,debug=.true.)

    ! (-1) t [a b i k] * F [k j] =+ Omega [a b i j]
    !if(me==0) call array_convert(t2,w_o2v2,t2%nelms)
      !call array_gather(1.0E0_realk,t2,0.0E0_realk,w_o2v2,o2v2)

#ifdef VAR_MPI
      write(*,*) 'Johannes first two_dim_batch',infpar%lg_mynum,fai1,tl1
      write(*,*) 'Johannes first fri tri',infpar%lg_mynum,fri,tri
#endif

      call array_two_dim_1batch(t_par,[4,2,1,3],'g',w2,3,fai1,tl1,.false.,debug=.true.)
      write(*,*) 'Johannes after first two_dim_batch'

      call dgemm('n','n',tl1,no,no,-1.0E0_realk,w2,tl1,pfock%elm1,no,0.0E0_realk,w_o2v2,v2o)
      write(*,*) 'Johannes Debug after dgemm'

    call array_two_dim_1batch(omega2,[1,2,3,4],'a',w_o2v2,2,fri,tri,.false.,debug=.false.)
      write(*,*) 'Johannes Debug after first omega2 update'
    call array_two_dim_1batch(omega2,[3,4,1,2],'a',w_o2v2,2,fri,tri,.false.,debug=.false.)
      write(*,*) 'Johannes Debug after 2 omega2 update'
#ifdef VAR_MPI
    call lsmpi_barrier(infpar%lg_comm)
#endif

    call mem_dealloc(w_o2v2)
    call mem_dealloc(w2)
    !DO ALL THINGS DEPENDING ON 2
    call mem_alloc(w_o2v2,tl2*nv)
    call mem_alloc(w2,tl2*nv)

    ! F[a c] * t [c b i j] =+ Omega [a b i j]
      call array_two_dim_1batch(t_par,[4,3,2,1],'g',w2,3,fai2,tl2,.false.,debug=.false.)
      !call array_gather(1.0E0_realk,t2,0.0E0_realk,w_o2v2,o2v2)
      call dgemm('n','n',nv,tl2,nv,1.0E0_realk,qfock%elm1,nv,w3,nv,0.0E0_realk,w_o2v2(1+(fai2-1)*nv),nv)

    call array_two_dim_1batch(omega2,[4,3,2,1],'a',w_o2v2,2,fai2,tl2,.false.,debug=.false.)
    call array_two_dim_1batch(omega2,[3,4,1,2],'a',w_o2v2,2,fai2,tl2,.false.,debug=.false.)
#ifdef VAR_MPI
    call lsmpi_barrier(infpar%lg_comm)
#endif


!#endif


   call array_free(t_par)
   call mem_dealloc(w_o2v2)

    return
  end subroutine RPA_fock_para2




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
         1.0E0_realk,gmo,dim1,Sckdl%val,dim1,0.0E0_realk,Dckbj%val,dim1)

    call array4_reorder(omega2,[2,1,4,3])
    call dgemm('n','n',dim1,dim1,dim1, &
         2.0E0_realk,Sckdl%val,dim1,Dckbj%val,dim1,1.0E0_realk,omega2%val,dim1)

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
    !type(array), intent(inout) :: omega2,t2
    real(realk),pointer, intent(inout) :: gmo(:)
    !type(array), intent(inout) :: gmo
    type(array2), intent(inout) :: pfock,qfock
    !type(array) :: foo,fvv
    integer, intent(in) :: nocc,nvirt
    real(realk),pointer :: foo(:,:),fvv(:,:)
    !real(realk), intent(inout) :: pfock(nocc,nocc),qfock(nvirt,nvirt)
    !type(array4) :: tmp
    real(realk),pointer :: w2(:)!,foo(:,:),fvv(:,:)
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k
    real(realk) :: starttime,stoptime
    type(array) :: om_par,t_par,t_par1
    integer :: nnod,mynum,fai,tl
    integer :: nvir,noc
    logical :: master
    character(ARR_MSG_LEN) :: msg
    type(array) :: Sckdl
    

    nvir=nvirt
    noc=nocc

    !call mem_alloc(foo,noc,noc)
    !call mem_alloc(fvv,nvir,nvir)

    !foo = array_minit([nocc,nocc],2,atype='TDAR')
    !fvv = array_minit([nvirt,nvirt],2,atype='TDAR',local=.false.)
    !call array_convert(pfock%val,foo)
    !call array_convert(qfock%val,fvv)
    call mem_alloc(foo,nocc,nocc)
    call mem_alloc(fvv,nvirt,nvirt)
    foo(:,:)=pfock%val(:,:)
    fvv(:,:)=qfock%val(:,:)


    call cpu_time(starttime)

    t_par = array_minit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')

    call array4_reorder(t2,[1,3,2,4])
    call array_convert(t2%val,t_par)
    call array4_reorder(t2,[1,3,2,4])
    !write(*,*) 't2',t_par%elm1

    !write(*,*) 't_par dims',t_par%dims
    !write(*,*) 't_par tdim',t_par%tdim
    !write(*,*) 'init om_par'
    !om_par = array_minit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')
    !write(*,*) 'om_par dims',om_par%dims
    !write(*,*) 'om_par tdim',om_par%tdim

    call RPA_fock_partdeb(omega2,t2,pfock,qfock,nocc,nvirt)
    !call RPA_fock_para(omega2,t_par,foo,fvv,noc,nvir)


    !call mem_dealloc(foo)
    !call mem_dealloc(fvv)

    !call array4_reorder(omega2,[1,3,2,4])

    !write(*,*) 'tiles',om_par%ntiles
    !call array_gather(1.0E0_realk,om_par,0.0E0_realk,omega2%val,i8*nvirt*nocc*nvirt*nocc)

    !call array4_reorder(omega2,[1,3,2,4])


    Sckdl = array_minit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')
    !call copy_array(t2,Sckdl)
    call array4_reorder(t2,[1,3,2,4])
    call array_convert(t2%val,Sckdl)
    !Sckdl = array4_duplicate(t2)
    call array4_reorder(t2,[1,3,2,4])

    do a=1,nvirt
    do i=1,nocc
    Sckdl%elm4(a,a,i,i)=Sckdl%elm4(a,a,i,i)+1._realk
    !remember to change to a,a,i,i,when t2 = array
    enddo
    enddo

    !write(*,*) 'JOHANNES PAR addres'
    !call RPA_residual_addpar(omega2,Sckdl,gmo,noc,nvir)
    !write(*,*) 'JOHANNES added'

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
  subroutine RPA_residual_par(omega2,t2,gmo,pfock,qfock,nocc,nvirt)

    implicit none
    !type(array4), intent(inout) :: omega2,t2
    type(array), intent(inout) :: omega2
    type(array), intent(inout) :: t2
    integer, intent(in) :: nocc,nvirt
    real(realk),pointer, intent(inout) :: gmo(:)
    !type(array), intent(inout) :: gmo
    type(array),intent(inout) :: pfock,qfock
    !type(array2), intent(inout) :: pfock,qfock
    !real(realk),pointer,intent(inout) :: pfock(:),qfock(:)
    !real(realk), intent(inout) :: pfock(nocc,nocc),qfock(nvirt,nvirt)
    !type(array4) :: tmp
    real(realk),pointer :: w2(:)!,foo(:,:),fvv(:,:)
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k
    real(realk) :: starttime,stoptime
    type(array) :: om_par,t_par,t_par1
    integer :: nnod,mynum,fai,tl
    integer :: nvir,noc,dim1
    logical :: master
    character(ARR_MSG_LEN) :: msg
    type(array) :: Sckdl
    
    nvir=nvirt
    noc=nocc
    dim1=nocc*nvirt
  !  write(*,*) 'Johannes addr res_par',omega2%addr_p_arr
  !  write(*,*) 'Johannes addr res_par',t2%addr_p_arr

!    call mem_alloc(foo,nocc*nocc)
!    call mem_alloc(fvv,nvirt*nvirt)


    call cpu_time(starttime)

    !t_par = array_minit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')

    Sckdl = array_minit([nvirt,nvirt,nocc,nocc],4,local=.true.,atype='TDAR')
    call array_gather(1.0E0_realk,t2,0.0E0_realk,Sckdl%elm1,i8*dim1*dim1,oo=[1,2,3,4])

    !call RPA_fock_para2(omega2,Sckdl,pfock,qfock,noc,nvir)
    call RPA_fock_para(omega2,Sckdl,pfock,qfock,noc,nvir)
    msg = 'Norm of fockpart'
    call print_norm(omega2,msg)

    !msg = 'Norm of omega2 after fock para'
    !call print_norm(omega2,msg)

    !Sckdl = array_minit([nvirt,nvirt,nocc,nocc],4,local=.true.,atype='TDAR')
    !call array_gather(1.0E0_realk,t2,0.0E0_realk,Sckdl%elm1,i8*dim1*dim1,oo=[4,2,3,1])

    do a=1,nvirt
    do i=1,nocc
    Sckdl%elm4(a,a,i,i)=Sckdl%elm4(a,a,i,i)+1._realk
    enddo
    enddo
    !call RPA_residual_addpar(omega2,Sckdl,gmo,noc,nvir)
    call RPA_residual_par_add(omega2,Sckdl,gmo,noc,nvir)
    !write(*,*) 'JOHANNES added'

    call cpu_time(stoptime)


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
!!    t_par = array_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')
!!
!    call array4_reorder(t2,[1,3,2,4])
!!    call array_convert(t2%val,t_par)
!
!!    call array_two_dim_1batch(t_par,[1,3,2,4],'g',w2,2,fai,tl,.false.,debug=.true.)
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
  subroutine RPA_residual_par_add(omega2,u2,gmo,nocc,nvirt)

    implicit none
    !type(array4), intent(inout) :: omega2!,u2
    type(array), intent(inout) :: u2 ,omega2
    real(realk), intent(inout),pointer :: gmo(:)
    !type(array), intent(inout):: gmo
    integer,intent(inout) :: nocc,nvirt
    !type(array4) :: Sckdl,Dckbj
    type(array) :: Sckdl,Dckbj
    type(array) :: t_par,omegaw1
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k,dim1
    integer :: fai,tl,mynum,nnod
    real(realk) :: starttime,stoptime
    real(realk),pointer :: w2(:),w3(:),omegw(:),w4(:)
    character(ARR_MSG_LEN) :: msg
    logical :: master

    master=.true.
#ifdef VAR_MPI
    master        = .false.
    mynum         = infpar%lg_mynum
    master        = (infpar%lg_mynum == infpar%master)
    nnod          = infpar%lg_nodtot
#endif


    dim1=nocc*nvirt
#ifdef VAR_MPI
    StartUpSlaves: if(master .and. infpar%lg_nodtot>1) then
      call ls_mpibcast(RPAGETRESIDUAL,infpar%master,infpar%lg_comm)
      call rpa_res_communicate_data(gmo,u2,omega2,nvirt,nocc)
    endif StartUpSlaves
#endif


    call mo_work_dist(nvirt*nocc,fai,tl)
    !call mem_alloc(w2,tl*nocc*nvirt)
!#ifdef VAR_MPI
!    write(*,*) 'size of tile',tl,fai,infpar%lg_mynum,nvirt*nocc!*nvirt*nocc
!#endif
    call mem_alloc(w2,tl*nocc*nvirt)
    call mem_alloc(w3,tl*nocc*nvirt)
    call mem_alloc(w4,nocc*nvirt*nocc*nvirt)
    call mem_alloc(omegw,tl*nocc*nvirt)
    !call mem_alloc(omegw,tl*tl)

    t_par = array_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR',local=.false.)
    call array_convert(u2%elm4,t_par)

    call array_two_dim_1batch(t_par,[4,2,3,1],'g',w2,2,fai,tl,.false.,debug=.true.)


    omegw=0.0_realk

    !When fock part is parallelized instead of zero 1.0_realk
    call dgemm('n','n',tl,dim1,dim1, &
         1.0E0_realk,w2,tl,gmo,dim1,0.0E0_realk,w3,tl)
!#ifdef VAR_MPI
!    write(msg,*) 'Norm of w3',infpar%lg_mynum
!#else
!    write(msg,*) 'Norm of w3'
!#endif
!    call print_norm(w3,i8*tl*nvirt*nocc,msg)

!#ifdef VAR_MPI
    !call array_gather(1.0E0_realk,u2,0.0E0_realk,w4,i8*dim1*dim1,oo=[1,3,2,4])
    !call array_gather(1.0E0_realk,t_par,0.0E0_realk,w4,i8*dim1*dim1,oo=[1,3,2,4])
    call array_gather(1.0E0_realk,t_par,0.0E0_realk,w4,i8*dim1*dim1,oo=[4,2,3,1])
!    write(msg,*) 'Norm of w4',infpar%lg_mynum
!#else
!    write(msg,*) 'Norm of w4'

!#endif
    !call sleep(1)
    !call lsmpi_barrier(infpar%lg_comm)
!    call print_norm(w4,i8*dim1*dim1,msg)
    !call sleep(1)
    !call lsmpi_barrier(infpar%lg_comm)
    !stop
    
    call dgemm('n','n',tl,dim1,dim1, &
         2.0E0_realk,w3,tl,w4,dim1,0.0E0_realk,omegw,tl)

!#ifdef VAR_MPI
!    write(msg,*) 'Norm of omegw',infpar%lg_mynum
    !call sleep(1)
    !call lsmpi_barrier(infpar%lg_comm)
!    call print_norm(omegw,i8*tl*nvirt*nocc,msg)
    !call sleep(1)
    !call lsmpi_barrier(infpar%lg_comm)
    !stop
!#endif
    

    !call array_two_dim_1batch(omegaw1,[1,3,2,4],'a',omegw,2,fai,tl,.false.,debug=.true.)
    call array_two_dim_1batch(omega2,[4,2,3,1],'a',omegw,2,fai,tl,.false.,debug=.true.)

#ifdef VAR_MPI

    call lsmpi_barrier(infpar%lg_comm)

#endif

!#ifdef VAR_MPI
!    if(master) then

!      call array_gather(1.0E0_realk,omegaw1,0.0E0_realk,omega2%val,i8*nvirt*nocc*nvirt*nocc)
      
!      call array4_reorder(omega2,[1,3,2,4])

     ! write(msg,*) 'Norm of omega2'
     ! call print_norm(omega2%val,i8*dim1*dim1,msg)
!      write(*,*) 'checkpoint 3',infpar%lg_mynum
      !call array4_reorder(omega2,[2,1,4,3])
!    endif


    call array_free(t_par)
    call array_free(u2)
    !call array_free(omegaw1)
    call mem_dealloc(w2)
    call mem_dealloc(w3)
    call mem_dealloc(w4)
    call mem_dealloc(omegw)


  end subroutine RPA_residual_par_add


  !\brief Calculate additional linear and quadratic terms of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
!  subroutine RPA_residual_addpar(omega2,u2,gmo,nocc,nvirt)
!
!    implicit none
!    type(array4), intent(inout) :: omega2!,u2
!    type(array), intent(inout) :: u2 !,omega2
!    real(realk), intent(inout),pointer :: gmo(:)
!    !type(array), intent(inout):: gmo
!    integer,intent(inout) :: nocc,nvirt
!    !type(array4) :: Sckdl,Dckbj
!    type(array) :: Sckdl,Dckbj
!    type(array) :: t_par,omegaw1
!    integer, dimension(4) :: tmp_dims
!    integer :: a,b,c,i,j,k,dim1
!    integer :: fai,tl,mynum,nnod
!    real(realk) :: starttime,stoptime
!    real(realk),pointer :: w2(:),w3(:),omegw(:),w4(:)
!    character(ARR_MSG_LEN) :: msg
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
!    omegaw1 = array_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR',local=.false.)
!    !write(*,*) 'omegaw1 tdim=',omegaw1%tdim
!    !write(*,*) 'omegaw1 dims=',omegaw1%dims
!    !write(*,*) 'omegaw1 tile=',omegaw1%ntiles
!
!!#ifdef VAR_MPI
!!    if(master) then
!!      call array4_reorder(omega2,[1,3,2,4])
!!      call array_convert(omega2%val,omegaw1)
!!      !write(msg,*) 'Norm of omegaw1',infpar%lg_mynum
!!      write(*,*) 'omegaw1 tdim=',omegaw1%tdim
!!      !call print_norm(omegaw1,msg)
!!    endif
!!#else
!!#endif
!
!    call array4_reorder(omega2,[1,3,2,4])
!    call array_convert(omega2%val,omegaw1)
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
!    t_par = array_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR',local=.false.)
!    !call array4_reorder(u2,[1,3,2,4])
!    !call copy_array(u2,t_par)
!    call array_convert(u2%elm4,t_par)
!
!    call array_two_dim_1batch(t_par,[4,2,3,1],'g',w2,2,fai,tl,.false.,debug=.true.)
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
!    !call array_gather(1.0E0_realk,u2,0.0E0_realk,w4,i8*dim1*dim1,oo=[1,3,2,4])
!    !call array_gather(1.0E0_realk,t_par,0.0E0_realk,w4,i8*dim1*dim1,oo=[1,3,2,4])
!    call array_gather(1.0E0_realk,t_par,0.0E0_realk,w4,i8*dim1*dim1,oo=[4,2,3,1])
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
!    !call array_two_dim_1batch(omegaw1,[1,3,2,4],'a',omegw,2,fai,tl,.false.,debug=.true.)
!    call array_two_dim_1batch(omegaw1,[4,2,3,1],'a',omegw,2,fai,tl,.false.,debug=.true.)
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
!      call array_gather(1.0E0_realk,omegaw1,0.0E0_realk,omega2%val,i8*nvirt*nocc*nvirt*nocc)
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
!    call array_free(t_par)
!    call array_free(omegaw1)
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
         1.0E0_realk,gmo%val,dim1,t2%val,dim1,0.0E0_realk,omega2%val,dim1)

    call dgemm('n','n',dim1,dim1,dim1, &
         1.0E0_realk,t2%val,dim1,gmo%val,dim1,1.0E0_realk,omega2%val,dim1)
    
       !gmo*t2*m2
    call dgemm('n','n',dim1,dim1,dim1, &
         1.0E0_realk,gmo%val,dim1,t2_final%val,dim1,0.0E0_realk,Sckdl%val,dim1)

    call dgemm('n','n',dim1,dim1,dim1, &
         1.0E0_realk,Sckdl%val,dim1,t2%val,dim1,1.0E0_realk,omega2%val,dim1)

       !m2*t2*gmo
    call dgemm('n','n',dim1,dim1,dim1, &
         1.0E0_realk,t2%val,dim1,t2_final%val,dim1,0.0E0_realk,Sckdl%val,dim1)

    call dgemm('n','n',dim1,dim1,dim1, &
         1.0E0_realk,Sckdl%val,dim1,gmo%val,dim1,1.0E0_realk,omega2%val,dim1)

    call array4_free(Sckdl)

    return

  end subroutine RPA_multi_add
 

  function get_rpa_energy_arrnew(t2,gmo,nocc,nvirt) result(ecorr)

    implicit none
    type(array), intent(in) :: t2
    type(array), intent(inout) :: gmo
    integer, intent(in) :: nocc,nvirt
    real(realk) :: ecorr,ecorr_d
    integer :: a,i,b,j

    ecorr = 0.0E0_realk
    ecorr_d = 0.0E0_realk

    if(t2%itype==DENSE.and.gmo%itype==DENSE)then
      !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) PRIVATE(i,a,j,b) SHARED(nocc,&
      !$OMP nvirt,t2,gmo) REDUCTION(+:ecorr_d)
      do j=1,nocc
      do b=1,nvirt
      do i=1,nocc
      do a=1,nvirt
      ecorr_d = ecorr_d + t2%elm4(a,b,i,j)* &
        (1.0E0_realk*gmo%elm4(i,a,j,b))
      end do
      end do
      end do
      end do
      !$OMP END PARALLEL DO

      ecorr = ecorr_d

    else if(t2%itype==TILED_DIST.and.gmo%itype==TILED_DIST)then

      ecorr=get_rpa_energy_parallel(t2,gmo)

    endif


  end function get_rpa_energy_arrnew


  function get_sosex_cont_arrnew(t2,gmo,nocc,nvirt) result(ecorr)

    implicit none
    type(array), intent(in) :: t2
    type(array), intent(inout) :: gmo
    integer, intent(in) :: nocc,nvirt
    real(realk) :: ecorr,ecorr_d
    integer :: a,i,b,j

    ecorr = 0.0E0_realk
    ecorr_d = 0.0E0_realk

    if(t2%itype==DENSE.and.gmo%itype==DENSE)then
      !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) PRIVATE(i,a,j,b) SHARED(nocc,&
      !$OMP nvirt,t2,gmo) REDUCTION(+:ecorr_d)
      do j=1,nocc
      do b=1,nvirt
      do i=1,nocc
      do a=1,nvirt
      ecorr_d = ecorr_d + t2%elm4(a,b,i,j)* &
        (-0.5E0_realk*gmo%elm4(i,b,j,a))
      end do
      end do
      end do
      end do
      !$OMP END PARALLEL DO

      ecorr = ecorr_d

    else if(t2%itype==TILED_DIST.and.gmo%itype==TILED_DIST)then

      ecorr=get_sosex_cont_parallel(t2,gmo)

    endif


  end function get_sosex_cont_arrnew


  function get_rpa_energy_parallel(t2,gmo) result(Ec)
    implicit none
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
      call pdm_array_sync(infpar%lg_comm,JOB_GET_MP2_ENERGY,t2,gmo)
    endif
    call memory_allocate_array_dense(gmo)
    call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)

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
      !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(gmo,o,t,&
      !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) REDUCTION(+:E1,E2) COLLAPSE(3)
      do j=1,dj
        do i=1,di
          do b=1,db
            do a=1,da
     
              E2 = E2 + t(a,b,i,j)*&
              & (1.0E0_realk*  gmo%elm4(i+o(3),a+o(1),j+o(4),b+o(2)))
   
            enddo 
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      nullify(t)
    enddo

    call arr_deallocate_dense(gmo)
    
    call lsmpi_local_reduction(E2,infpar%master)

    Ec = E2
#else
    Ec = 0.0E0_realk
#endif
  end function get_rpa_energy_parallel


  function get_sosex_cont_parallel(t2,gmo) result(Ec)
    implicit none
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
      call pdm_array_sync(infpar%lg_comm,JOB_GET_MP2_ENERGY,t2,gmo)
    endif
    call memory_allocate_array_dense(gmo)
    call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)

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
      !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(gmo,o,t,&
      !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) REDUCTION(+:E1,E2) COLLAPSE(3)
      do j=1,dj
        do i=1,di
          do b=1,db
            do a=1,da
     
              E2 = E2 + t(a,b,i,j)*&
              & (-0.5E0_realk* gmo%elm4(i+o(3),b+o(2),j+o(4),a+o(1)) )
   
            enddo 
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      nullify(t)
    enddo

    call arr_deallocate_dense(gmo)
    
    call lsmpi_local_reduction(E2,infpar%master)

    Ec = E2
#else
    Ec = 0.0E0_realk
#endif
  end function get_sosex_cont_parallel



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
subroutine rpa_res_slave()
  use dec_typedef_module
  use typedeftype,only:lsitem
  use tensor_interface_module
  use decmpi_module,only:rpa_res_communicate_data
  use infpar_module
  use rpa_module
  implicit none
  !> number of orbitals:
  type(array) :: omega2,t2
  !type(array4) :: omega2!,t2
  real(realk),pointer :: gmo(:)!,t2(:,:,:,:)
  type(array2)  :: pfock,qfock
  integer :: nbas, nocc, nvirt
  !> how to pack integrals:
  
  print*, infpar%lg_mynum,'rpa_res_slave'
  call rpa_res_communicate_data(gmo,t2,omega2,nvirt,nocc)
  call RPA_residual_par_add(omega2,t2,gmo,nocc,nvirt)

end subroutine rpa_res_slave

subroutine rpa_fock_slave()
  use dec_typedef_module
  use typedeftype,only:lsitem
  use tensor_interface_module
  use decmpi_module,only:rpa_fock_communicate_data
  use infpar_module
  use rpa_module
  implicit none
  !> number of orbitals:
  type(array) :: t2
  type(array) :: omega2!,t2
  type(array)  :: pfock,qfock
  real(realk),pointer :: gmo(:)
  !real(realk),pointer  :: pfock(:),qfock(:)
  integer :: nbas, nocc, nvirt
  !> how to pack integrals:
  print*, infpar%lg_mynum,'rpa_fock_slave'
  call rpa_fock_communicate_data(t2,omega2,pfock,qfock,nocc,nvirt)
  write(*,*) 'slaves, calling fock_para'
  call RPA_fock_para2(omega2,t2,pfock,qfock,nocc,nvirt)
  !call RPA_fock_para(omega2,t2,pfock,qfock,nocc,nvirt)

end subroutine rpa_fock_slave
#endif
