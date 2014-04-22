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
              &rpa_residualdeb,rpa_residualpar,rpa_residual_addpar
    

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


    call cpu_time(starttime)

    !get the MP2 part of the residual
    !call getDoublesResidualMP2_simple(Omega2,t2,gmo,pfock,qfock, &
    !           & nocc,nvirt)
    write(*,*) 'I am now in residualdeb, and everything is ok'

    call RPA_fock_partdeb(omega2,t2,pfock,qfock,nocc,nvirt)
    call RPA_residual_adddeb(omega2,t2,gmo,nocc,nvirt)

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
    integer, intent(in) :: nocc,nvirt
    type(array4) :: tmp
    real(realk),pointer :: w2(:)
    integer, dimension(4) :: tmp_dims
    integer :: a,b,c,i,j,k
    real(realk) :: starttime,stoptime
    type(array) :: g_par,t_par,t_par1
    integer :: nnod,mynum,fai,tl
    integer :: nvir,noc
    logical :: master
    character(ARR_MSG_LEN) :: msg
    type(array) :: Sckdl
    

   nvir=nvirt
   noc=nocc

    call cpu_time(starttime)

    call RPA_fock_partdeb(omega2,t2,pfock,qfock,nocc,nvirt)
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

    write(*,*) 'JOHANNES PAR addres'
    call RPA_residual_addpar(omega2,Sckdl,gmo,noc,nvir)
    write(*,*) 'JOHANNES added'

    ! MPI: here you should start the slaves!!

    !call lsmpi_barrier(infpar%lg_comm)
    !call lsmpi_barrier(infpar%lg_comm)
    !print*, infpar%lg_mynum, 'par residual'
    !call lsmpi_barrier(infpar%lg_comm)
    




    call cpu_time(stoptime)

    !call lsquit('RPA_residual: Needs implementation',-1)

  end subroutine RPA_residualpar

  !\brief Calculate fock matrix part of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_fock_partpar(omega2,t2,pfock,qfock,nocc,nvirt)

!    implicit none
!    type(array4), intent(inout) :: omega2,t2
!    type(array2), intent(inout) :: pfock,qfock
!    integer, intent(in) :: nocc,nvirt
!    type(array4) :: tmp
!    integer, dimension(4) :: tmp_dims
!    integer :: a,b,c,i,j,k,noc,nvir
!
!
!#ifdef VAR_MPI
!    StartUpSlaves: if(master .and. infpar%lg_nodtot>1) then
!      call ls_mpibcast(RPAGETRESIDUAL,infpar%master,infpar%lg_comm)
!      call rpa_fock_communicate_data(t2,omega2,pfock,qfock,nvir,noc)
!    endif StartUpSlaves
!#endif
!
!    call mo_work_dist(nvirt*nocc,fai,tl)
!    call mem_alloc(w2,tl*nocc*nvirt*nocc)
!
!    t_par = array_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')
!
!    call array4_reorder(t2,[1,3,2,4])
!    call array_convert(t2%val,t_par)
!
!    call array_two_dim_1batch(t_par,[1,3,2,4],'g',w2,2,fai,tl,.false.,debug=.true.)
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
  end subroutine RPA_fock_partpar

  !\brief Calculate additional linear and quadratic terms of the RPA residual 
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  subroutine RPA_residual_addpar(omega2,u2,gmo,nocc,nvirt)

    implicit none
    type(array4), intent(inout) :: omega2!,u2
    type(array), intent(inout) :: u2 !,omega2
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
    master        = (infpar%lg_mynum == 0)
    nnod          = infpar%lg_nodtot
#endif


    dim1=nocc*nvirt
#ifdef VAR_MPI
    StartUpSlaves: if(master .and. infpar%lg_nodtot>1) then
      call ls_mpibcast(RPAGETRESIDUAL,infpar%master,infpar%lg_comm)
      call rpa_res_communicate_data(gmo,u2,omega2,nvirt,nocc)
    endif StartUpSlaves
#endif
   

    omegaw1 = array_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR',local=.false.)

!#ifdef VAR_MPI
!    if(master) then
!      call array4_reorder(omega2,[1,3,2,4])
!      call array_convert(omega2%val,omegaw1)
!      !write(msg,*) 'Norm of omegaw1',infpar%lg_mynum
!      write(*,*) 'itype omegaw1=',omegaw1%itype
!      !call print_norm(omegaw1,msg)
!    endif
!#else

    call array4_reorder(omega2,[1,3,2,4])
    call array_convert(omega2%val,omegaw1)
!#endif


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
    !call array4_reorder(u2,[1,3,2,4])
    !call copy_array(u2,t_par)
    call array_convert(u2%elm4,t_par)

!#ifdef VAR_MPI
!    write(msg,*) 'Norm of t_par',infpar%lg_mynum
!    call print_norm(t_par,msg)
!     write(*,*) 'itype t_par=',t_par%itype
!#endif

    !write(*,*) 'in residue dim t_par', t_par%tdim

    
!#ifdef VAR_MPI
    !In this u2 is array4, hence I use t_par
!    write(*,*) 'before array two_dim',infpar%lg_mynum
!    call lsmpi_barrier(infpar%lg_comm)

    !call array_two_dim_1batch(t_par,[1,3,2,4],'g',w2,2,fai,tl,.false.,debug=.true.)
    call array_two_dim_1batch(t_par,[4,2,3,1],'g',w2,2,fai,tl,.false.,debug=.true.)

    !In this u2 is array and no need for t_par
    !call array_two_dim_1batch(u2,[1,3,2,4],'g',w2,2,fai,tl,.false.,debug=.true.)

!    if(master) then
!      write(msg,*) 'Norm of t_par',infpar%lg_mynum
!      call print_norm(t_par,msg)
!    endif
!#endif
!#ifdef VAR_MPI
!    write(msg,*) 'Norm of w2',infpar%lg_mynum
!#else
!    write(msg,*) 'Norm of w2'
!#endif
!    call print_norm(w2,i8*dim1*tl,msg)
!#ifdef VAR_MPI
!    call lsmpi_barrier(infpar%lg_comm)
!#endif

 !   write(*,*) 'checkpoint 1',infpar%lg_mynum
  !  call sleep(1)
  !  call lsmpi_barrier(infpar%lg_comm)
     
    !call print_norm(gmo,i8*nocc*nvirt*nocc*nvirt)
    !
    !call print_norm(u2)
    !call lsmpi_barrier(infpar%lg_comm)
    !write(*,*) 'checkpoint 1',infpar%lg_mynum
    !call sleep(1)
    !call lsmpi_barrier(infpar%lg_comm)


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
    call array_two_dim_1batch(omegaw1,[4,2,3,1],'a',omegw,2,fai,tl,.false.,debug=.true.)
#ifdef VAR_MPI

    call lsmpi_barrier(infpar%lg_comm)
    !write(msg,*) 'Norm of omegaw1',infpar%lg_mynum
   ! if(master) write(*,*) 'omegaw1'
   ! if(master) write(*,*) omegaw1%elm4
    !call sleep(1)
    !call lsmpi_barrier(infpar%lg_comm)
    !call print_norm(omegaw1,msg)
    !call sleep(1)
    !call lsmpi_barrier(infpar%lg_comm)
    !stop
#endif

!#ifdef VAR_MPI
    if(master) then
      !call array_convert(omegaw1,omega2%val)
      !call array4_reorder(omega2,[1,3,2,4])
      !call array4_reorder(omega2,[1,3,2,4])
      !write(*,*) 'order of omega2',size(omega2%val(:,1,1,1))
      !write(*,*) 'order of omega2',size(omega2%val(1,:,1,1))
      !write(*,*) 'order of omega2',size(omega2%val(1,1,:,1))
      write(*,*) 'Master writes this'
      !call array4_reorder(omega2,[2,1,4,3])
      call array_gather(1.0E0_realk,omegaw1,0.0E0_realk,omega2%val,i8*nvirt*nocc*nvirt*nocc)

     ! write(msg,*) 'Norm of omega2'
     ! call print_norm(omega2%val,i8*dim1*dim1,msg)
!      write(*,*) 'checkpoint 3',infpar%lg_mynum
      !call array4_reorder(omega2,[2,1,4,3])
      call array4_reorder(omega2,[1,3,2,4])
    endif
!#else
!
!    call array_gather(1.0E0_realk,omegaw1,0.0E0_realk,omega2%val,i8*nvirt*nocc*nvirt*nocc)
!
!    call array4_reorder(omega2,[1,3,2,4])
!
!#endif


    call array_free(t_par)
    call array_free(omegaw1)
    call mem_dealloc(w2)
    call mem_dealloc(w3)
    call mem_dealloc(w4)
    call mem_dealloc(omegw)


  end subroutine RPA_residual_addpar


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



  !\brief Calculate RPA energy
  !> \author Johannes Rekkedal and Thomas Bondo
  !> \date March 2013
  function RPA_energy(t2,gmo) result(energy)
    implicit none
    type(array4) :: J,X,gmo
    type(array4), intent(in) :: t2
    real(realk) :: energy

    write(*,*) 'In rpa_energy'
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
  type(array) :: t2
  type(array4) :: omega2!,t2
  real(realk),pointer :: gmo(:)
  type(array2)  :: pfock,qfock
  integer :: nbas, nocc, nvirt
  !> how to pack integrals:
  
  print*, infpar%lg_mynum,'rpa_res_slave'
  call rpa_res_communicate_data(gmo,t2,omega2,nvirt,nocc)
  call RPA_residual_addpar(omega2,t2,gmo,nocc,nvirt)

end subroutine rpa_res_slave
#endif
