!> @file
!> Two-electron integrals and amplitudes for RI-MP2.
!> \ author Thomas Kjaergaard

module rimp2_module

#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
  use precision
  use lstiming!, only: lstimer
  use lowdin_module
  use screen_mod!, only: DECscreenITEM
  use dec_typedef_module
  use typedeftype!, only: Lsitem,lssetting
  use BUILDAOBATCH!,only:build_batchesofaos,determine_maxbatchorbitalsize,&
 !      & determine_MaxOrbitals
  use typedef!, only: typedef_free_setting,copy_setting
  use memory_handling
  use screen_mod!,only: free_decscreen, DECSCREENITEM
  use lsparameters
  use IntegralInterfaceMod!, only: II_getBatchOrbitalInfo
  use IntegralInterfaceDEC!, only: II_precalc_DECScreenMat,&
!       & II_getBatchOrbitalScreen, II_GET_DECPACKED4CENTER_J_ERI
  use IntegralInterfaceModuleDF
  use IchorErimoduleHost
  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
  use cc_tools_module
#ifdef VAR_MPI
      use decmpi_module !, only: mpi_communicate_mp2_int_and_amp
#endif
  use dec_workarounds_module

  use dec_fragment_utils!,only: calculate_fragment_memory, &
!       & dec_simple_dgemm_update,start_flop_counter,&
!       & end_flop_counter, dec_simple_dgemm, mypointer_init, &
!       & get_currently_available_memory, atomic_fragment_free
  use array2_simple_operations!, only: array2_free, array2_extract_EOS, &
!       & get_mp2_integral_transformation_matrices, get_mp2_integral_transformation_matrices_fc, &
 !      & extract_occupied_eos_mo_indices, extract_virtual_EOS_MO_indices,array2_init,array2_print
  use array4_simple_operations!, only: array4_delete_file, array4_init_file, &
!       & array4_init_standard, array4_free, array4_reorder, array4_init, &
!       & array4_contract1, array4_open_file, array4_write_file_type2, &
!       & array4_close_file, array4_write_file_type1, mat_transpose, &
 !     & array4_read_file_type2
  use ri_util_module
  use iso_c_binding
#ifdef VAR_OPENACC
  use openacc
#endif
  use gpu_interfaces

contains
!> Purpose: Wrapper routine to get RI-MP2 amplitudes and integrals
!           in comnnection with DECNP
!> Author:  Pablo Baudin
!> Date:    April 2015
subroutine decnp_RIMP2_integrals_and_amplitudes(frag,VOVO,t2)
   implicit none
  !> Atomic fragment
  type(decfrag), intent(inout) :: frag
  !> full AOS integrals
  type(tensor),intent(inout) :: VOVO
  !> full AOS amplitudes
  type(tensor),intent(inout) :: t2
  !> dummy tensors should not be allocated in RIMP2 routine
  type(tensor) :: tdum, gdum
  !integer :: dims(4)

  !dims = [frag%nvirtAOS,frag%noccAOS,frag%nvirtAOS,frag%noccAOS]   ! Output order
  !call tensor_ainit(VOVO,dims,4)
  !call tensor_ainit(t2,dims,4)
  !call tensor_zero(VOVO)
  !call tensor_zero(t2)

  call RIMP2_integrals_and_amplitudes(frag,VOVO,t2,tdum,gdum)

end subroutine decnp_RIMP2_integrals_and_amplitudes

!> \brief Calculate EOS integrals and EOS amplitudes for RI-MP2 calculation -
!> both for occupied and virtual partitioning schemes.
!> \author Thomas Kjaergaard
!> \date Marts 2014
subroutine RIMP2_integrals_and_amplitudes(MyFragment,&
     & goccEOS, toccEOS,gvirtEOS, tvirtEOS, djik,blad)

  implicit none
  !FIXME : Memory Usage RIMP2MEM, 
  !FIXME : nAux distribution better - for better load balancing - split atoms

  !> Atomic fragment (or pair fragment)
  type(decfrag), intent(inout) :: MyFragment
  !> Integrals for occ EOS: (d j|c i) in the order (d,j,c,i) 
  type(tensor),intent(inout) :: goccEOS
  !> Amplitudes for occ EOS in the order (d,j,c,i) 
  type(tensor),intent(inout) :: toccEOS
  !> Integrals for virt EOS: (b l|a k) in the order (b,l,a,k) 
  type(tensor),intent(inout) :: gvirtEOS
  !> Amplitudes for virt EOS in the order (b,l,a,k) 
  type(tensor),intent(inout) :: tvirtEOS
  !> Occ EOS integrals (d j | i k) in the order (d,j,i,k)  
  type(tensor),optional :: djik  !V_aos*O_eos*O_eos*O_aos
  !> Virt EOS integrals (b l | a d) in the order (b,l,a,d)  
  type(tensor),optional :: blad  !V_eos*O_aos*V_eos*V_aos
  !local variables
  type(mp2_batch_construction) :: bat
  type(array2) :: CDIAGocc, CDIAGvirt, Uocc, Uvirt
  type(array2) :: LoccEOS,LvirtEOS, tmparray2, LoccTALL,CDIAGoccTALL,UoccALL
  real(realk), pointer :: EVocc(:), EVvirt(:)
  integer :: nbasis,nocc,nvirt, noccEOS, nvirtEOS,nocctot,ncore,nCoccTmp
  integer :: alpha,gamma,beta,delta,info,mynum,numnodes,nb
  integer :: IDIAG,JDIAG,ADIAG,BDIAG,ALPHAAUX,myload,nb2,natomsAux
  integer :: ILOC,JLOC,ALOC,BLOC,M,N,K,nAtoms,nbasis2,nbasisAux
  logical :: fc,ForcePrint,master,wakeslave
  logical :: CollaborateWithSlaves
  real(realk),pointer :: Calpha(:),Calpha2(:),Calpha3(:),CalphaVV(:)
  real(realk),pointer :: UoccEOS(:,:),UvirtEOS(:,:),Ctmp(:),Ctmp2(:)
  real(realk),pointer :: tocc(:),UoccEOST(:,:),UvirtT(:,:),tocc3(:)
  real(realk),pointer :: toccTMP(:,:),TMPAlphaBeta_minus_sqrt(:,:),tocc2(:)
  real(realk),pointer :: tvirtTMP(:,:),tvirt(:),UoccT(:,:),UvirtEOST(:,:)
  real(realk),pointer :: tvirt2(:),tvirt3(:),CoccEOS(:,:)
  real(realk),pointer :: UoccallT(:,:),CalphaOcc(:),tocc2TMP(:)
  real(realk) :: deltaEPS,goccAIBJ,goccBIAJ,Gtmp,Ttmp,Eocc,TMP,Etmp,twmpi2
  real(realk) :: gmocont,Gtmp1,Gtmp2,Eocc2,TMP1,flops,tmpidiff,EnergyMPI(2)
  real(realk) :: tcpu, twall,tcpu1,twall1,tcpu2,twall2,tcmpi1,tcmpi2,twmpi1
  real(realk) :: Evirt,Evirt2,dummy(2),MemInGBCollected,gpuflops
  real(realk) :: maxsize
  Integer :: iAtomA,nBastLocA,startRegA,endRegA,nAuxA,startAuxA,endAuxA,lupri
  integer :: MynAtomsMPI,startA2,StartA,B,I,startB2,iAtomB,StartB,node,myOriginalRank
  Integer :: OriginalRanknbasisAuxMPI,NBA,dimocc(4),dimvirt(4)
  real(realk) :: time_i,time_c,time_w
  real(realk),pointer :: OccContribsFull(:),VirtContribsFull(:),Calpha_debug(:,:,:)
  real(realk),pointer :: occ_tmp(:),virt_tmp(:),ABdecomp(:,:),CDIAGoccALL(:,:)
  real(realk),pointer :: CvirtAOS(:,:),CvirtEOS(:,:),CoccTmp(:,:),CDIAGoccALLcf(:,:)
  real(realk),pointer :: TauOcc(:,:),TauVirt(:,:)
  logical :: ABdecompCreate
  integer,pointer :: IPVT(:)
  integer,pointer :: nbasisAuxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
  TYPE(MOLECULARORBITALINFO) :: orbitalInfo
  real(realk), pointer   :: work1(:),Etmp2222(:)
  real(realk)            :: RCOND
  integer(kind=ls_mpik)  :: COUNT,TAG,IERR,request,Receiver,sender,J,COUNT2
  real(realk) :: TS,TE,TS2,TE2,TS3,TE3,CPUTIME,WALLTIMESTART,WALLTIMEEND,WTIME
  real(realk) :: tcpu_start,twall_start, tcpu_end,twall_end,MemEstimate
  real(realk) :: MemStep1,MemStep2,MemStep3,TS4,TE4
  integer ::CurrentWait(2),nAwaitDealloc,iAwaitDealloc,oldAORegular,oldAOdfAux
  integer :: MaxVirtSize,nTiles,offsetV,offset,MinAuxBatch
  integer :: noccOut,nvirtOut
  logical :: useAlphaCD5,useAlphaCD6,ChangedDefault,first_order,PerformTiling
  logical :: use_bg_buf
  integer(kind=ls_mpik)  :: request5,request6
  real(realk) :: phase_cntrs(nphases),bytes_to_alloc,MinMem
  integer(kind=long) :: nSize,nsize1,nsize2,nsize3,nbasisAux8
  integer(kind=long) :: nocc8,nvirt8,noccEOS8,nbasis8,nvirtEOS8,nocctot8
  character :: intspec(5)
  TYPE(MoleculeInfo),pointer :: molecule1,molecule2,molecule3,molecule4
  !Laplace values
  integer,parameter :: nLaplace=10
  real(realk),parameter,dimension(10) :: LaplaceAmp = (/ -0.003431, &
       & -0.023534, -0.088984, -0.275603, -0.757121, -1.906218, -4.485611, &
       & -10.008000, -21.491075, -45.877205 /)
  real(realk),parameter,dimension(10) :: LaplaceW = (/ 0.009348, &
       & 0.035196, 0.107559, 0.293035, 0.729094, 1.690608, 3.709278, &
       & 7.810243, 16.172017, 35.929402 /)
  
  ! cublas stuff
  type(c_ptr) :: cublas_handle
  integer*4 :: stat
  !> async handles
  integer :: num_ids
#ifdef VAR_OPENACC
  integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
#ifdef VAR_PGF90
  integer*4, external :: acc_set_cuda_stream
#endif
  integer(c_size_t) :: total_gpu,free_gpu ! total and free gpu mem in bytes
#else
  integer, pointer, dimension(:) :: async_id
#endif
#ifdef VAR_PAPI
  integer(8) :: papiflops
  integer :: eventset2
#endif
#ifdef VAR_MPI
  INTEGER(kind=ls_mpik) :: HSTATUS
  CHARACTER*(MPI_MAX_PROCESSOR_NAME) ::  HNAME
  TAG = 1319
#endif  

#ifdef VAR_PAPI
  CALL LS_GETTIM(CPUTIME,WALLTIMESTART)
  call myPAPI_start(eventset2)
#endif
  ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
  num_ids = 7
  call mem_alloc(async_id,num_ids)

#ifdef VAR_OPENACC

  if (DECinfo%acc_sync) then
     async_id = acc_async_sync
  else
     do m = 1,num_ids
        async_id(m) = int(m,kind=acc_handle_kind)
     enddo
  endif

#else

  if (DECinfo%acc_sync) then
     async_id = 0
  else
     do m = 1,num_ids
        async_id(m) = -m
     enddo
  endif

#endif

#ifdef VAR_CUBLAS

  ! initialize the CUBLAS context
  stat = cublasCreate_v2(cublas_handle)
  print*,'cublasCreate_v2 gives stat=',stat
  stat = cudaSetDevice(int(0,kind=4))
  print*,'cudaSetDevice(0) gives stat=',stat
  print*,'known results are:'
  !print*,'cudaSuccess',cudaSuccess
  !print*,'cudaErrorInvalidDevice',cudaErrorInvalidDevice
  !print*,'cudaErrorDeviceAlreadyInUse ',cudaErrorDeviceAlreadyInUse 

  ! set the cublas handle to match the synchronous openacc handle 
  stat = acc_set_cuda_stream(acc_async_sync,cublas_handle)

#endif

  IF(present(djik))THEN
     IF(present(blad))THEN
        first_order=.TRUE. !first order integrals are required
     ELSE
        call lsquit('RIMP2 prog error: djik is present but blad is not',-1)
     ENDIF
  ELSE
     first_order=.FALSE.     
  ENDIF

#ifdef VAR_TIME
  ForcePrint = .TRUE.
#else
  IF(LSTIME_PRINT)THEN
     ForcePrint = .TRUE.
  ELSE
     ForcePrint = .FALSE.
  ENDIF
#endif

  call LSTIMER('START ',TS,TE,DECinfo%output,ForcePrint)
  LUPRI = DECinfo%output
  CALL LSTIMER('START ',TS2,TE2,LUPRI)
  ChangedDefault = .FALSE.
  !The 3 Options 
!  call time_start_phase(PHASE_WORK)   
!  call time_start_phase( PHASE_COMM )
!  call time_start_phase( PHASE_IDLE )
  call time_start_phase(PHASE_WORK)   
#ifdef VAR_TIME
  call time_phases_get_current(current_wt=phase_cntrs)
#endif

  myload = 0
  ! If MPI is not used, consider the single node to be "master"
  master=.true.
#ifdef VAR_MPI
  master= (infpar%lg_mynum == infpar%master)
  IF(.NOT.master) LUPRI = 6 !standard Output
#endif
  call LSTIMER('START',tcmpi1,twmpi1,DECinfo%output)
  if(.not. master) then  ! flop counting for slaves
     call start_flop_counter()
  end if
  ! Initialize stuff

  natoms = MyFragment%natoms
  nbasis = MyFragment%nbasis
  nocc = MyFragment%noccAOS        ! occupied AOS (only valence for frozen core)
  nvirt = MyFragment%nvirtAOS     ! virtual AOS
  noccEOS = MyFragment%noccEOS     ! occupied EOS
  nvirtEOS = MyFragment%nvirtEOS  ! virtual EOS
  nocctot = MyFragment%nocctot     ! total occ: core+valence (identical to nocc without frozen core)
  ncore = MyFragment%ncore         ! number of core orbitals
  call determine_maxBatchOrbitalsize(DECinfo%output,MyFragment%mylsitem%SETTING,MinAuxBatch,'D')
  offset = 0 
  nocc8 = nocc
  nvirt8 = nvirt
  noccEOS8 = noccEOS
  nbasis8 = nbasis
  nvirtEOS8 = nvirtEOS
  nocctot8 = nocctot

  ! For frozen core energy calculation, we never need core orbitals
  ! (but we do if first order integrals are required)
  if(DECinfo%frozencore .and. (.not. first_order)) nocctot = nocc

  ! In general, for frozen core AND first order integrals, special care must be taken
  ! No frozen core OR frozen core calculation for just energy uses the same
  ! code from now on because the frozen core approximation is "built into" the fragment,
  if(DECinfo%frozencore .and. first_order) then
     fc = .true.
     offset = ncore
  else
     fc =.false.
  end if
  if(first_order) then
   IF(DECinfo%RIMP2_Laplace)THEN
    if(DECinfo%PL>0) write(DECinfo%output,*) &
         & 'Calculating Laplace RIMP2 integrals (both energy and density) and RIMP2 amplitudes...'
   ELSE
    if(DECinfo%PL>0) write(DECinfo%output,*) 'Calculating RIMP2 integrals (both energy and density) and RIMP2 amplitudes...'
   ENDIF
  else
   IF(DECinfo%RIMP2_Laplace)THEN
    if(DECinfo%PL>0) write(DECinfo%output,*) 'Calculating Laplace RIMP2 integrals (only energy) and RIMP2 amplitudes...'
   ELSE
    if(DECinfo%PL>0) write(DECinfo%output,*) 'Calculating RIMP2 integrals (only energy) and RIMP2 amplitudes...'
   ENDIF
  end if

  !==================================================================
  ! Background memory buffering 
  !==================================================================
  use_bg_buf = .FALSE.
#ifdef VAR_MPI
  IF(DECinfo%use_bg_buffer) use_bg_buf = mem_is_background_buf_init()
#endif
  

  IF(DECinfo%AuxAtomicExtent)THEN
     call getMolecularDimensions(MyFragment%mylsitem%INPUT%AUXMOLECULE,nAtomsAux,nBasis2,nBasisAux)
  ELSE
     call getMolecularDimensions(MyFragment%mylsitem%SETTING%MOLECULE(1)%p,nAtomsAux,nBasis2,nBasisAux)
     if(natoms.NE.natomsAux)call lsquit('Error in RIMP2 natoms dim mismatch',-1)
  ENDIF
  nbasisAux8 = nbasisAux
  IF(nBasisAux.EQ.0)THEN
     WRITE(DECinfo%output,'(1X,A)')'RIMP2MEM: Warning no Aux basis have been chosen for RIMP2, Using Regular'
     ChangedDefault = .TRUE.
     call get_default_AOs(oldAORegular,oldAOdfAux) !the current values for Regular and Aux Basis 
     call set_default_AOs(oldAORegular,oldAORegular) !change to use Regular for Aux 
     call getMolecularDimensions(MyFragment%mylsitem%SETTING%MOLECULE(1)%p,nAtoms,nBasis2,nBasisAux)
  ENDIF
  if(DECinfo%PL>0)THEN
     if(master) then
        WRITE(*,'(A,2I3,4I6,1I5)')'RIMP2: DIM(nocc,noccEOS,nvirt,nvirtEOS,nbasis,nBasisAux,natoms)=',&
             & nocc,noccEOS,nvirt,nvirtEOS,nbasis,nBasisAux,natoms
     endif
  endif
!#ifndef VAR_MPI
!  do not use the Calpha or the alphaCD as a vector at any point
!  call Test_if_64bit_integer_required(nBasisAux,nocc,nvirt)
!#endif

  call Test_if_64bit_integer_required(nvirt,noccEOS,nvirt,noccEOS)
  call Test_if_64bit_integer_required(nvirtEOS,nocc,nvirtEOS,nocc)
  if(first_order) then
     call Test_if_64bit_integer_required(nvirtEOS,nocc,nvirtEOS,nvirt)
     call Test_if_64bit_integer_required(nvirt,noccEOS,noccEOS,nocctot)
  endif

  IF(.NOT.master)THEN
     IF(DECinfo%MemDebugPrint)THEN
        nsize=nocc+nvirt+nocc*noccEOS+nvirt*nvirt+nvirt*nvirtEOS+nocc*nocc
        IF(first_order)THEN
           nsize = nsize + nbasis*nvirt + nbasis*nvirtEOS + nbasis*noccEOS + nbasis*nocctot
        ENDIF
        if(fc) then
           nsize = nsize + nocctot*nocctot+nbasis*nocctot
        endif
        nsize = nsize + nbasisAux8*nbasisAux8
!        print*,'SLAVE Before Allocation of EVocc etc. nsize=',nsize
!        call stats_globalmem(6)
     ENDIF
     IF(use_bg_buf)THEN
!        IF(DECinfo%MemDebugPrint)call printBGinfo()
!        IF(DECinfo%MemDebugPrint)print*,'BG: alloc EVocc etc. ',nsize
        IF(DECinfo%RIMP2_Laplace)THEN
           call mem_pseudo_alloc(TauVirt,nvirt8,nLaplace*i8)
           call mem_pseudo_alloc(TauOcc,nocc8,nLaplace*i8)
        ENDIF
        call mem_pseudo_alloc(EVocc,nocc8)
        call mem_pseudo_alloc(EVvirt,nvirt8)
        call mem_pseudo_alloc(UoccEOST,nocc8,noccEOS8) 
        call mem_pseudo_alloc(UvirtT,nvirt8,nvirt8) 
        call mem_pseudo_alloc(UvirtEOST,nvirt8,nvirtEOS8)
        call mem_pseudo_alloc(UoccT,nocc8,nocc8) 
        IF(first_order)THEN
           call mem_pseudo_alloc(CvirtAOS,nbasis8,nvirt8)
           call mem_pseudo_alloc(CvirtEOS,nbasis8,nvirtEOS8)
           call mem_pseudo_alloc(CoccEOS,nbasis8,noccEOS8)
        ENDIF
        if(fc) then
           call mem_pseudo_alloc(UoccallT,nocctot8,nocctot8) 
           call mem_pseudo_alloc(CDIAGoccALL,nbasis8,nocctot8) 
        endif
        IF(first_order)THEN
           call mem_pseudo_alloc(CoccTmp,nbasis8,nocctot8) 
        ENDIF
        call mem_pseudo_alloc(ABdecomp,nbasisAux8,nbasisAux8)
     ELSE     
        IF(DECinfo%RIMP2_Laplace)THEN
           call mem_alloc(TauVirt,nvirt,nLaplace)
           call mem_alloc(TauOcc,nocc,nLaplace)
        ENDIF
        call mem_alloc(EVocc,nocc)
        call mem_alloc(EVvirt,nvirt)
        call mem_alloc(UoccEOST,nocc,noccEOS) 
        call mem_alloc(UvirtT,nvirt,nvirt) 
        call mem_alloc(UvirtEOST,nvirt,nvirtEOS)
        call mem_alloc(UoccT,nocc,nocc) 
        IF(first_order)THEN
           call mem_alloc(CvirtAOS,nbasis,nvirt)
           call mem_alloc(CvirtEOS,nbasis,nvirtEOS)
           call mem_alloc(CoccEOS,nbasis,noccEOS)
        ENDIF
        if(fc) then
           call mem_alloc(UoccallT,nocctot,nocctot) 
           call mem_alloc(CDIAGoccALL,nbasis,nocctot) 
        endif
        IF(first_order)THEN
           call mem_alloc(CoccTmp,nbasis,nocctot) 
        ENDIF
        call mem_alloc(ABdecomp,nbasisAux,nbasisAux)
     ENDIF
!     IF(DECinfo%MemDebugPrint)THEN
!        print*,'SLAVE After Allocation of EVocc etc. nsize=',nsize
!        call stats_globalmem(6)
!     ENDIF
  ENDIF

  IF(use_bg_buf)THEN
     !Due to the push pull mechanisme we must deallocate in the 
     !reverse order we allocate - which means I need to allocate these
     !quantities before I allocate anything else in the background buffer

     !deallocated in fragment_energy.F90 line 663 (April 2015)
     if(first_order) then
        dimvirt = [nvirtEOS,nocc,nvirtEOS,nvirt]   
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: alloc tensor blad(',dimvirt(1)*dimvirt(2)*dimvirt(3)*dimvirt(4),')'
        call tensor_ainit(blad,dimvirt,4,bg=.TRUE.)
        dimvirt = [nvirt,noccEOS,noccEOS,nocctot]   
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: alloc tensor djik(',dimvirt(1)*dimvirt(2)*dimvirt(3)*dimvirt(4),')'
        call tensor_ainit(djik,dimvirt,4,bg=.TRUE.)
     endif
     !deallocated in fragment_energy.F90 line 677 (April 2015)
     dimvirt = [nvirtEOS,nocc,nvirtEOS,nocc]    
     IF(DECinfo%MemDebugPrint)call printBGinfo()
     IF(DECinfo%MemDebugPrint)print*,'BG: alloc tensor tvirtEOS(',dimvirt(1)*dimvirt(2)*dimvirt(3)*dimvirt(4),')'
     call tensor_ainit(tvirtEOS,dimvirt,4,bg=.TRUE.)
     dimvirt = [nvirtEOS,nocc,nvirtEOS,nocctot] 
     IF(DECinfo%MemDebugPrint)call printBGinfo()
     IF(DECinfo%MemDebugPrint)print*,'BG: alloc tensor gvirtEOS(',dimvirt(1)*dimvirt(2)*dimvirt(3)*dimvirt(4),')'
     call tensor_ainit(gvirtEOS,dimvirt,4,bg=.TRUE.)
     dimocc = [nvirt,noccEOS,nvirt,noccEOS]     
     IF(DECinfo%MemDebugPrint)call printBGinfo()
     IF(DECinfo%MemDebugPrint)print*,'BG: alloc tensor toccEOS(',dimocc(1)*dimocc(2)*dimocc(3)*dimocc(4),')'
     call tensor_ainit(toccEOS,dimocc,4,bg=.TRUE.)
     IF(DECinfo%MemDebugPrint)call printBGinfo()
     IF(DECinfo%MemDebugPrint)print*,'BG: alloc tensor toccEOS(',dimocc(1)*dimocc(2)*dimocc(3)*dimocc(4),')'
     call tensor_ainit(goccEOS,dimocc,4,bg=.TRUE.)
  ENDIF

  IF(master)THEN
     IF(DECinfo%MemDebugPrint)THEN
        nsize=nocc+nvirt+nocc*noccEOS+nvirt*nvirt+nvirt*nvirtEOS+nocc*nocc
        IF(first_order)THEN
           nsize = nsize + nbasis*nvirt + nbasis*nvirtEOS + nbasis*noccEOS + nbasis*nocctot
        ENDIF
        if(fc) then
           nsize = nsize + nocctot*nocctot+nbasis*nocctot
        endif
        nsize = nsize + nbasisAux8*nbasisAux8
!        print*,'MASTER Before Allocation of EVocc etc. nsize=',nsize
!        call stats_globalmem(6)
     ENDIF
     IF(use_bg_buf)THEN
 !       IF(DECinfo%MemDebugPrint)call printBGinfo()
 !       IF(DECinfo%MemDebugPrint)print*,'BG: alloc EVocc etc. ',nsize
        IF(DECinfo%RIMP2_Laplace)THEN
           call mem_pseudo_alloc(TauVirt,nvirt8,nLaplace*i8)
           call mem_pseudo_alloc(TauOcc,nocc8,nLaplace*i8)
        ENDIF
        call mem_pseudo_alloc(EVocc,nocc8)
        call mem_pseudo_alloc(EVvirt,nvirt8)
        call mem_pseudo_alloc(UoccEOST,nocc8,noccEOS8) 
        call mem_pseudo_alloc(UvirtT,nvirt8,nvirt8) 
        call mem_pseudo_alloc(UvirtEOST,nvirt8,nvirtEOS8)
        call mem_pseudo_alloc(UoccT,nocc8,nocc8) 
        IF(first_order)THEN
           call mem_pseudo_alloc(CvirtAOS,nbasis8,nvirt8)
           call mem_pseudo_alloc(CvirtEOS,nbasis8,nvirtEOS8)
           call mem_pseudo_alloc(CoccEOS,nbasis8,noccEOS8)
        ENDIF
        if(fc) then
           call mem_pseudo_alloc(UoccallT,nocctot8,nocctot8) 
           call mem_pseudo_alloc(CDIAGoccALL,nbasis8,nocctot8) 
        endif
        IF(first_order)THEN
           call mem_pseudo_alloc(CoccTmp,nbasis8,nocctot8) 
        ENDIF
        call mem_pseudo_alloc(ABdecomp,nbasisAux8,nbasisAux8)
     ELSE     
        IF(DECinfo%RIMP2_Laplace)THEN
           call mem_alloc(TauVirt,nvirt,nLaplace)
           call mem_alloc(TauOcc,nocc,nLaplace)
        ENDIF
        call mem_alloc(EVocc,nocc)
        call mem_alloc(EVvirt,nvirt)
        call mem_alloc(UoccEOST,nocc,noccEOS) 
        call mem_alloc(UvirtT,nvirt,nvirt) 
        call mem_alloc(UvirtEOST,nvirt,nvirtEOS)
        call mem_alloc(UoccT,nocc,nocc) 
        IF(first_order)THEN
           call mem_alloc(CvirtAOS,nbasis,nvirt)
           call mem_alloc(CvirtEOS,nbasis,nvirtEOS)
           call mem_alloc(CoccEOS,nbasis,noccEOS)
        ENDIF
        if(fc) then
           call mem_alloc(UoccallT,nocctot,nocctot) 
           call mem_alloc(CDIAGoccALL,nbasis,nocctot) 
        endif
        IF(first_order)THEN
           call mem_alloc(CoccTmp,nbasis,nocctot) 
        ENDIF
        call mem_alloc(ABdecomp,nbasisAux,nbasisAux)
     ENDIF
!     IF(DECinfo%MemDebugPrint)THEN
!        print*,'MASTER After Allocation of EVocc etc. nsize=',nsize
!        call stats_globalmem(6)
!     ENDIF
  ENDIF

  CALL LSTIMER('DECRIMP2: INIT ',TS2,TE2,LUPRI,FORCEPRINT)

  ! *************************************
  ! Get arrays for transforming integrals: Cocc,Cvirt,UoccEOST,UvirtT,UvirtEOST,UoccT
  ! *************************************
  ! CDIAGocc, CDIAGvirt:  MO coefficients for basis where Fock matrix is diagonal
  ! Uocc, Uvirt: Transform from diagonal basis to local basis (and vice versa)
  ! Note: Uocc and Uvirt have indices (local,diagonal)
  
  if(fc) then
     ! Special routine for MP2 density/gradient using frozen core approx.
     ! Note: The occupied coefficient matrices CDIAGocc and Uocc and the eigenvalues EVocc
     !       refer only to the valence space, while
     !       LoccTALL, CDIAGoccTALL, and UoccALL refer to both core+valence!
     ! Special note: CDIAGoccTALL have valence orbitals BEFORE core orbitals, see
     !               get_MP2_integral_transformation_matrices_fc!
     call get_MP2_integral_transformation_matrices_fc(MyFragment,CDIAGocc, CDIAGvirt, Uocc, Uvirt, &
          & LoccTALL, CDIAGoccTALL, UoccALL, EVocc, EVvirt)
     call array2_free(LoccTALL)
  else  ! not frozen core or simple energy calculation
     call get_MP2_integral_transformation_matrices(MyFragment,CDIAGocc, CDIAGvirt, Uocc, Uvirt, &
          & EVocc, EVvirt)
  end if

  IF(DECinfo%RIMP2_Laplace)THEN
     !  tau(a,l) = exp(epsilon_A*amp_l)   !l is the laplace points
     call BuildTauVirt(TauVirt,nvirt,nLaplace,EVvirt,LaplaceAmp)
     !  tau(i,l) = exp(-epsilon_I*amp_l)   !l is the laplace points
     call BuildTauOcc(TauOcc,nocc,nLaplace,EVOcc,LaplaceAmp)
  ENDIF

  ! Extract occupied EOS indices from rows of Uocc
  call array2_extract_EOS(Uocc,MyFragment,'O','R',tmparray2)
  !make UoccEOS(noccEOS,nocc)
  M = noccEOS   !row of Input Matrix
  N = nocc      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,tmparray2%val,0.0E0_realk,UoccEOST)
  call array2_free(tmparray2)

  M = nvirt      !row of Input Matrix
  N = nvirt      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,Uvirt%val,0.0E0_realk,UvirtT)

  ! Extract virtual EOS indices from rows of Uvirt
  call array2_extract_EOS(Uvirt,MyFragment,'V','R',tmparray2)
  call array2_free(Uvirt)
  M = nvirtEOS   !row of Input Matrix
  N = nvirt      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,tmparray2%val,0.0E0_realk,UvirtEOST)
  call array2_free(tmparray2)

  M = nocc      !row of Input Matrix
  N = nocc      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,Uocc%val,0.0E0_realk,UoccT)
  call array2_free(Uocc)

!$acc enter data copyin(EVocc,EVvirt,UoccEOST,UvirtEOST,UoccT,UvirtT) async(async_id(1)) if(.not. fc)

  IF(first_order)THEN
     !CvirtAOS(nbasis,nvirt) = CDIAGvirt%val(nbasis,nvirt)*UvirtT(nvirt,nvirt)     
     M = nbasis  !rows of Output Matrix
     N = nvirt   !columns of Output Matrix
     K = nvirt   !summation dimension
     call dgemm('N','N',M,N,K,1.0E0_realk,CDIAGvirt%val,M,UvirtT,K,0.0E0_realk,CvirtAOS,M)

     !CvirtEOS(nbasis,nvirtEOS) = CDIAGvirt%val(nbasis,nvirt)*UvirtEOST(nvirt,nvirtEOS)     
     M = nbasis  !rows of Output Matrix
     N = nvirtEOS   !columns of Output Matrix
     K = nvirt   !summation dimension
     call dgemm('N','N',M,N,K,1.0E0_realk,CDIAGvirt%val,M,UvirtEOST,K,0.0E0_realk,CvirtEOS,M)

     !CoccEOS(nbasis,noccEOS) = CDIAGocc%val(nbasis,nocc)*UoccEOST(nocc,noccEOS)     
     M = nbasis    !rows of Output Matrix
     N = noccEOS   !columns of Output Matrix
     K = nocc      !summation dimension
     call dgemm('N','N',M,N,K,1.0E0_realk,CDIAGocc%val,M,UoccEOST,K,0.0E0_realk,CoccEOS,M)
  ENDIF

  if(fc) then
     M = nocctot      !row of Input Matrix
     N = nocctot      !columns of Input Matrix
     call mat_transpose(M,N,1.0E0_realk,Uoccall%val,0.0E0_realk,UoccallT)
     call array2_free(Uoccall)
     M = nocctot      !row of Input Matrix
     N = nbasis       !columns of Input Matrix
     call mat_transpose(M,N,1.0E0_realk,CDIAGoccTALL%val,0.0E0_realk,CDIAGoccALL)
     call array2_free(CDIAGoccTALL)
  endif

  if(first_order)then
     !Diagonal Ccoeff with Core First
     IF(DECinfo%frozencore)THEN
        IF(use_bg_buf)THEN
           call mem_pseudo_alloc(CDIAGoccALLcf,nocctot8,nbasis8)
        ELSE
           call mem_alloc(CDIAGoccALLcf,nocctot,nbasis)
        ENDIF
        call PlaceCoreOrbFirst2(CDIAGoccALL,nbasis,nocctot,ncore,nocc,CDIAGoccALLcf)     
        M = nbasis           !rows of Output Matrix
        N = nocctot          !columns of Output Matrix
        K = nocctot          !summation dimension
        !CDIAGoccALLcf(nbasis,nocctot)*UoccallT(nocctot,nocctot)
        call dgemm('N','N',M,N,K,1.0E0_realk,CDIAGoccALLcf,M,UoccallT,K,0.0E0_realk,CoccTmp,M)
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(CDIAGoccALLcf) 
        ELSE
           call mem_dealloc(CDIAGoccALLcf)            
        ENDIF
     ELSE
        M = nbasis        !rows of Output Matrix
        N = nocc          !columns of Output Matrix
        K = nocc          !summation dimension
        !CDIAGocc%val(nbasis,nocc)*UoccT(nocc,nocc)
        call dgemm('N','N',M,N,K,1.0E0_realk,CDIAGocc%val,M,UoccT,K,0.0E0_realk,CoccTmp,M)
     ENDIF
     nCoccTmp = nocctot
  endif
  IF(fc)call array2_free(CDIAGocc)

  CALL LSTIMER('DECRIMP2: TransMats ',TS2,TE2,LUPRI,FORCEPRINT)

!$acc enter data copyin(EVocc,EVvirt,UoccEOST,UvirtEOST,UoccT,UvirtT,UoccallT) async(async_id(1)) if(fc)

  ! *************************************************************
  ! *                    Start up MPI slaves                    *
  ! *************************************************************

#ifdef VAR_MPI

  ! Only use slave helper if there is at least one local slave available.
  if(infpar%lg_nodtot.GT.1) then
     wakeslave=.true.
     CollaborateWithSlaves=.true.
  else
     wakeslave=.false.
     CollaborateWithSlaves=.false.
  end if

  ! Master starts up slave
  StartUpSlaves: if(wakeslave .and. master) then

     bat%MaxAllowedDimAlpha = 0
     bat%MaxAllowedDimGamma = 0
     bat%virtbatch = 0
     bat%size1=0
     bat%size2=0
     bat%size3=0

     ! Sanity check
     if(.not. MyFragment%BasisInfoIsSet) then
        call lsquit('MP2_RI_EnergyContributions: &
             & Basis info for master is not set!',-1)
     end if

     call time_start_phase( PHASE_COMM )
     ! Wake up slaves to do the job: slaves awoken up with (RIMP2INAMP)
     ! and call MP2_RI_EnergyContribution_slave which calls
     ! mpi_communicate_mp2_int_and_amp and then MP2_RI_EnergyContribution.
     call ls_mpibcast(RIMP2INAMP,infpar%master,infpar%lg_comm)
     ! Communicate fragment information to slaves
     call mpi_communicate_mp2_int_and_amp(MyFragment,bat,first_order,.true.)
     call time_start_phase( PHASE_WORK )
  endif StartUpSlaves
  CALL LSTIMER('DECRIMP2: WakeSlaves ',TS2,TE2,LUPRI,FORCEPRINT)
  mynum = infpar%lg_mynum
  numnodes = infpar%lg_nodtot
  HSTATUS = 80
  CALL MPI_GET_PROCESSOR_NAME(HNAME,HSTATUS,IERR)
#else
  mynum = 0
  numnodes = 1
  wakeslave = .false.
  CollaborateWithSlaves = .false.
#endif

  CALL LSTIMER('START ',TS2,TE2,LUPRI)
  ABdecompCreate = .TRUE.
  IF(fc)THEN
     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 3
     intspec(3) = 'R' !Regular AO basis function on center 4
     intspec(4) = 'C' !Coulomb Operator
     intspec(5) = 'C' !Coulomb Operator
     call Build_CalphaMO2(MyFragment%mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
          & FORCEPRINT,CollaborateWithSlaves,CDIAGvirt%val,nvirt,CDIAGoccALL,nocctot,&
          & mynum,numnodes,Calpha,NBA,ABdecomp,ABdecompCreate,intspec,use_bg_buf)
  ELSE
     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 3
     intspec(3) = 'R' !Regular AO basis function on center 4
     intspec(4) = 'C' !Coulomb Operator
     intspec(5) = 'C' !Coulomb Operator     
     call Build_CalphaMO2(MyFragment%mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
          & FORCEPRINT,CollaborateWithSlaves,CDIAGvirt%val,nvirt,CDIAGocc%val,nocc,&
          & mynum,numnodes,Calpha,NBA,ABdecomp,ABdecompCreate,intspec,use_bg_buf)
     call array2_free(CDIAGocc)
  ENDIF

  CALL LSTIMER('DECRIMP2: CalphaMO',TS2,TE2,LUPRI,FORCEPRINT)
  call array2_free(CDIAGvirt)
  IF(first_order)THEN
     ABdecompCreate = .FALSE. !do not need to create again
  ENDIF
  !At this point we have the Calpha in the diagonal basis 

! here: wait for async updates on async handle 1
!$acc wait(async_id(1))

  !=====================================================================================
  !  Major Step 5: Generate toccEOS(nvirt,noccEOS,nvirt,noccEOS)
  !=====================================================================================
  noccOut = noccEOS
  if (DECinfo%DECNP) noccOut = nocc
  dimocc = [nvirt,noccOut,nvirt,noccOut]   ! Output order

  IF(NBA.GT.0)THEN
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)     
     IF(DECinfo%RIMP2_Laplace)THEN
        !toccEOS(a,i,b,j) = sum_l w_l*TauVirt(A,l)*TauVirt(B,l)*TauOcc(I,l)*TauOcc(J,l)*C(alpha,A,I)*C(alpha,B,J)*U(A,a)*U(B,b)*U(I,i)*U(J,j)
        !toccEOS(a,i,b,j) = sum_l w_l*Ctmp2(alpha,a,i,l)*Ctmp2(alpha,b,j,l)
        !Ctmp2(alpha,a,i,l) = TauVirt(A,l)*Ctmp(alpha,A,i,l)*U(A,a)
        !Ctmp(alpha,A,i,l) = TauOcc(I,l)*C(alpha,A,I)*U(I,i)
        CALL LSTIMER('START ',TS4,TE4,LUPRI,FORCEPRINT)
        nsize1 = noccOut*(nvirt*i8)*NBA*nLaplace
        nsize2 = noccOut*(nvirt*i8)*NBA*nLaplace
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: alloc Ctmp2(',nsize1,')'
           call mem_pseudo_alloc(Ctmp2,nsize1)
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: alloc Ctmp(',nsize2,')'
           call mem_pseudo_alloc(Ctmp,nsize2)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc Ctmp2(',nsize1,')'
           call mem_alloc(Ctmp2,nsize1)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc Ctmp(',nsize2,')'
           call mem_alloc(Ctmp,nsize2)
        ENDIF
        !Ctmp(alpha,A,i,l) = TauOcc(I,l)*C(alpha,A,I)*U(I,i)
        if (DECinfo%DECNP) then
           call BuildCtmpLaplace(Calpha,NBA,nvirt,nocc,noccOut,TauOcc,nLaplace,Ctmp,UoccT)
        else
           call BuildCtmpLaplace(Calpha,NBA,nvirt,nocc,noccOut,TauOcc,nLaplace,Ctmp,UoccEOST)
        endif
        CALL LSTIMER('RIMP2: Ctmp1o ',TS4,TE4,LUPRI,FORCEPRINT)
        !Ctmp2(alpha,a,i,l) = TauVirt(A,l)*Ctmp(alpha,A,i,l)*U(A,a)
        call BuildCtmp2Laplace(Ctmp,NBA,nvirt,nvirt,noccOut,TauVirt,nLaplace,Ctmp2,UvirtT)
        CALL LSTIMER('RIMP2: Ctmp2v ',TS4,TE4,LUPRI,FORCEPRINT)
        IF(.NOT.use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc tensor toccEOS(',dimocc(1)*dimocc(2)*dimocc(3)*dimocc(4),')'
           call tensor_ainit(toccEOS,dimocc,4)
        ENDIF
        !toccEOS(a,i,b,j) = sum_l w_l*Ctmp2(alpha,a,i,l)*Ctmp2(alpha,b,j,l)
        call BuildTampLaplace(Ctmp2,NBA,nvirt,noccOut,toccEOS%elm1,nLaplace,LaplaceW)
        CALL LSTIMER('RIMP2: TampLaplaceOcc',TS4,TE4,LUPRI,FORCEPRINT)
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(Ctmp)
           call mem_pseudo_dealloc(Ctmp2)
        ELSE
           call mem_dealloc(Ctmp)
           call mem_dealloc(Ctmp2)
        ENDIF
     ELSE
        !NON LAPLACE VERSION
        !Perform Tiling if tocc(nocc,noccOut,nvirt,nvirt) does not fit in memory
        IF(use_bg_buf)THEN
           MemInGBCollected = mem_get_bg_buf_free()*8.0E-9_realk
           IF(DECinfo%MemDebugPrint)THEN
              print*,'mem_get_bg_buf_free=',mem_get_bg_buf_free()
              print*,'MemInGBCollected',MemInGBCollected,' GB'
              call printBGinfo()
           ENDIF
        ELSE
           MemInGBCollected = 0.0E0_realk
           call get_currently_available_memory(MemInGBCollected)
           MemInGBCollected = MemInGBCollected*0.80E0_realk !80%
        ENDIF
        IF(use_bg_buf)THEN
           MemStep1 = (nocc*noccOut*nvirt*nvirt+noccOut*noccOut*nvirt*nvirt)*8.0E-9_realk !tocc,tocc2
           MemStep2 = (2*noccOut*noccOut*nvirt*nvirt)*8.0E-9_realk                        !tocc2,tocc3
           Maxsize = MAX(MemStep1,MemStep2)
        ELSE
           MemStep1 = (nocc*noccOut*nvirt*nvirt+noccOut*noccOut*nvirt*nvirt)*8.0E-9_realk                 !tocc,tocc2
           MemStep2 = (2*noccOut*noccOut*nvirt*nvirt)*8.0E-9_realk                                        !tocc2,tocc3
           Maxsize = MAX(MemStep1,MemStep2)
        ENDIF
        PerformTiling = MaxSize.GT.MemInGBCollected
        IF(PerformTiling)THEN 
           IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then
              WRITE(DECinfo%output,'(A,F10.2,A,F10.2,A)')'DECRIMP2: Performing Tiling MaxSize=',&
                   &MaxSize,' GB > memory available = ',MemInGBCollected,' GB'
           endif
           IF(MemInGBCollected.LT.(2*noccOut*noccOut*nvirt*nvirt)*8.0E-9_realk)THEN
              call lsquit('RIMP2: Not enough memory for tiling in rimp2',-1)
           ENDIF
           !When Performing tiling we need 3 intermediates of
           !nsize1 = noccOut*noccOut*nvirt*nvirt        (tocc2)
           !nsize2 = nocc*noccOut*nvirt*MaxVirtSize     
           !nsize3 = noccOut*noccOut*nvirt*MaxVirtSize
           !followed by a transformation using tocc2 and tocc3 of (noccOut*noccOut*nvirt*nvirt)
           !resulting in Memreq = noccOut*noccOut*nvirt*nvirt+(nocc+noccOut)*noccOut*nvirt*MaxVirtSize
           !MaxVirtSize = (Memreq-noccOut*noccOut*nvirt*nvirt)/((nocc+noccOut)*noccOut*nvirt)
           MaxVirtSize = MIN(nvirt,FLOOR((MemInGBCollected-(noccOut*noccOut*nvirt*nvirt)*8.0E-9_realk) &
                & /((nocc+noccOut)*noccOut*nvirt*8.0E-9_realk)))
           IF(MaxVirtSize.LT.1)call lsquit('Not enough memory for tiling in rimp2',-1)
           IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then        
              WRITE(DECinfo%output,'(A,I10)')'DECRIMP2: MaxVirtSize =',MaxVirtSize 
           endif
           nTiles =  nvirt/MaxVirtSize 
           IF(nTiles.EQ.0)PerformTiling = .FALSE.
        ELSE
           IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then
              WRITE(DECinfo%output,'(A,F10.2,A,F10.2,A)')'DECRIMP2: No Tiling MaxSize=',&
                   &MaxSize,' GB < memory available = ',MemInGBCollected,' GB'
           ENDIF
        ENDIF
#if defined(VAR_OPENACC) && defined(VAR_CUDA)
        !In case of GPU usage tocc must also fit on device memory
        call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
        IF(PerformTiling)THEN 
           !check that tilesize determine accoriding to CPU memory is valid for gpu
           IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then        
              WRITE(DECinfo%output,'(A,I12)')'DECRIMP2: The CPU requires tiling in step 5  MaxVirtSize=',MaxVirtSize        
           ENDIF
           MaxSize = (noccOut*noccOut*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccOut)*8.0E0_realk+&
                & MaxVirtSize*((nocc+noccOut)*noccOut*nvirt)*8.0E0_realk
           IF(Maxsize .GT. free_gpu)THEN
              !reduce MaxVirtSize
              MaxVirtSize = MIN(nvirt,FLOOR( (free_gpu*0.80E0_realk-(noccOut*noccOut*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccOut)*&
                   & 8.0E0_realk)/(((nocc+noccOut)*noccOut*nvirt)*8.0E0_realk))) 
              IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then
                 WRITE(DECinfo%output,'(A,I12)')'DECRIMP2: The GPU requires a smaller tiling in step 5 New MaxVirtSize=',MaxVirtSize        
              ENDIF
              nTiles =  nvirt/MaxVirtSize
              IF(nTiles.EQ.0)PerformTiling = .FALSE.
           ENDIF
        ELSE
           !determine if GPU requires tiling even if CPU does not
           MaxSize = (nocc*noccOut*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccOut)*8.0E0_realk  !in BYTES
           PerformTiling = MaxSize.GT.free_gpu*0.80E0_realk
        ENDIF
        IF(PerformTiling)THEN 
           maxsize = (noccOut*noccOut*nvirt*nvirt+NBA*nvirt*nocc+noccOut*nocc+nocc*noccOut*nvirt+noccOut*noccOut*nvirt)*8.0E0_realk
           IF(Maxsize .GT. free_gpu)THEN
              print*,'Calpha requires',NBA*nvirt*nocc*8,'Bytes'
              print*,'U requires',noccOut*nocc*8,'Bytes'
              print*,'tocc2 requires',noccOut*noccOut*nvirt*nvirt*8,'Bytes'
              print*,'tocc which requires at least ',nocc*noccOut*nvirt*8,'Bytes'
              print*,'tocc2TMP which requires at least',noccOut*noccOut*nvirt*8,'Bytes'
              print*,'In total',MaxSize,'Bytes'
              print*,'Free on the GPU: ',free_gpu,'Bytes'
              call lsquit('GPU memory cannot hold required objects')
           ENDIF
           MaxVirtSize = MIN(nvirt,FLOOR( (free_gpu*0.80E0_realk-(noccOut*noccOut*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccOut)*&
                & 8.0E0_realk)/(((nocc+noccOut)*noccOut*nvirt)*8.0E0_realk))) 
           nTiles =  nvirt/MaxVirtSize
           IF(nTiles.EQ.0)PerformTiling = .FALSE.
        ENDIF
#endif
        
#if defined(VAR_OPENACC) && defined(VAR_CUDA)
        IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then
           !In case of GPU usage tocc must also fit on device memory
           call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
           WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Memory available on device (step 5)     ',free_gpu*1.0E0_realk,' Bytes'
           IF(PerformTiling)THEN 
              WRITE(DECinfo%output,'(A,I12)')'DECRIMP2: MaxVirtSize',MaxVirtSize
              WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Memory required in Step 5 using tiling  ',&
                   & (noccOut*noccOut*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccOut)*8.0E0_realk+MaxVirtSize*((nocc+noccOut)*noccOut*nvirt)*8.0E0_realk,' Bytes'
           ELSE
              WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Memory required in Step 5 without tiling',&
                   & (noccOut*noccOut*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccOut)*8.0E0_realk+nocc*noccOut*nvirt*nvirt*8.0E0_realk,' Bytes'
           ENDIF
        endif
#endif
        if (DECinfo%RIMP2_tiling)THEN
           PerformTiling = .TRUE. ! enforce tiling
           MaxVirtSize = 1
           nTiles =  nvirt/MaxVirtSize
        ENDIF
        IF(PerformTiling)THEN
           nsize1 = noccOut*(noccOut*i8)*nvirt*(nvirt*i8)
           nsize2 = nocc*(noccOut*i8)*nvirt*(MaxVirtSize*i8)
           nsize3 = noccOut*noccOut*(nvirt*MaxVirtSize*i8)
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc tocc2(',nsize1,')'
              call mem_pseudo_alloc(tocc2,nsize1)
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc tocc(',nsize2,')'
              call mem_pseudo_alloc(tocc,nsize2)
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc tocc2TMP(',nsize3,')'
              call mem_pseudo_alloc(tocc2TMP,nsize3)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc tocc2(',nsize1,')'
              call mem_alloc(tocc2,nsize1)
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc tocc(',nsize2,')'
              call mem_alloc(tocc,nsize2)
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc tocc2TMP(',nsize3,')'
              call mem_alloc(tocc2TMP,nsize3)
           ENDIF
           !$acc enter data create(tocc,tocc2,tocc2TMP) copyin(Calpha) 
           DO I=1,nTiles
              offsetV = (I-1)*MaxVirtSize
              if (DECinfo%DECNP) then
                 call RIMP2_calc_toccA(nvirt,nocc,noccOut,NBA,Calpha,EVocc,EVvirt,tocc,UoccT,&
                      & MaxVirtSize,offsetV)
              else
                 call RIMP2_calc_toccA(nvirt,nocc,noccOut,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,&
                      & MaxVirtSize,offsetV)
              endif
              !Transform second occupied index (IDIAG,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BDIAG)
              M = noccOut              !rows of Output Matrix
              N = noccOut*nvirt*MaxVirtSize  !columns of Output Matrix
              K = nocc                 !summation dimension
              IF(DECinfo%DECNP)THEN
                 call ls_dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccT,K,tocc,K,0.0E0_realk,tocc2TMP,M,&
                      & int((i8*K)*M,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
              ELSE
                 call ls_dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2TMP,M,&
                      & int((i8*K)*M,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
              ENDIF
              call PlugInTotocc2(tocc2,noccOut,nvirt,tocc2TMP,MaxVirtSize,offsetV)
           ENDDO
           IF(MOD(nvirt,MaxVirtSize).NE.0)THEN !Remainder
              offsetV = nTiles*MaxVirtSize
              MaxVirtSize = MOD(nvirt,MaxVirtSize)
              IF(DECinfo%DECNP)THEN
                 call RIMP2_calc_toccA(nvirt,nocc,noccOut,NBA,Calpha,EVocc,EVvirt,tocc,UoccT,&
                      & MaxVirtSize,offsetV)
              ELSE
                 call RIMP2_calc_toccA(nvirt,nocc,noccOut,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,&
                      & MaxVirtSize,offsetV)
              ENDIF
              !Transform second occupied index (IDIAG,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BDIAG)
              M = noccOut                    !rows of Output Matrix
              N = noccOut*nvirt*MaxVirtSize  !columns of Output Matrix
              K = nocc                       !summation dimension
              IF(DECinfo%DECNP)THEN
                 call ls_dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccT,K,tocc,K,0.0E0_realk,tocc2TMP,M,&
                      & int((i8*K)*M,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
              ELSE
                 call ls_dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2TMP,M,&
                      & int((i8*K)*M,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
              ENDIF
              call PlugInTotocc2(tocc2,noccOut,nvirt,tocc2TMP,MaxVirtSize,offsetV)
           ENDIF
           !$acc exit data delete(tocc,tocc2TMP)
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tocc2TMP(',size(tocc2TMP),')'
              call mem_pseudo_dealloc(tocc2TMP)
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tocc(',size(tocc),')'
              call mem_pseudo_dealloc(tocc)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc tocc2TMP(',size(tocc2TMP),')'
              call mem_dealloc(tocc2TMP)
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc tocc(',size(tocc),')'
              call mem_dealloc(tocc)
           ENDIF
        ELSE
           !No Tiling
           !Calculate and partial transform to local basis:
           !transform 1 occupied indices (IDIAG,JLOC,ADIAG,BDIAG)
           offsetV=0
           nsize1 = noccOut*(noccOut*i8)*nvirt*(nvirt*i8)
           nsize2 = nocc*(noccOut*i8)*nvirt*(nvirt*i8)
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc tocc2(',nsize1,')'
              call mem_pseudo_alloc(tocc2,nsize1)
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc tocc(',nsize2,')'
              call mem_pseudo_alloc(tocc,nsize2)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc tocc2(',nsize1,')'
              call mem_alloc(tocc2,nsize1)
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc tocc(',nsize2,')'
              call mem_alloc(tocc,nsize2)
           ENDIF
           !$acc enter data create(tocc) copyin(Calpha)
           IF(DECinfo%DECNP)THEN
              call RIMP2_calc_toccA(nvirt,nocc,noccOut,NBA,Calpha,EVocc,EVvirt,tocc,UoccT,nvirt,offsetV)
           ELSE
              call RIMP2_calc_toccA(nvirt,nocc,noccOut,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,nvirt,offsetV)
           ENDIF
           !Transform second occupied index (IDIAG,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BDIAG)
           M = noccOut              !rows of Output Matrix
           N = noccOut*nvirt*nvirt  !columns of Output Matrix
           K = nocc                 !summation dimension
           !$acc enter data create(tocc2)
           IF(DECinfo%DECNP)THEN
              call ls_dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccT,K,tocc,K,0.0E0_realk,tocc2,M,&
                   & int((i8*K)*M,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
           ELSE
              call ls_dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2,M,&
                   & int((i8*K)*M,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
           ENDIF
           !$acc exit data delete(tocc)
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tocc(',size(tocc),')'
              call mem_pseudo_dealloc(tocc)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc tocc(',size(tocc),')'
              call mem_dealloc(tocc)
           ENDIF
        ENDIF
        !Transform first Virtual index (ILOC,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BLOC)
        M = noccOut*noccOut*nvirt  !rows of Output Matrix
        N = nvirt                  !columns of Output Matrix
        K = nvirt                  !summation dimension
        nsize = nvirt*(nvirt*i8)*noccOut*(i8*noccOut)
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: alloc tocc3(',nsize,')'
           call mem_pseudo_alloc(tocc3,nsize)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc tocc3(',nsize,')'
           call mem_alloc(tocc3,nsize)
        ENDIF
        !$acc enter data create(tocc3)
        call ls_dgemm_acc('N','N',M,N,K,1.0E0_realk,tocc2,M,UvirtT,K,0.0E0_realk,tocc3,M,&
             & int((i8*M)*K,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
        !$acc exit data delete(tocc2)
        !Final virtual transformation and reorder to dimocc
        IF(.NOT.use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc tocc2(',size(tocc2),')'
           call mem_dealloc(tocc2)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc tensor toccEOS(',dimocc(1)*dimocc(2)*dimocc(3)*dimocc(4),')'
           call tensor_ainit(toccEOS,dimocc,4)
        ENDIF
        !$acc enter data create(toccEOS%elm1)
        call RIMP2_calc_toccB(nvirt,noccOut,tocc3,UvirtT,toccEOS%elm1)
        !$acc exit data copyout(toccEOS%elm1) async(async_id(2))
        !$acc exit data delete(tocc3)
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tocc3(',size(tocc3),')'
           call mem_pseudo_dealloc(tocc3)     
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tocc2(',size(tocc2),')'
           call mem_pseudo_dealloc(tocc2)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc tocc3(',size(tocc3),')'
           call mem_dealloc(tocc3)     
        ENDIF
     ENDIF
     CALL LSTIMER('RIMP2: toccEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     IF(.NOT.use_bg_buf)call tensor_ainit(toccEOS,dimocc,4)
     nsize = nvirt*noccOut*nvirt*noccOut
     call ls_dzero8(toccEOS%elm1,nsize)
  ENDIF
  CALL LSTIMER('DECRIMP2: tocc          ',TS2,TE2,LUPRI,FORCEPRINT)

  !=====================================================================================
  !  Major Step 6: Generate tvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
  !=====================================================================================

  nvirtOut = nvirtEOS
  dimvirt = [nvirtOut,nocc,nvirtOut,nocc]   ! Output order
  ! For DECNP, the amplitudes are return in the full AOS space
  ! this mean that virtual and occupied partitioning are equivalent
  ! at this stage. Therefore we return only the occupied amplitudes
  not_DECNP_1: if (.not.DECinfo%DECNP) then
     IF(NBA.GT.0)THEN
        CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
        IF(DECinfo%RIMP2_Laplace)THEN
           !tvirtEOS(a,i,b,j) = sum_l w_l*TauVirt(A,l)*TauVirt(B,l)*TauOcc(I,l)*TauOcc(J,l)*C(alpha,A,I)*C(alpha,B,J)*U(A,a)*U(B,b)*U(I,i)*U(J,j)
           !tvirtEOS(a,i,b,j) = sum_l w_l*Ctmp2(alpha,a,i,l)*Ctmp2(alpha,b,j,l)
           !Ctmp2(alpha,a,i,l) = TauOcc(I,l)*Ctmp(alpha,a,I,l)*U(I,i)
           !Ctmp(alpha,a,I,l) = TauVirt(A,l)*C(alpha,A,I)*U(A,a)
           CALL LSTIMER('START ',TS4,TE4,LUPRI,FORCEPRINT)
           nsize2 = nvirtOut*(nocc*i8)*NBA*nLaplace
           nsize1 = nvirtOut*(nocc*i8)*NBA*nLaplace
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc Ctmp2(',nsize1,')'
              call mem_pseudo_alloc(Ctmp2,nsize1)
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc Ctmp(',nsize2,')'
              call mem_pseudo_alloc(Ctmp,nsize2)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc Ctmp2(',nsize1,')'
              call mem_alloc(Ctmp2,nsize1)
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc Ctmp(',nsize2,')'
              call mem_alloc(Ctmp,nsize2)
           ENDIF
           !Ctmp(alpha,a,I,l) = TauVirt(A,l)*C(alpha,A,I)*U(A,a)
           if (DECinfo%DECNP) then
              call BuildCtmpVLaplace(Calpha,NBA,nvirt,nocc,nvirtOut,TauVirt,nLaplace,Ctmp,UvirtT)
           else
              call BuildCtmpVLaplace(Calpha,NBA,nvirt,nocc,nvirtOut,TauVirt,nLaplace,Ctmp,UvirtEOST)
           endif
           CALL LSTIMER('RIMP2: Ctmp1v ',TS4,TE4,LUPRI,FORCEPRINT)
           !Ctmp2(alpha,a,i,l) = TauOcc(I,l)*Ctmp(alpha,a,I,l)*U(I,i)
           call BuildCtmpVLaplace2(Ctmp,NBA,nvirtOut,nocc,TauOcc,nLaplace,Ctmp2,UoccT)
           CALL LSTIMER('RIMP2: Ctmp2o ',TS4,TE4,LUPRI,FORCEPRINT)
           !toccEOS(a,i,b,j) = sum_l w_l*Ctmp2(alpha,a,i,l)*Ctmp2(alpha,b,j,l)
           IF(.NOT.use_bg_buf)call tensor_ainit(tvirtEOS,dimvirt,4)
           call BuildTampLaplace(Ctmp2,NBA,nvirtOut,nocc,tvirtEOS%elm1,nLaplace,LaplaceW)
           CALL LSTIMER('RIMP2: TampLaplaceVirt',TS4,TE4,LUPRI,FORCEPRINT)
           IF(use_bg_buf)THEN
              call mem_pseudo_dealloc(Ctmp)
              call mem_pseudo_dealloc(Ctmp2)
           ELSE
              call mem_dealloc(Ctmp)
              call mem_dealloc(Ctmp2)
           ENDIF
        ELSE

#if defined(VAR_OPENACC) && defined(VAR_CUDA)
           !In case of GPU usage tvirt must also fit on device memory
           call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
           IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then
              WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Memory available on device (step 6)     ',free_gpu*1.0E0_realk,' Bytes'
              WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Additional memory requirement in Step 6 ',&
                   & (nvirtEOS*nvirt*nocc*nocc+nocc*nocc*nvirtEOS*nvirtEOS)*8.0E0_realk,' Bytes'     
           endif
           IF((nvirtEOS*nvirt*nocc*nocc+nocc*nocc*nvirtEOS*nvirtEOS)*8.0E0_realk.GT.free_gpu*1.0E0_realk)THEN
              print*,'DECRIMP2: Memory available on device (step 6)     ',free_gpu*1.0E0_realk,' Bytes'
              print*,'DECRIMP2: Additional memory requirement in Step 6 ',&
                   & (nvirtEOS*nvirt*nocc*nocc+nocc*nocc*nvirtEOS*nvirtEOS)*8.0E0_realk,' Bytes'     
              call lsquit('DECRIMP2: Not enough memory on the device for step 6')
           ENDIF
#endif

           !Calculate and partial transform to local basis - transform occupied indices
           nsize1 = nocc*nocc*(nvirtOut*i8)*nvirt
           nsize2 = nocc*nocc*(nvirtOut*i8)*nvirtOut
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc tvirt2(',nsize2,')'
              call mem_pseudo_alloc(tvirt2,nsize2)
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc tvirt(',nsize1,')'
              call mem_pseudo_alloc(tvirt,nsize1) !IDIAG,JDIAG,ALOC,BDIAG        
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc tvirt2(',nsize2,')'
              call mem_alloc(tvirt2,nsize2)
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc tvirt(',nsize1,')'
              call mem_alloc(tvirt,nsize1) !IDIAG,JDIAG,ALOC,BDIAG
           ENDIF
           !$acc enter data create(tvirt)
           IF(DECinfo%DECNP)THEN
              call RIMP2_calc_tvirtA(nvirt,nocc,nvirtOut,NBA,Calpha,EVocc,EVvirt,tvirt,UvirtT)
           ELSE
              call RIMP2_calc_tvirtA(nvirt,nocc,nvirtOut,NBA,Calpha,EVocc,EVvirt,tvirt,UvirtEOST)
           ENDIF
           !$acc exit data delete(EVocc,EVvirt)
           
           !Transform first Virtual index (IDIAG,JDIAG,ALOC,BDIAG) => (IDIAG,JDIAG,ALOC,BLOC)
           M = nocc*nocc*nvirtOut     !rows of Output Matrix
           N = nvirtOut               !columns of Output Matrix
           K = nvirt                  !summation dimension
           !$acc enter data create(tvirt2)
           IF(DECinfo%DECNP)THEN
              call ls_dgemm_acc('N','N',M,N,K,1.0E0_realk,tvirt,M,UvirtT,K,0.0E0_realk,tvirt2,M,&
                   & int((i8*M)*K,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
           ELSE
              call ls_dgemm_acc('N','N',M,N,K,1.0E0_realk,tvirt,M,UvirtEOST,K,0.0E0_realk,tvirt2,M,&
                   & int((i8*M)*K,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
           ENDIF
           !$acc exit data delete(tvirt)
           nsize = nocc*nocc*(nvirtOut*i8)*nvirtOut
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tvirt(',size(tvirt),')'
              call mem_pseudo_dealloc(tvirt)
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc tvirt3(',nsize,')'
              call mem_pseudo_alloc(tvirt3,nsize)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc tvirt(',size(tvirt),')'
              call mem_dealloc(tvirt)
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc tvirt3(',nsize,')'
              call mem_alloc(tvirt3,nsize)
           ENDIF
           !Transform first occupied index (IDIAG,JDIAG,ALOC,BLOC) => (ILOC,JDIAG,ALOC,BLOC)
           M = nocc                    !rows of Output Matrix
           N = nocc*nvirtOut*nvirtOut  !columns of Output Matrix
           K = nocc                    !summation dimension
           !$acc enter data create(tvirt3)
           call ls_dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccT,K,tvirt2,M,0.0E0_realk,tvirt3,M,&
                & int((i8*K)*M,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
           !$acc exit data delete(tvirt2)
           !transform last occ index to local basis and reorder 
           IF(.NOT.use_bg_buf)call tensor_ainit(tvirtEOS,dimvirt,4)
           !$acc enter data create(tvirtEOS%elm1)
           call RIMP2_calc_tvirtB(nvirtOut,nocc,tvirt3,UoccT,tvirtEOS%elm1)
           !$acc exit data copyout(tvirtEOS%elm1) async(async_id(3))
           !$acc exit data delete(tvirt3)
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tvirt3(',size(tvirt3),')'
              call mem_pseudo_dealloc(tvirt3)
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tvirt2(',size(tvirt2),')'
              call mem_pseudo_dealloc(tvirt2)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc tvirt3(',size(tvirt3),')'
              call mem_dealloc(tvirt3)
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc tvirt2(',size(tvirt2),')'
              call mem_dealloc(tvirt2)
           ENDIF
        ENDIF
        CALL LSTIMER('RIMP2: tvirtEOS',TS3,TE3,LUPRI,FORCEPRINT)
     ELSE
        IF(.NOT.use_bg_buf)Then
           call tensor_ainit(tvirtEOS,dimvirt,4)
        ENDIF
        nSize = nvirtOut*nocc*nvirtOut*nocc
        call ls_dzero8(tvirtEOS%elm1,nSize)
     ENDIF
#if defined(VAR_OPENACC) && defined(VAR_CUDA)
     call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
     IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then
        WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Memory available on device (step 7)     ',free_gpu*1.0E0_realk,' Bytes'
        WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Additional memory requirement in Step 7 ',&
             & (nba*nvirt*noccEOS*2)*8.0E0_realk,' Bytes'     
     endif
     IF((nba*nvirt*noccEOS*2)*8.0E0_realk.GT.free_gpu*1.0E0_realk)THEN
        print*,'DECRIMP2: Memory available on device (step 7)     ',free_gpu*1.0E0_realk,' Bytes'
        print*,'DECRIMP2: Additional memory requirement in Step 7 ',(nba*nvirt*noccEOS*2)*8.0E0_realk,' Bytes'     
        call lsquit('DECRIMP2: Not enough memory on the device for step 6')
     ENDIF
#endif
  end if not_DECNP_1

  !=====================================================================================
  !  Major Step 7: Generate goccEOS(nvirt,noccEOS,nvirt,noccEOS)
  !=====================================================================================

  IF(NBA.GT.0)THEN
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
     nsize = nba*nvirt*noccOut
     IF(use_bg_buf)THEN
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha2(',nsize,')'
        call mem_pseudo_alloc(Calpha2,nsize)
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha3(',nsize,')'
        call mem_pseudo_alloc(Calpha3,nsize)
     ELSE
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha2(',nsize,')'
        call mem_alloc(Calpha2,nsize)
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha3(',nsize,')'
        call mem_alloc(Calpha3,nsize)
     ENDIF
     ! Transform Calpha(ALPHA,a,i) to local occupied index and local Virt
     ! Transform index delta to local occupied index 
     !(alphaAux;gamma,Jloc) = (alphaAux;gamma,J)*U(J,Jloc)     UoccEOST(iDIAG,iLOC)
     M = nba*nvirt  !rows of Output Matrix
     N = noccOut          !columns of Output Matrix
     K = nocc             !summation dimension
!$acc enter data create(Calpha2)
     IF(DECinfo%DECNP)THEN
        call ls_dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccT,K,0.0E0_realk,Calpha2,M,&
             & int((i8*M)*K,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
!$acc exit data delete(UoccT) if(.not. first_order)
     ELSE
        call ls_dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccEOST,K,0.0E0_realk,Calpha2,M,&
             & int((i8*M)*K,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
!$acc exit data delete(UoccEOST) if(.not. first_order)
     ENDIF

!$acc enter data create(Calpha3)
     call RIMP2_TransAlpha1(nvirt,noccOut,nba,UvirtT,Calpha2,Calpha3)

     IF(DECinfo%DECNP)THEN
!$acc exit data delete(Calpha2)
     ELSE
!$acc exit data delete(Calpha2,UvirtT) if(.NOT.first_order)
!$acc exit data delete(Calpha2) if(first_order)
     ENDIF
     IF(.NOT.use_bg_buf) call tensor_ainit(goccEOS,dimocc,4)
!$acc enter data create(goccEOS%elm1)
     call RIMP2_calc_gocc(nvirt,noccOut,NBA,Calpha3,goccEOS%elm1)
!$acc exit data copyout(goccEOS%elm1) async(async_id(4))
!$acc exit data delete(Calpha3)
     IF(use_bg_buf)THEN
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha3(',size(Calpha3),')'
        call mem_pseudo_dealloc(Calpha3)
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha2(',size(Calpha2),')'
        call mem_pseudo_dealloc(Calpha2)
     ELSE
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha3(',size(Calpha3),')'
        call mem_dealloc(Calpha3)
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha2(',size(Calpha2),')'
        call mem_dealloc(Calpha2)
     ENDIF
     CALL LSTIMER('RIMP2: goccEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     IF(.NOT.use_bg_buf)call tensor_ainit(goccEOS,dimocc,4)
     nSize = nvirt*noccOut*nvirt*noccOut
     call ls_dzero8(goccEOS%elm1,nsize)
  ENDIF

  !=====================================================================================
  !  Major Step 8: Generate gvirtEOS(nvirtEOS,nocc,nvirtEOS,nocctot)
  !=====================================================================================
  not_DECNP_2: if (.not. DECinfo%DECNP) then
#if defined(VAR_OPENACC) && defined(VAR_CUDA)
     call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
     IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then
        WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Memory available on device (step 8)     ',free_gpu*1.0E0_realk,' Bytes'
        WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Additional memory requirement in Step 8 ',&
             & (nba*nvirt*nocctot*2)*8.0E0_realk,' Bytes'     
     endif
     IF((nba*nvirt*nocctot*2)*8.0E0_realk.GT.free_gpu*1.0E0_realk)THEN
        print*,'DECRIMP2: Memory available on device (step 8)     ',free_gpu*1.0E0_realk,' Bytes'
        print*,'DECRIMP2: Additional memory requirement in Step 8 ',(nba*nvirt*nocctot*2)*8.0E0_realk,' Bytes'     
        call lsquit('DECRIMP2: Not enough memory on the device for step 6')
     ENDIF
#endif
   
     dimvirt = [nvirtOut,nocc,nvirtOut,nocctot]   ! Output order
     IF(NBA.GT.0)THEN
        CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
        IF(fc)THEN
           nsize2 = nba*nvirt*nocctot
           nsize3 = nba*nvirt*nocctot
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha2(',nsize2,')'
              call mem_pseudo_alloc(Calpha2,nsize2)
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha3(',nsize3,')'
              call mem_pseudo_alloc(Calpha3,nsize3)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha2(',nsize2,')'
              call mem_alloc(Calpha2,nsize2)
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha3(',nsize3,')'
              call mem_alloc(Calpha3,nsize3)
           ENDIF
           !Look at the MP2 code for discussion on frozen core and first_order_integrals
           !and the order of core and valence in nocctot
!$acc enter data create(Calpha3)
           call PlaceCoreOrbFirst(Calpha,NBA,nvirt,nocctot,ncore,nocc,Calpha3)
!$acc exit data delete(Calpha) if(fc .and. (.not. first_order))
           ! Transform index delta to local occupied index 
           !(alphaAux;gamma,Jloc) = (alphaAux;gamma,J)*U(J,Jloc)    
           M = nba*nvirt        !rows of Output Matrix
           N = nocctot          !columns of Output Matrix
           K = nocctot          !summation dimension
!$acc enter data create(Calpha2)
           call ls_dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha3,M,UoccallT,K,0.0E0_realk,Calpha2,M,&
                & int((i8*M)*K,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
!$acc exit data delete(UoccallT,Calpha3) if(fc .and. (.not. first_order))
!$acc exit data delete(Calpha3) if(fc .and. first_order)
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha3(',size(Calpha3),')'
              call mem_pseudo_dealloc(Calpha3)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha3(',size(Calpha3),')'
              call mem_dealloc(Calpha3)
           ENDIF
        ELSE
           nsize = nba*nvirt*nocc
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha2(',nsize,')'
              call mem_pseudo_alloc(Calpha2,nsize)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha2(',nsize,')'
              call mem_alloc(Calpha2,nsize)
           ENDIF
           ! Transform index delta to local occupied index 
           !(alphaAux;gamma,Jloc) = (alphaAux;gamma,J)*U(J,Jloc)     UoccEOST(iDIAG,iLOC)
           M = nba*nvirt        !rows of Output Matrix
           N = nocc             !columns of Output Matrix
           K = nocc             !summation dimension
!$acc enter data create(Calpha2)
           call ls_dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccT,K,0.0E0_realk,Calpha2,M,&
                & int((i8*M)*K,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
!$acc exit data delete(UoccT,Calpha) if(.not. first_order)
        ENDIF
        nsize = nba*nvirtOut*nocctot
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha3(',nsize,')'
           call mem_pseudo_alloc(Calpha3,nsize)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha3(',nsize,')'
           call mem_alloc(Calpha3,nsize)
        ENDIF
!$acc enter data create(Calpha3)
        IF(DECinfo%DECNP)THEN
           call RIMP2_TransAlpha2(nocctot,nvirt,nvirtOut,nba,UvirtT,Calpha2,Calpha3)
!$acc exit data delete(UvirtT) if(.NOT.first_order)
        ELSE
           call RIMP2_TransAlpha2(nocctot,nvirt,nvirtOut,nba,UvirtEOST,Calpha2,Calpha3)
        ENDIF
!$acc exit data delete(Calpha2,UvirtEOST) if(.not. first_order)
!$acc exit data delete(Calpha2) if(first_order)
        IF(.NOT.use_bg_buf)call tensor_ainit(gvirtEOS,dimvirt,4)
!$acc enter data create(gvirtEOS%elm1)
        call RIMP2_calc_gvirt(nvirtOut,nocctot,NBA,nocc,Calpha3,gvirtEOS%elm1,offset)
!$acc exit data copyout(gvirtEOS%elm1) async(async_id(5))
!$acc exit data delete(Calpha3)
   
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha3(',size(Calpha3),')'
           call mem_pseudo_dealloc(Calpha3)
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha2(',size(Calpha2),')'
           call mem_pseudo_dealloc(Calpha2)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha3(',size(Calpha3),')'
           call mem_dealloc(Calpha3)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha2(',size(Calpha2),')'
           call mem_dealloc(Calpha2)
        ENDIF
        IF(.NOT.first_order)THEN
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha(',size(Calpha),')'
              call mem_pseudo_dealloc(Calpha)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha(',size(Calpha),')'
              call mem_dealloc(Calpha)
           ENDIF
        ENDIF   
        CALL LSTIMER('RIMP2: gvirtEOS',TS3,TE3,LUPRI,FORCEPRINT)
     ELSE
        IF(.NOT.use_bg_buf)call tensor_ainit(gvirtEOS,dimvirt,4)
        nSize = nvirtOut*nocc*nvirtOut*nocctot
        call ls_dzero8(gvirtEOS%elm1,nsize)
     ENDIF
  else ! DECNP: dealloc stuff
     if (.not.first_order) then
        if (NBA > 0) then
           if(use_bg_buf)then
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha(',size(Calpha),')'
              call mem_pseudo_dealloc(Calpha)
           else
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha(',size(Calpha),')'
              call mem_dealloc(Calpha)
           endif
        endif
     end if
  end if not_DECNP_2

  !=====================================================================================
  !  Major Step 9: Collect toccEOS, tvirtEOS, goccEOS, and gvirtEOS
  !=====================================================================================

! here: wait for async copyouts on async handles 2-5
!$acc wait(async_id(2),async_id(3),async_id(4),async_id(5))

#ifdef VAR_MPI
  IF(CollaborateWithSlaves) then
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     ! occupied part
     nSize = nvirt*noccEOS*nvirt*noccEOS
     call lsmpi_reduction(toccEOS%elm1,nsize,infpar%master,infpar%lg_comm)
     nSize = nvirt*noccEOS*nvirt*noccEOS
     call lsmpi_reduction(goccEOS%elm1,nSize,infpar%master,infpar%lg_comm)
     ! virtual part
     if (.not. DECinfo%DECNP) then
        nSize = nvirtEOS*nocc*nvirtEOS*nocc
        call lsmpi_reduction(tvirtEOS%elm1,nSize,infpar%master,infpar%lg_comm)
        nSize = nvirtEOS*nocc*nvirtEOS*nocctot
        call lsmpi_reduction(gvirtEOS%elm1,nsize,infpar%master,infpar%lg_comm)
     end if
     call time_start_phase(PHASE_WORK)
     if (.not. Master) then
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tensor goccEOS'
        call tensor_free(goccEOS)
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tensor toccEOS'
        call tensor_free(toccEOS)
        if (.not. DECinfo%DECNP) then
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tensor gvirtEOS'
           call tensor_free(gvirtEOS)
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc tensor tvirtEOS'
           call tensor_free(tvirtEOS)
        end if
     endif
  ENDIF
#endif
  IF(first_order)THEN

     !=====================================================================================
     !  first_order prop: Generate djik(nvirtAOS,noccEOS,noccEOS,noccAOS=nocctot)
     !=====================================================================================
     dimvirt = [nvirt,noccEOS,noccEOS,nocctot]   ! Output order
     IF(NBA.GT.0)THEN
        nsize2 = nba*MAX(nvirt*noccEOS,nocctot*nocctot)
        nsize3 = nba*MAX(nvirt*noccEOS,nocc*nocctot)
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha2(',nsize2,')'
           call mem_pseudo_alloc(Calpha2,nsize2)
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha3(',nsize3,')'
           call mem_pseudo_alloc(Calpha3,nsize3)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha2(',nsize2,')'
           call mem_alloc(Calpha2,nsize2)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha3(',nsize3,')'
           call mem_alloc(Calpha3,nsize3)
        ENDIF

        CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
        !(alphaAux;nvirt,JnoccEOS) = (alphaAux;nvirt,J)*UoccEOST(J,JnoccEOS)
        M = nba*nvirt        !rows of Output Matrix
        N = noccEOS          !columns of Output Matrix
        K = nocc             !summation dimension
!$acc enter data create(Calpha2)
        call ls_dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccEOST,K,0.0E0_realk,Calpha2,M,&
             & int((i8*M)*K,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
        !(alphaAux,nvirtAOS,noccEOS) = (alphaAux;nvirt,noccEOS)*UvirtT(nvirt,nvirtAOS)
!$acc enter data create(Calpha3)
        call RIMP2_TransAlpha2(noccEOS,nvirt,nvirt,nba,UvirtT,Calpha2,Calpha3)
        CALL LSTIMER('START ',TS2,TE2,LUPRI)
        intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
        intspec(2) = 'R' !Regular AO basis function on center 3
        intspec(3) = 'R' !Regular AO basis function on center 4
        intspec(4) = 'C' !Coulomb Operator
        intspec(5) = 'C' !Coulomb Operator
        call Build_CalphaMO2(MyFragment%mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
             & FORCEPRINT,CollaborateWithSlaves,CoccEOS,noccEOS,CoccTmp,nCoccTmp,&
             & mynum,numnodes,CalphaOcc,NBA,ABdecomp,ABdecompCreate,intspec,use_bg_buf)
        CALL LSTIMER('DECRIMP2: CalphaOO',TS2,TE2,LUPRI,FORCEPRINT)
        IF(nCoccTmp.NE.nocctot)call lsquit('Error in djik dim4',-1)
        !  djikEOS(nvirtAOS,noccEOS,noccEOS,noccAOS)
        IF(.NOT.use_bg_buf)call tensor_ainit(djik,dimvirt,4)
!$acc enter data create(djik%elm1) copyin(CalphaOcc)
        call RIMP2_calc_gen4DimFO(NBA,Calpha3,nvirt,noccEOS,CalphaOcc,noccEOS,nocctot,djik%elm1)
!$acc exit data delete(Calpha3,CalphaOcc)
!$acc exit data copyout(djik%elm1) async(async_id(6))
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc CalphaOcc(',size(CalphaOcc),')'
           call mem_pseudo_dealloc(CalphaOcc)
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha3(',size(Calpha3),')'
           call mem_pseudo_dealloc(Calpha3)
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha2(',size(Calpha2),')'
           call mem_pseudo_dealloc(Calpha2)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc CalphaOcc(',size(CalphaOcc),')'
           call mem_dealloc(CalphaOcc)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha3(',size(Calpha3),')'
           call mem_dealloc(Calpha3)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha2(',size(Calpha2),')'
           call mem_dealloc(Calpha2)
        ENDIF
        CALL LSTIMER('RIMP2: djik',TS3,TE3,LUPRI,FORCEPRINT)
     ELSE
        IF(.NOT.use_bg_buf)call tensor_ainit(djik,dimvirt,4)
        nSize = nvirt*noccEOS*noccEOS*nocctot
        call ls_dzero8(djik%elm1,nsize)
     ENDIF

     !=====================================================================================
     !  first_order prop: Generate blad(nvirtEOS,noccAOS,nvirtEOS,nvirtAOS)
     !=====================================================================================  

     dimvirt = [nvirtEOS,nocc,nvirtEOS,nvirt]   ! Output order
     IF(NBA.GT.0)THEN
        CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
        nsize = nba*nvirt*MAX(nocc,nvirt)
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha2(',nsize,')'
           call mem_pseudo_alloc(Calpha2,nsize)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha2(',nsize,')'
           call mem_alloc(Calpha2,nsize)
        ENDIF
        !(alphaAux;nvirt,noccAOS) = (alphaAux;nvirt,nocc)*U(nocc,noccAOS)
        M = nba*nvirt        !rows of Output Matrix
        N = nocc             !columns of Output Matrix
        K = nocc             !summation dimension
!$acc enter data create(Calpha2)
        call ls_dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccT,K,0.0E0_realk,Calpha2,M,&
             & int((i8*M)*K,kind=8),int(K*(N*i8),kind=8),int(M*(N*i8),kind=8),async_id(1),cublas_handle)
!$acc exit data delete(UoccT) 
        
        !(alphaAux,nvirtEOS,noccAOS) = (alphaAux;nvirt,noccAOS)*Uvirt(nvirt,nvirtEOS)
        nsize = nba*nvirtEOS*nocc
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha3(',nsize,')'
           call mem_pseudo_alloc(Calpha3,nsize)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha3(',nsize,')'
           call mem_alloc(Calpha3,nsize)
        ENDIF
!$acc enter data create(Calpha3)
        call RIMP2_TransAlpha2(nocc,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha3)   

        CALL LSTIMER('START ',TS2,TE2,LUPRI)
        intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
        intspec(2) = 'R' !Regular AO basis function on center 3
        intspec(3) = 'R' !Regular AO basis function on center 4
        intspec(4) = 'C' !Coulomb Operator
        intspec(5) = 'C' !Coulomb Operator
        call Build_CalphaMO2(MyFragment%mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
             & FORCEPRINT,CollaborateWithSlaves,CvirtEOS,nvirtEOS,CvirtAOS,nvirt,&
             & mynum,numnodes,CalphaVV,NBA,ABdecomp,ABdecompCreate,intspec,use_bg_buf)
        CALL LSTIMER('DECRIMP2: CalphaVV',TS2,TE2,LUPRI,FORCEPRINT)
        !generate blad(nvirtEOS,noccAOS,nvirtEOS,nvirtAOS)
        IF(.NOT.use_bg_buf)call tensor_ainit(blad,dimvirt,4)
!$acc enter data create(blad%elm1) copyin(CalphaVV)
        call RIMP2_calc_gen4DimFO(NBA,Calpha3,nvirtEOS,nocc,CalphaVV,nvirtEOS,nvirt,blad%elm1)
!$acc exit data delete(Calpha3,CalphaVV)
!$acc exit data copyout(blad%elm1) async(async_id(7))
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc CalphaVV(',size(CalphaVV),')'
           call mem_pseudo_dealloc(CalphaVV)
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha3(',size(Calpha3),')'
           call mem_pseudo_dealloc(Calpha3)
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha2(',size(Calpha2),')'
           call mem_pseudo_dealloc(Calpha2)
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha(',size(Calpha),')'
           call mem_pseudo_dealloc(Calpha)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc CalphaVV(',size(CalphaVV),')'
           call mem_dealloc(CalphaVV)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha3(',size(Calpha3),')'
           call mem_dealloc(Calpha3)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha2(',size(Calpha2),')'
           call mem_dealloc(Calpha2)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha(',size(Calpha),')'
           call mem_dealloc(Calpha)
        ENDIF
        CALL LSTIMER('RIMP2: blad',TS3,TE3,LUPRI,FORCEPRINT)
     ELSE
        IF(.NOT.use_bg_buf)call tensor_ainit(blad,dimvirt,4)
        nSize = nvirtEOS*nocc*nvirtEOS*nvirt
        call ls_dzero8(blad%elm1,nsize)
     ENDIF

     !=====================================================================================
     !  first_order reduction: djik and blad 
     !=====================================================================================

! here: wait for async copyouts on async handles 6 and 7
!$acc wait(async_id(6),async_id(7))

#ifdef VAR_MPI
     IF(CollaborateWithSlaves) then
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        nSize = nvirt*noccEOS*noccEOS*nocctot
        call lsmpi_reduction(djik%elm1,nsize,infpar%master,infpar%lg_comm)
        nSize = nvirtEOS*nocc*nvirtEOS*nvirt
        call lsmpi_reduction(blad%elm1,nsize,infpar%master,infpar%lg_comm)
        call time_start_phase( PHASE_WORK )
        if (.not. Master) then
           call tensor_free(djik)
           call tensor_free(blad)
        endif
     ENDIF
#endif

  ENDIF

  IF(use_bg_buf)THEN
     IF(DECinfo%MemDebugPrint)call printBGinfo()
     IF(DECinfo%MemDebugPrint)print*,'BG: dealloc ABdecomp(',size(ABdecomp),')'
     call mem_pseudo_dealloc(ABdecomp) 
     IF(DECinfo%MemDebugPrint)call printBGinfo()
     IF(DECinfo%MemDebugPrint)print*,'BG: Before dealloc Evocc etc'
     !Order is important
     IF(first_order)THEN
        call mem_pseudo_dealloc(CoccTmp) 
     ENDIF
     if(fc) then
        call mem_pseudo_dealloc(CDIAGoccALL) 
        call mem_pseudo_dealloc(UoccallT) 
     endif
     IF(first_order)THEN
        call mem_pseudo_dealloc(CoccEOS)
        call mem_pseudo_dealloc(CvirtEOS)
        call mem_pseudo_dealloc(CvirtAOS)
     ENDIF
     call mem_pseudo_dealloc(UoccT) 
     call mem_pseudo_dealloc(UvirtEOST)
     call mem_pseudo_dealloc(UvirtT) 
     call mem_pseudo_dealloc(UoccEOST) 
     call mem_pseudo_dealloc(EVvirt)
     call mem_pseudo_dealloc(EVocc)
     IF(DECinfo%RIMP2_Laplace)THEN
        call mem_pseudo_dealloc(TauOcc)
        call mem_pseudo_dealloc(TauVirt)
     ENDIF
  ELSE
     IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
     IF(DECinfo%MemDebugPrint)print*,'STD: dealloc ABdecomp(',size(ABdecomp),')'
     call mem_dealloc(ABdecomp) 
     IF(DECinfo%MemDebugPrint)print*,'STD: Before dealloc EVocc etc'
     IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
     !Order is not important
     IF(DECinfo%RIMP2_Laplace)THEN
        call mem_dealloc(TauVirt)
        call mem_dealloc(TauOcc)
     ENDIF
     call mem_dealloc(EVocc)
     call mem_dealloc(EVvirt)
     call mem_dealloc(UoccEOST) 
     call mem_dealloc(UvirtT) 
     call mem_dealloc(UvirtEOST)
     call mem_dealloc(UoccT) 
     IF(first_order)THEN
        call mem_dealloc(CvirtAOS)
        call mem_dealloc(CvirtEOS)
        call mem_dealloc(CoccEOS)
        call mem_dealloc(CoccTmp) 
     ENDIF
     if(fc) then
        call mem_dealloc(UoccallT) 
        call mem_dealloc(CDIAGoccALL) 
     endif
     IF(DECinfo%MemDebugPrint)print*,'STD: After dealloc EVocc etc'
     IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
  ENDIF

#ifdef VAR_CUBLAS

  ! Destroy the CUBLAS context
  stat = cublasDestroy_v2(cublas_handle)

#endif

  ! release async handles array
  call mem_dealloc(async_id)

  call LSTIMER('START',tcmpi2,twmpi2,DECinfo%output)
  tmpidiff = twmpi2-twmpi1
#ifdef VAR_MPI
  if(DECinfo%PL>0) write(DECinfo%output,'(a,i12,g18.8)') 'RANK, TIME(s) ',infpar%mynum,tmpidiff
  if(master) write(DECinfo%output,'(1X,a,g18.8)') 'TIME INTEGRALLOOP(s) = ', tmpidiff
#endif
  if(.not. master) then
     ! effective time for slaves
     MyFragment%slavetime_work(MODEL_RIMP2) = tmpidiff
     ! FLOP count for integral loop for slaves
     call end_flop_counter(flops,gpuflops)
  end if

#ifdef VAR_MPI
  ! If slaves were not invoked
  ! then we of course skip the reduction.
  MPIcollect: if(wakeslave) then

     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     ! FLOP counting
     if(master) then
        flops=0.0E0_realk  ! we want to count only flops from slaves (these were set above)
        gpuflops = 0.0E0_realk ! we want to count only gpu flops from slaves (these were set above)
        ! Total time for all slaves (not local master itself)
        MyFragment%slavetime_work(MODEL_RIMP2)=0.0E0_realk
     end if

     call lsmpi_reduction(flops,infpar%master,infpar%lg_comm)
     call lsmpi_reduction(gpuflops,infpar%master,infpar%lg_comm)
     if(master)MyFragment%flops_slaves=flops !save flops for local slaves (not local master)
     if(master)MyFragment%gpu_flops_slaves=gpuflops !save flops for local slaves (not local master)

     ! Total time for all slaves (not local master itself)
     if(master) MyFragment%slavetime_work(MODEL_RIMP2)=0.0E0_realk
     call lsmpi_reduction(MyFragment%slavetime_work(MODEL_RIMP2),infpar%master,infpar%lg_comm)
     call time_start_phase( PHASE_WORK )
     if(.not. master) then ! SLAVE: Done with arrays and fragment
        call atomic_fragment_free(MyFragment)
     end if
  end if MPIcollect

  ! Number of MPI tasks (Could change to nAuxBasis)
  MyFragment%ntasks = nAtomsAux
#endif

  IF(DECinfo%MemDebugPrint.OR.DECinfo%PL>2)then
     WRITE(DECinfo%output,'(1X,A)')'MP2MEM: RIMP2_integrals_and_amplitudes_workhorse:'
     WRITE(DECinfo%output,'(1X,A)')'MP2MEM: Memory Statistics at the end of the subroutine'
     call stats_globalmem(DECinfo%output)
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
     FLUSH(DECinfo%output)
#endif
  endif

  CALL LSTIMER('DECRIMP2: Finalize',TS2,TE2,LUPRI,FORCEPRINT)
  call LSTIMER('DEC-RIMP2 ',TS,TE,DECinfo%output,ForcePrint)
  IF(ChangedDefault)THEN
     call set_default_AOs(oldAORegular,oldAOdfAux) !revert Changes
  ENDIF
#ifdef VAR_TIME
  call time_phases_get_diff(current_wt=phase_cntrs)
  time_w = phase_cntrs( PHASE_WORK_IDX )
  time_c = phase_cntrs( PHASE_COMM_IDX )
  time_i = phase_cntrs( PHASE_IDLE_IDX )  
  write(*,'(A,g10.3,A)')"DECRIMP2 time WORK",time_w," seconds"
  write(*,'(A,g10.3,A)')"DECRIMP2 time COMM",time_c," seconds"
  write(*,'(A,g10.3,A)')"DECRIMP2 time IDLE",time_i," seconds"
#ifdef VAR_PAPI
  CALL LS_GETTIM(CPUTIME,WALLTIMEEND)
  WTIME = WALLTIMEEND-WALLTIMESTART  
  CALL ls_TIMTXT('>>>  WALL Time used in RIMP2_integrals_and_amp is ',WTIME,LUPRI)
  papiflops=0 ! zero flops (this is probably redundant)
  call myPAPI_stop(eventset2,papiflops)
  write(LUPRI,*) 'FLOPS for RIMP2_integrals_and_amplitudes   = ', papiflops
  write(LUPRI,*) 'FLOPS/s for RIMP2_integrals_and_amplitudes = ', papiflops/WTIME
#endif
#endif

end subroutine RIMP2_integrals_and_amplitudes

subroutine BuildTauVirt(Tvirt,nvirt,nLaplace,EpsVirt,LaplaceAmp)
  implicit none
  integer,intent(in) :: nvirt,nLaplace
  real(realk),intent(in) :: EpsVirt(nvirt),LaplaceAmp(nLaplace)
  real(realk),intent(inout) :: Tvirt(nvirt,nLaplace)
  !
  integer :: L,A
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(L,&
  !$OMP A) SHARED(nvirt,nLaplace,EpsVirt,LaplaceAmp,Tvirt)
  do L=1,nLaplace
     do A=1,nvirt
        Tvirt(A,L) = exp(EpsVirt(A)*LaplaceAmp(L))
     enddo
  enddo
  !$OMP END PARALLEL DO
end subroutine BuildTauVirt

subroutine BuildTauOcc(Tocc,nocc,nLaplace,EpsOcc,LaplaceAmp)
  implicit none
  integer,intent(in) :: nocc,nLaplace
  real(realk),intent(in) :: EpsOcc(nocc),LaplaceAmp(nLaplace)
  real(realk),intent(inout) :: Tocc(nocc,nLaplace)
  !
  integer :: L,I
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(L,&
  !$OMP I) SHARED(nocc,nLaplace,EpsOcc,LaplaceAmp,Tocc)
  do L=1,nLaplace
     do I=1,nocc
        Tocc(I,L) = exp(-EpsOcc(I)*LaplaceAmp(L))
     enddo
  enddo
  !$OMP END PARALLEL DO
end subroutine BuildTauOcc

subroutine PlaceCoreOrbFirst(Calpha,NBA,nvirtEOS,nocctot,ncore,nocc,Calpha3)
  implicit none
  integer,intent(in) :: NBA,nvirtEOS,nocctot,ncore,nocc
  real(realk),intent(in) :: Calpha(NBA,nvirtEOS,nocctot)
  real(realk),intent(inout) :: Calpha3(NBA,nvirtEOS,nocctot)
  integer :: I,J,K
#ifdef VAR_OPENACC
  !$ACC PARALLEL DEFAULT(none) PRIVATE(I,K,J) &
  !$ACC COPYIN(ncore,NBA,nocc,nocctot,nvirtEOS) present(Calpha,Calpha3)
  !$ACC LOOP COLLAPSE(3) 
#else
  !$OMP PARALLEL DEFAULT(none) PRIVATE(I,K,J) &
  !$OMP SHARED(ncore,NBA,nocc,nocctot,nvirtEOS,Calpha,Calpha3)
  !$OMP DO COLLAPSE(3) 
#endif
  DO K=1,ncore
     DO J=1,nvirtEOS
        DO I=1,NBA
           Calpha3(I,J,K) = Calpha(I,J,K+nocc)
        ENDDO
     ENDDO
  ENDDO
#ifdef VAR_OPENACC
  !$ACC LOOP COLLAPSE(3) 
#else
  !$OMP END DO NOWAIT
  !$OMP DO COLLAPSE(3)
#endif
  DO K=1,nocc
     DO J=1,nvirtEOS
        DO I=1,NBA
           Calpha3(I,J,K+ncore) = Calpha(I,J,K)
        ENDDO
     ENDDO
  ENDDO
#ifdef VAR_OPENACC
  !$ACC END PARALLEL
#else
  !$OMP END DO
  !$OMP END PARALLEL
#endif
end subroutine PlaceCoreOrbFirst

subroutine PlaceCoreOrbFirst2(CC,nbasis,nocctot,ncore,nocc,CC3)
  implicit none
  integer,intent(in) :: nbasis,nocctot,ncore,nocc
  real(realk),intent(in) :: CC(nbasis,nocctot)
  real(realk),intent(inout) :: CC3(nbasis,nocctot)
  integer :: J,K
  !$OMP PARALLEL DEFAULT(none) PRIVATE(K,J) &
  !$OMP SHARED(ncore,nocc,nocctot,nbasis,CC,CC3)
  !$OMP DO COLLAPSE(2) 
  DO K=1,ncore
     DO J=1,nbasis
        CC3(J,K) = CC(J,K+nocc)
     ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  !$OMP DO COLLAPSE(2)
  DO K=1,nocc
     DO J=1,nbasis
        CC3(J,K+ncore) = CC(J,K)
     ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine PlaceCoreOrbFirst2

subroutine PlugInTotocc2(tocc2,noccEOS,nvirt,tocc2TMP,MaxVirtSize,offsetV)
  implicit none
  integer,intent(in) :: noccEOS,nvirt,MaxVirtSize,offsetV
  real(realk),intent(inout) :: tocc2(noccEOS*noccEOS*nvirt,nvirt)
  real(realk),intent(in) :: tocc2TMP(noccEOS*noccEOS*nvirt,MaxVirtSize)
  !local variables
  integer :: I,B
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC PRIVATE(I,B) &
  !$acc firstprivate(nvirt,noccEOS,MaxVirtSize,offsetV) &
  !$ACC present(tocc2,tocc2TMP)
#else
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
  !$OMP PRIVATE(I,B) &
  !$OMP SHARED(nvirt,noccEOS,MaxVirtSize,offsetV,tocc2,tocc2TMP)
#endif
  DO B=1,MaxVirtSize
     DO I=1,noccEOS*noccEOS*nvirt
        tocc2(I,B+offsetV) = tocc2TMP(I,B)
     ENDDO
  ENDDO
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
#else
  !$OMP END PARALLEL DO
#endif
end subroutine PlugInTotocc2

!alphaCD(NBA,nvirt,nocc) is in the diagonal basis 
subroutine RIMP2_calc_toccA(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,nvirt2,offset2)
  implicit none
  integer,intent(in) :: nvirt,nocc,noccEOS,NBA,nvirt2,offset2
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(in) :: EVocc(nocc),EVvirt(nvirt),UoccEOST(nocc,noccEOS)
  real(realk),intent(inout) :: tocc(nocc,noccEOS,nvirt,nvirt2)
  !
  integer :: BDIAG,ADIAG,IDIAG,JDIAG,ALPHAAUX,ILOC,JLOC
  real(realk) :: gmocont,deltaEPS,TMP,gpuflops
  real(realk) :: toccTMP(nocc)
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(3)&
  !$ACC PRIVATE(BDIAG,ADIAG,IDIAG,JDIAG,&
  !$ACC         ALPHAAUX,ILOC,JLOC,gmocont,deltaEPS,toccTMP,TMP) &
  !$acc firstprivate(nvirt,nocc,noccEOS,NBA,nvirt2,offset2) &
  !$acc present(tocc,Calpha,UoccEOST,EVocc,EVvirt)
#else
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
  !$OMP PRIVATE(BDIAG,ADIAG,IDIAG,JDIAG,&
  !$OMP         ALPHAAUX,ILOC,JLOC,gmocont,deltaEPS,toccTMP,TMP) &
  !$OMP SHARED(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,nvirt2,offset2)
#endif
  do BDIAG=1,nvirt2
     do ADIAG=1,nvirt
        do IDIAG=1,nocc
#ifdef VAR_OPENACC
           !$ACC loop seq
#endif 
           do JDIAG=1,nocc
              gmocont = 0.0E0_realk  
#ifdef VAR_OPENACC
              !$ACC loop seq
#endif 
              do ALPHAAUX=1,nba  
                 gmocont = gmocont + Calpha(ALPHAAUX,ADIAG,IDIAG)*Calpha(ALPHAAUX,offset2+BDIAG,JDIAG)
              enddo
              deltaEPS = EVocc(IDIAG)+EVocc(JDIAG)-EVvirt(offset2+BDIAG)-EVvirt(ADIAG)
              toccTMP(JDIAG)=gmocont/deltaEPS                
           enddo
#ifdef VAR_OPENACC
           !$ACC loop seq
#endif 
           do jLOC=1,noccEOS
              TMP = 0.0E0_realk
#ifdef VAR_OPENACC
              !$ACC loop seq
#endif 
              do JDIAG=1,nocc
                 TMP = TMP + toccTMP(JDIAG)*UoccEOST(jDIAG,jLOC)
              enddo
              tocc(IDIAG,JLOC,ADIAG,BDIAG) = TMP
           enddo
        enddo
     enddo
  enddo
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
!  gpuflops = NBA*nocc*nocc*nvirt*nvirt + nocc*nocc*nvirt*nvirt*noccEOS
!  call AddFLOP_FLOPonGPUaccouting(gpuflops)
#else
  !$OMP END PARALLEL DO
#endif
END subroutine RIMP2_calc_toccA

subroutine RIMP2_calc_tvirtA(nvirt,nocc,nvirtEOS,NBA,Calpha,EVocc,EVvirt,tvirt,UvirtEOST)
  implicit none
  integer,intent(in) :: nvirt,nocc,nvirtEOS,NBA
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(in) :: EVocc(nocc),EVvirt(nvirt),UvirtEOST(nvirt,nvirtEOS)
  real(realk),intent(inout) :: tvirt(nocc,nocc,nvirtEOS,nvirt)
  !
  integer :: BDIAG,ADIAG,IDIAG,JDIAG,ALPHAAUX,ALOC,BLOC
  real(realk) :: gmocont,deltaEPS,TMP,tvirtTMP(nvirt),gpuflops
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(3) &
  !$ACC& PRIVATE(BDIAG,ADIAG,IDIAG,JDIAG,&
  !$ACC&         ALPHAAUX,ALOC,BLOC,gmocont,deltaEPS,tvirtTMP,TMP)&
  !$acc& firstprivate(nvirt,nocc,nvirtEOS,NBA)&
  !$ACC& present(tvirt,Calpha,UvirtEOST,EVocc,EVvirt)
#else
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
  !$OMP PRIVATE(BDIAG,ADIAG,IDIAG,JDIAG,ALPHAAUX,ALOC,BLOC,gmocont,deltaEPS,TMP,tvirtTMP) &
  !$OMP SHARED(nvirt,nocc,nvirtEOS,NBA,Calpha,EVocc,EVvirt,tvirt,UvirtEOST)
#endif
  do JDIAG=1,nocc
     do IDIAG=1,nocc
        do BDIAG=1,nvirt
#ifdef VAR_OPENACC
           !$ACC loop seq
#endif 
           do ADIAG=1,nvirt
              gmocont = 0.0E0_realk  
#ifdef VAR_OPENACC
              !$ACC loop seq
#endif 
              do ALPHAAUX=1,nba  
                 gmocont = gmocont + Calpha(ALPHAAUX,ADIAG,IDIAG)*Calpha(ALPHAAUX,BDIAG,JDIAG)
              enddo
              deltaEPS = EVocc(IDIAG)+EVocc(JDIAG)-EVvirt(BDIAG)-EVvirt(ADIAG)
              tvirtTMP(ADIAG)=gmocont/deltaEPS                
           enddo
#ifdef VAR_OPENACC
           !$ACC loop seq
#endif 
           do ALOC=1,nvirtEOS
              TMP = 0.0E0_realk
#ifdef VAR_OPENACC
              !$ACC loop seq
#endif 
              do ADIAG=1,nvirt
                 TMP = TMP + tvirtTMP(ADIAG)*UvirtEOST(ADIAG,ALOC)
              enddo
              tvirt(IDIAG,JDIAG,ALOC,BDIAG) = TMP
           enddo
        enddo
     enddo
  enddo
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
!  gpuflops = NBA*nocc*nocc*nvirt*nvirt + nocc*nocc*nvirt*nvirt*nvirtEOS
!  call AddFLOP_FLOPonGPUaccouting(gpuflops)
#else
  !$OMP END PARALLEL DO
#endif
END subroutine RIMP2_calc_tvirtA

!tocc(occLOC,occLOC,virtDIAG,virtLOC)=(I,J,A,B) !Transform A
subroutine RIMP2_calc_toccB(nvirt,noccEOS,tocc,UvirtT,toccEOS)
  implicit none
  integer,intent(in) :: nvirt,noccEOS
  real(realk),intent(in) :: tocc(noccEOS,noccEOS,nvirt,nvirt),UvirtT(nvirt,nvirt)
  real(realk),intent(inout) :: toccEOS(nvirt,noccEOS,nvirt,noccEOS)
  !local variables
  integer :: BLOC,JLOC,ILOC,ALOC,ADIAG
  real(realk) :: TMP,gpuflops
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(4) &
  !$ACC PRIVATE(BLOC,JLOC,ILOC,ALOC,ADIAG,TMP) &
  !$acc firstprivate(nvirt,noccEOS) &
  !$acc present(tocc,UvirtT,toccEOS)
  !dir$ noblocking
#else
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(BLOC,JLOC,ILOC,ALOC,ADIAG,TMP) &
  !$OMP SHARED(nvirt,noccEOS,tocc,UvirtT,toccEOS)
#endif
  do bLOC=1,nvirt
     do jLOC=1,noccEOS
        do aLOC=1,nvirt
#ifdef VAR_OPENACC
           !dir$ noblocking
#endif
           do iLOC=1,noccEOS
              TMP = 0.0E0_realk
#ifdef VAR_OPENACC
              !$ACC loop seq
#endif 
              do ADIAG=1,nvirt
                 TMP = TMP + tocc(ILOC,JLOC,ADIAG,BLOC)*UvirtT(ADIAG,aLOC)
              enddo
              toccEOS(ALOC,ILOC,BLOC,JLOC) = TMP              
           enddo
        enddo
     enddo
  enddo
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
!  gpuflops = noccEOS*noccEOS*nvirt*nvirt*nvirt
!  call AddFLOP_FLOPonGPUaccouting(gpuflops)
#else
  !$OMP END PARALLEL DO
#endif
end subroutine RIMP2_calc_toccB

subroutine RIMP2_calc_tvirtB(nvirtEOS,nocc,tvirt,UoccT,tvirtEOS)
  implicit none
  integer,intent(in) :: nvirtEOS,nocc
  real(realk),intent(in) :: tvirt(nocc,nocc,nvirtEOS,nvirtEOS),UoccT(nocc,nocc)
  real(realk),intent(inout) :: tvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
  !local variables
  integer :: BLOC,JLOC,ILOC,ALOC,JDIAG
  real(realk) :: TMP,gpuflops
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(4) &
  !$ACC PRIVATE(BLOC,JLOC,ILOC,ALOC,JDIAG,TMP) &
  !$ACC firstprivate(nocc,nvirtEOS) &
  !$acc present(tvirt,UoccT,tvirtEOS)
  !dir$ noblocking
#else
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(BLOC,JLOC,ILOC,ALOC,JDIAG,TMP) &
  !$OMP SHARED(nocc,nvirtEOS,tvirt,UoccT,tvirtEOS)
#endif
  do jLOC=1,nocc
     do bLOC=1,nvirtEOS
        do iLOC=1,nocc
#ifdef VAR_OPENACC
           !dir$ noblocking
#endif
           do aLOC=1,nvirtEOS
              TMP = 0.0E0_realk
#ifdef VAR_OPENACC
              !$ACC loop seq
#endif
              do JDIAG=1,nocc
                 TMP = TMP + tvirt(ILOC,JDIAG,ALOC,BLOC)*UoccT(JDIAG,JLOC)
              enddo
              tvirtEOS(ALOC,ILOC,BLOC,JLOC) = TMP
           enddo
        enddo
     enddo
  enddo
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
!  gpuflops = nvirtEOS*nocc*nvirtEOS*nocc*nvirt
!  call AddFLOP_FLOPonGPUaccouting(gpuflops)
#else
  !$OMP END PARALLEL DO
#endif
end subroutine RIMP2_calc_tvirtB

subroutine RIMP2_TransAlpha1(nvirt,noccEOS,nba,UvirtT,AlphaCD4,AlphaCD5)
  implicit none
  integer,intent(in) :: nvirt,noccEOS,nba
  real(realk),intent(in) :: UvirtT(nvirt,nvirt),AlphaCD4(nba,nvirt,noccEOS) 
  real(realk),intent(inout) :: AlphaCD5(nba,nvirt,noccEOS)
  !local variables
  integer :: BLOC,JLOC,BDIAG,ALPHAAUX
  real(realk) :: TMP,gpuflops
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(3) &
  !$ACC PRIVATE(BLOC,JLOC,BDIAG,ALPHAAUX,TMP) &
  !$ACC COPYIN(nvirt,noccEOS,NBA) &
  !$acc present(AlphaCD4,AlphaCD5,UvirtT) 
#else
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(BLOC,JLOC,BDIAG,ALPHAAUX,TMP) &
  !$OMP SHARED(nvirt,noccEOS,nba,UvirtT,AlphaCD4,AlphaCD5)
#endif
  do JLOC = 1,noccEOS
     do BLOC = 1,nvirt
        do ALPHAAUX = 1,nba
           TMP = 0.0E0_realk
#ifdef VAR_OPENACC
           !$ACC loop seq
#endif 
           do BDIAG = 1,nvirt
              TMP = TMP + UvirtT(BDIAG,BLOC)*AlphaCD4(ALPHAAUX,BDIAG,JLOC)
           enddo
           AlphaCD5(ALPHAAUX,BLOC,JLOC) = TMP
        enddo
     enddo
  enddo
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
!  gpuflops = NBA*nvirt*nvirt*noccEOS
!  call AddFLOP_FLOPonGPUaccouting(gpuflops)
#else
  !$OMP END PARALLEL DO
#endif
end subroutine RIMP2_TransAlpha1

!AlphaCD5(NBA,n3,n2) = UvirtEOST(n1,n3)*AlphaCD4(NBA,n1,n2)
subroutine RIMP2_TransAlpha2(n2,n1,n3,nba,UvirtEOST,AlphaCD4,AlphaCD5)
  implicit none
  integer,intent(in) :: nba,n1,n2,n3
  real(realk),intent(in) :: UvirtEOST(n1,n3)
  real(realk),intent(in) :: AlphaCD4(NBA,n1,n2)
  real(realk),intent(inout) :: AlphaCD5(NBA,n3,n2)
  !
  integer :: JLOC,BLOC,ALPHAAUX,BDIAG
  real(realk) :: TMP
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(3) &
  !$ACC PRIVATE(BLOC,JLOC,BDIAG,ALPHAAUX,TMP) &
  !$acc firstprivate(n1,n2,n3,NBA) &
  !$acc present(AlphaCD4,AlphaCD5,UvirtEOST) 
#else
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(BLOC,JLOC,BDIAG,ALPHAAUX,TMP) &
  !$OMP SHARED(n1,n2,n3,nba,UvirtEOST,AlphaCD4,AlphaCD5)
#endif
  do JLOC = 1,n2
     do BLOC = 1,n3
        do ALPHAAUX = 1,nba
           TMP = 0.0E0_realk
#ifdef VAR_OPENACC
           !$ACC loop seq
#endif 
           do BDIAG = 1,n1
              TMP = TMP + UvirtEOST(BDIAG,BLOC)*AlphaCD4(ALPHAAUX,BDIAG,JLOC)
           enddo
           AlphaCD5(ALPHAAUX,BLOC,JLOC) = TMP
        enddo
     enddo
  enddo
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
!  gpuflops = NBA*nvirtEOS*nocctot*nvirt
!  call AddFLOP_FLOPonGPUaccouting(gpuflops)
#else
  !$OMP END PARALLEL DO
#endif
end subroutine RIMP2_TransAlpha2

subroutine RIMP2_calc_gocc(nvirt,noccEOS,NBA,Calpha3,goccEOS)
  implicit none
  integer,intent(in) :: nvirt,noccEOS,NBA
  real(realk),intent(in) :: Calpha3(NBA,nvirt,noccEOS)
  real(realk),intent(inout) :: goccEOS(nvirt,noccEOS,nvirt,noccEOS)
  !local variables
  integer :: BLOC,JLOC,ILOC,ALOC,ALPHAAUX
  real(realk) :: TMP
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(4) &
  !$ACC PRIVATE(BLOC,JLOC,ILOC,ALOC,ALPHAAUX,TMP) &
  !$ACC COPYIN(nvirt,noccEOS,NBA) &
  !$acc present(Calpha3,goccEOS)
  !dir$ noblocking
#else
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(BLOC,JLOC,ILOC,ALOC,ALPHAAUX,TMP) &
  !$OMP SHARED(nvirt,noccEOS,NBA,Calpha3,goccEOS)
#endif
  do jLOC=1,noccEOS
     do bLOC=1,nvirt
        do iLOC=1,noccEOS
           do aLOC=1,nvirt
              TMP = 0.0E0_realk
#ifdef VAR_OPENACC
              !$ACC loop seq
#endif 
              do ALPHAAUX = 1,nba
                 tmp = tmp + Calpha3(alphaAUX,ALOC,ILOC)*Calpha3(alphaAUX,BLOC,JLOC) 
              enddo
              goccEOS(ALOC,ILOC,BLOC,JLOC) = tmp
           enddo
        enddo
     enddo
  enddo
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
#else
  !$OMP END PARALLEL DO
#endif
end subroutine RIMP2_calc_gocc

subroutine RIMP2_calc_gvirt(nvirtEOS,nocctot,NBA,nocc,Calpha3,gvirtEOS,offset)
  implicit none
  integer,intent(in) :: nvirtEOS,nocctot,NBA,nocc,offset
  real(realk),intent(in) :: Calpha3(NBA,nvirtEOS,nocctot)
  real(realk),intent(inout) :: gvirtEOS(nvirtEOS,nocc,nvirtEOS,nocctot)
  !local variables
  integer :: BLOC,JLOC,ILOC,ALOC,ALPHAAUX
  real(realk) :: TMP
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(4) &
  !$ACC PRIVATE(BLOC,JLOC,ILOC,ALOC,ALPHAAUX,TMP) &
  !$acc firstprivate(nvirtEOS,nocc,NBA,nocctot,offset) &
  !$acc present(Calpha3,gvirtEOS)
#else
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(BLOC,JLOC,ILOC,ALOC,ALPHAAUX,TMP) &
  !$OMP SHARED(nvirtEOS,nocc,NBA,Calpha3,gvirtEOS,nocctot,offset)
#endif
  do jLOC=1,nocctot
     do bLOC=1,nvirtEOS
        do iLOC=1,nocc
           do aLOC=1,nvirtEOS
              TMP = 0.0E0_realk
#ifdef VAR_OPENACC
              !$ACC loop seq
#endif 
              do ALPHAAUX = 1,nba
                 tmp = tmp + Calpha3(alphaAUX,ALOC,offset+ILOC)*Calpha3(alphaAUX,BLOC,JLOC) 
              enddo
              gvirtEOS(ALOC,ILOC,BLOC,JLOC) = tmp
           enddo
        enddo
     enddo
  enddo
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
#else
  !$OMP END PARALLEL DO
#endif
end subroutine RIMP2_calc_gvirt

subroutine RIMP2_calc_gen4DimFO(NBA,Calpha3,n1,n2,Calpha4,n3,n4,djik)
  implicit none
  integer,intent(in) :: NBA,n1,n2,n3,n4
  real(realk),intent(in) :: Calpha3(NBA,n1,n2),Calpha4(NBA,n3,n4)
  real(realk),intent(inout) :: djik(n1,n2,n3,n4)
  !local variables
  integer :: A,B,C,D,ALPHAAUX
  real(realk) :: TMP
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(4) &
  !$ACC PRIVATE(a,b,c,d,ALPHAAUX,TMP) &
  !$acc firstprivate(n1,n2,n3,n4,NBA) &
  !$acc present(Calpha3,Calpha4,djik)
#else
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(A,B,C,D,ALPHAAUX,TMP) &
  !$OMP SHARED(n1,n2,n3,n4,NBA,Calpha3,Calpha4,djik)
#endif
  do d=1,n4
     do c=1,n3
        do b=1,n2
           do a=1,n1
              TMP = 0.0E0_realk
#ifdef VAR_OPENACC
              !$acc loop seq
#endif
              do ALPHAAUX = 1,nba
                 tmp = tmp + Calpha3(alphaAUX,A,B)*Calpha4(alphaAUX,C,D)
              enddo
              djik(A,B,C,D) = tmp
           enddo
        enddo
     enddo
  enddo
#ifdef VAR_OPENACC
  !$ACC END PARALLEL LOOP
#else
  !$OMP END PARALLEL DO
#endif
end subroutine RIMP2_calc_gen4DimFO

!Ctmp2(alpha,a,i,l) = TauOcc(I,l)*Ctmp(alpha,a,I,l)*U(I,i)
subroutine BuildCtmpVLaplace2(Ctmp,NBA,nvirtEOS,nocc,TauOcc,nLaplace,Ctmp2,UoccT)        
  implicit none
  integer,intent(in) ::  NBA,nvirtEOS,nocc,nLaplace
  real(realk),intent(in) :: Ctmp(NBA,nvirtEOS,nocc,nLaplace)
  real(realk),intent(in) :: TauOcc(nocc,nLaplace)
  real(realk),intent(in) :: UoccT(nocc,nocc)
  real(realk),intent(inout) :: Ctmp2(NBA,nvirtEOS,nocc,nLaplace)
  !local variables
  integer :: l,ILOC,IDIAG,A,ALPHA
  real(realk) :: TMP
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(l,TMP,ILOC,IDIAG,A,&
  !$OMP ALPHA) SHARED(Ctmp,NBA,nvirtEOS,nocc,TauOcc,nLaplace,Ctmp2,UoccT)
  DO l = 1,nLaplace
     DO ILOC=1,nocc
        DO A=1,nvirtEOS
           DO ALPHA=1,NBA
              Ctmp2(ALPHA,A,ILOC,l) = 0.0E0_realk
           ENDDO
           DO IDIAG=1,nocc
              TMP = UoccT(IDIAG,ILOC)*TauOcc(IDIAG,l)
              DO ALPHA=1,NBA
                 Ctmp2(ALPHA,A,ILOC,l) = Ctmp2(ALPHA,A,ILOC,l) + Ctmp(ALPHA,A,IDIAG,l)*TMP
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
end subroutine BuildCtmpVLaplace2

!Ctmp(alpha,a,I,l) = TauVirt(A,l)*C(alpha,A,I)*U(A,a)
subroutine BuildCtmpVLaplace(Calpha,NBA,nvirt,nocc,nvirtEOS,TauVirt,nLaplace,Ctmp,UvirtEOST)
  implicit none
  integer,intent(in) ::  NBA,nvirt,nocc,nvirtEOS,nLaplace
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(in) :: TauVirt(nvirt,nLaplace)
  real(realk),intent(in) :: UvirtEOST(nvirt,nvirtEOS)
  real(realk),intent(inout) :: Ctmp(NBA,nvirtEOS,nocc,nLaplace)
  !local variables
  integer :: l,ALOC,ADIAG,I,ALPHA
  real(realk) :: TMP
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(l,ALOC,ADIAG,I,TMP,&
  !$OMP ALPHA) SHARED(Calpha,NBA,nvirt,nocc,nvirtEOS,TauVirt,nLaplace,Ctmp,UvirtEOST)
  DO l = 1,nLaplace
     DO I=1,nocc
        DO ALOC=1,nvirtEOS
           DO ALPHA=1,NBA
              Ctmp(ALPHA,ALOC,I,l) = 0.0E0_realk
           ENDDO
           DO ADIAG=1,nvirt
              TMP = UvirtEOST(ADIAG,ALOC)*TauVirt(ADIAG,l)
              DO ALPHA=1,NBA
                 Ctmp(ALPHA,ALOC,I,l) = Ctmp(ALPHA,ALOC,I,l) + Calpha(ALPHA,ADIAG,I)*TMP
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
end subroutine BuildCtmpVLaplace

!Ctmp(alpha,A,i,l) = TauOcc(I,l)*C(alpha,A,I)*U(I,i)
subroutine BuildCtmpLaplace(Calpha,NBA,nvirt,nocc,noccEOS,TauOcc,nLaplace,Ctmp,UoccEOST)        
  implicit none
  integer,intent(in) ::  NBA,nvirt,nocc,noccEOS,nLaplace
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(in) :: TauOcc(nocc,nLaplace)
  real(realk),intent(in) :: UoccEOST(nocc,noccEOS)
  real(realk),intent(inout) :: Ctmp(NBA,nvirt,noccEOS,nLaplace)
  !local variables
  integer :: l,ILOC,IDIAG,A,ALPHA
  real(realk) :: TMP
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(l,TMP,ILOC,IDIAG,A,&
  !$OMP ALPHA) SHARED(Calpha,NBA,nvirt,nocc,noccEOS,TauOcc,nLaplace,Ctmp,UoccEOST)
  DO l = 1,nLaplace
     DO ILOC=1,noccEOS
        DO A=1,nvirt
           DO ALPHA=1,NBA
              Ctmp(ALPHA,A,ILOC,l) = 0.0E0_realk
           ENDDO
           DO IDIAG=1,nocc
              TMP = UoccEOST(IDIAG,ILOC)*TauOcc(IDIAG,l)
              DO ALPHA=1,NBA
                 Ctmp(ALPHA,A,ILOC,l) = Ctmp(ALPHA,A,ILOC,l) + Calpha(ALPHA,A,IDIAG)*TMP
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  !Could try to do a UoccTauOcc(IDIAG,ILOC,l) = UoccEOST(IDIAG,ILOC)*TauOcc(IDIAG,l)
  ! and then call a DGEMM
end subroutine BuildCtmpLaplace

!Ctmp2(alpha,a,i,l) = TauVirt(A,l)*C(alpha,A,i,l)*U(A,a)
subroutine BuildCtmp2Laplace(Ctmp,NBA,nvirt,nvirtEOS,nocc,TauVirt,nLaplace,Ctmp2,UvirtEOST)
  implicit none
  integer,intent(in) ::  NBA,nvirt,nvirtEOS,nocc,nLaplace
  real(realk),intent(in) :: Ctmp(NBA,nvirt,nocc,nLaplace)
  real(realk),intent(in) :: TauVirt(nvirt,nLaplace)
  real(realk),intent(in) :: UvirtEOST(nvirt,nvirtEOS)
  real(realk),intent(inout) :: Ctmp2(NBA,nvirtEOS,nocc,nLaplace)
  !local variables
  integer :: l,ILOC,ADIAG,ALOC,ALPHA
  real(realk) :: TMP
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(l,TMP,ILOC,ADIAG,ALOC,&
  !$OMP ALPHA) SHARED(Ctmp,NBA,nvirt,nvirtEOS,nocc,TauVirt,nLaplace,Ctmp2,UvirtEOST)
  DO l = 1,nLaplace
     DO ALOC=1,nvirtEOS
        DO ILOC=1,nocc
           DO ALPHA=1,NBA
              Ctmp2(ALPHA,ALOC,ILOC,l) = 0.0E0_realk
           ENDDO
           DO ADIAG=1,nvirt
              TMP = UvirtEOST(ADIAG,ALOC)*TauVirt(ADIAG,l)
              DO ALPHA=1,NBA
                 Ctmp2(ALPHA,ALOC,ILOC,l) = Ctmp2(ALPHA,ALOC,ILOC,l) + Ctmp(ALPHA,ADIAG,ILOC,l)*TMP
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
end subroutine BuildCtmp2Laplace

!toccEOS(a,i,b,j) = sum_l w_l*Ctmp2(alpha,a,i,l)*Ctmp2(alpha,b,j,l)
subroutine BuildTampLaplace(Ctmp2,NBA,nvirt,noccEOS,toccEOS,nLaplace,LaplaceW)
  implicit none
  integer,intent(in) :: NBA,nvirt,noccEOS,nLaplace
  real(realk),intent(in) :: LaplaceW(nLaplace)
  real(realk),intent(in) :: Ctmp2(NBA,nvirt*noccEOS,nLaplace)
  real(realk),intent(inout) :: toccEOS(nvirt*noccEOS,nvirt*noccEOS)
  !local variables
  integer :: AI,BJ,ALPHA,L
  real(realk) :: TMP
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(TMP,AI,BJ,ALPHA,&
  !$OMP L) SHARED(Ctmp2,NBA,nvirt,noccEOS,toccEOS,nLaplace,LaplaceW)
  DO BJ=1,nvirt*noccEOS
     !L=1
     !Travel sequential through Ctmp2 LHS for each BJ while repeating RHS Ctmp2(1:ALPHA) nvirt*noccEOS times for each BJ
     DO AI=1,nvirt*noccEOS
        TMP = 0.0E0_realk
        DO ALPHA=1,NBA 
           TMP = TMP + Ctmp2(ALPHA,AI,1)*Ctmp2(ALPHA,BJ,1)*LaplaceW(1)
        ENDDO
        toccEOS(AI,BJ) = - TMP
     ENDDO
     DO L=2,nLaplace
        DO AI=1,nvirt*noccEOS
           TMP = 0.0E0_realk
           DO ALPHA=1,NBA 
              TMP = TMP + Ctmp2(ALPHA,AI,L)*Ctmp2(ALPHA,BJ,L)*LaplaceW(L)
           ENDDO
           toccEOS(AI,BJ) = toccEOS(AI,BJ) - TMP
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  
end subroutine BuildTampLaplace

end module rimp2_module


#ifdef VAR_MPI
!> \brief MPI Slave routine for RIMP2_integrals_and_amplitudes_energy.
!> The slave gets fragment information and other information from master rank,
!> then calls MP2_RI_EnergyContribution its specific components
!> \author Thomas Kjaergaard
!> \March 2014
subroutine RIMP2_integrals_and_amplitudes_slave()
  use precision
  use infpar_module
  use dec_typedef_module
  use lstiming
  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
  use lsmpi_type, only: ls_mpibcast
  use decmpi_module, only: mpi_communicate_mp2_int_and_amp
  use rimp2_module,only: RIMP2_integrals_and_amplitudes
  use tensor_type_def_module, only: tensor
  implicit none
  !> Fragment information
  type(decfrag) :: MyFragment
  !> Batch sizes
  type(mp2_batch_construction) :: bat
  !> Integrals for occ EOS: (d j|c i) in the order (d,j,c,i) [see notation inside]
  type(tensor) :: goccEOS
  !> Amplitudes for occ EOS in the order (d,j,c,i) [see notation inside]
  type(tensor) :: toccEOS
  !> Integrals for virt EOS: (b l|a k) in the order (b,l,a,k) [see notation inside]
  type(tensor) :: gvirtEOS
  !> Amplitudes for virt EOS in the order (b,l,a,k) [see notation inside]
  type(tensor) :: tvirtEOS
  !> Occ EOS integrals (d j | i k) in the order (d,j,i,k)  
  type(tensor) :: djik  !V_aos*O_eos*O_eos*O_aos
  !> Virt EOS integrals (b l | a d) in the order (b,l,a,d)  
  type(tensor) :: blad  !V_eos*O_aos*V_eos*V_aos
  !> Calculate intgrals for first order MP2 properties?
  logical :: first_order

  ! Receive fragment structure and other information from master rank
  ! *****************************************************************
  call time_start_phase( PHASE_COMM)
  call mpi_communicate_mp2_int_and_amp(MyFragment,bat,first_order,.true.)
  call time_start_phase( PHASE_WORK)
  IF(first_order)THEN
     call RIMP2_integrals_and_amplitudes(MyFragment,&
          & goccEOS,toccEOS,gvirtEOS,tvirtEOS,djik,blad)
  ELSE
     ! Calculate contribution to integrals/amplitudes for slave
     ! ********************************************************
     call RIMP2_integrals_and_amplitudes(MyFragment,&
          & goccEOS,toccEOS,gvirtEOS,tvirtEOS)
  ENDIF

end subroutine RIMP2_integrals_and_amplitudes_slave
#endif
