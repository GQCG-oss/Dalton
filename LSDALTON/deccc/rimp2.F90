!> @file
!> Two-electron integrals and amplitudes for MP2.
!> \ author Kasper Kristensen

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
  use iso_c_binding
#ifdef VAR_OPENACC
  use openacc
#endif
#if defined(VAR_CUDA) || defined(VAR_OPENACC)
  use gpu_interfaces
#endif

contains
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
  integer :: nbasis,nocc,nvirt, noccEOS, nvirtEOS,nocctot,ncore
  integer :: alpha,gamma,beta,delta,info,mynum,numnodes,nb
  integer :: IDIAG,JDIAG,ADIAG,BDIAG,ALPHAAUX,myload,nb2,natomsAux
  integer :: ILOC,JLOC,ALOC,BLOC,M,N,K,nAtoms,nbasis2,nbasisAux
  logical :: fc,ForcePrint,master,wakeslave
  logical :: CollaborateWithSlaves
  real(realk),pointer :: Calpha(:),Calpha2(:),Calpha3(:),CalphaVV(:)
  real(realk),pointer :: UoccEOS(:,:),UvirtEOS(:,:)
  real(realk),pointer :: tocc(:),UoccEOST(:,:),UvirtT(:,:),tocc3(:)
  real(realk),pointer :: toccTMP(:,:),TMPAlphaBeta_minus_sqrt(:,:),tocc2(:)
  real(realk),pointer :: tvirtTMP(:,:),tvirt(:),UoccT(:,:),UvirtEOST(:,:)
  real(realk),pointer :: tvirt2(:),tvirt3(:),Calpha4(:)
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
  logical :: ABdecompCreate
  integer,pointer :: IPVT(:)
  integer,pointer :: nbasisAuxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
  TYPE(MOLECULARORBITALINFO) :: orbitalInfo
  real(realk), pointer   :: work1(:),Etmp2222(:)
  real(realk)            :: RCOND
  integer(kind=ls_mpik)  :: COUNT,TAG,IERR,request,Receiver,sender,J,COUNT2
  real(realk) :: TS,TE,TS2,TE2,TS3,TE3
  real(realk) :: tcpu_start,twall_start, tcpu_end,twall_end,MemEstimate,memstep2
  integer ::CurrentWait(2),nAwaitDealloc,iAwaitDealloc,oldAORegular,oldAOdfAux
  integer :: MaxVirtSize,nTiles,offsetV,offset,MinAuxBatch
  logical :: useAlphaCD5,useAlphaCD6,ChangedDefault,first_order,PerformTiling
  logical :: use_bg_buf
  integer(kind=ls_mpik)  :: request5,request6
  real(realk) :: phase_cntrs(nphases),bytes_to_alloc,MinMem
  integer(kind=long) :: nSize,nsize1,nsize2,nsize3
  character :: intspec(4)
  TYPE(MoleculeInfo),pointer :: molecule1,molecule2,molecule3,molecule4
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
#ifdef VAR_MPI
  INTEGER(kind=ls_mpik) :: HSTATUS
  CHARACTER*(MPI_MAX_PROCESSOR_NAME) ::  HNAME
  TAG = 1319
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
     if(DECinfo%PL>0) write(DECinfo%output,*) 'Calculating RIMP2 integrals (both energy and density) and RIMP2 amplitudes...'
  else
     if(DECinfo%PL>0) write(DECinfo%output,*) 'Calculating RIMP2 integrals (only energy) and RIMP2 amplitudes...'
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

  CALL LSTIMER('DECRIMP2: INIT ',TS2,TE2,LUPRI,FORCEPRINT)

  ! *************************************
  ! Get arrays for transforming integrals: Cocc,Cvirt,UoccEOST,UvirtT,UvirtEOST,UoccT
  ! *************************************
  ! CDIAGocc, CDIAGvirt:  MO coefficients for basis where Fock matrix is diagonal
  ! Uocc, Uvirt: Transform from diagonal basis to local basis (and vice versa)
  ! Note: Uocc and Uvirt have indices (local,diagonal)
  call mem_alloc(EVocc,nocc)
  call mem_alloc(EVvirt,nvirt)
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

  ! Extract occupied EOS indices from rows of Uocc
  call array2_extract_EOS(Uocc,MyFragment,'O','R',tmparray2)

  !make UoccEOS(noccEOS,nocc)
  call mem_alloc(UoccEOST,nocc,noccEOS) 
  M = noccEOS   !row of Input Matrix
  N = nocc      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,tmparray2%val,0.0E0_realk,UoccEOST)
  call array2_free(tmparray2)

  call mem_alloc(UvirtT,nvirt,nvirt) 
  M = nvirt      !row of Input Matrix
  N = nvirt      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,Uvirt%val,0.0E0_realk,UvirtT)

  ! Extract virtual EOS indices from rows of Uvirt
  call array2_extract_EOS(Uvirt,MyFragment,'V','R',tmparray2)
  call array2_free(Uvirt)
  call mem_alloc(UvirtEOST,nvirt,nvirtEOS)
  M = nvirtEOS   !row of Input Matrix
  N = nvirt      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,tmparray2%val,0.0E0_realk,UvirtEOST)
  call array2_free(tmparray2)

  call mem_alloc(UoccT,nocc,nocc) 
  M = nocc      !row of Input Matrix
  N = nocc      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,Uocc%val,0.0E0_realk,UoccT)
  call array2_free(Uocc)

!$acc enter data copyin(EVocc,EVvirt,UoccEOST,UvirtEOST,UoccT,UvirtT) async(async_id(1)) if(.not. fc)

  if(fc) then
     call mem_alloc(UoccallT,nocctot,nocctot) 
     M = nocctot      !row of Input Matrix
     N = nocctot      !columns of Input Matrix
     call mat_transpose(M,N,1.0E0_realk,Uoccall%val,0.0E0_realk,UoccallT)
     call array2_free(Uoccall)

     call mem_alloc(CDIAGoccALL,nocctot,nbasis) 
     M = nocctot      !row of Input Matrix
     N = nbasis       !columns of Input Matrix
     call mat_transpose(M,N,1.0E0_realk,CDIAGoccTALL%val,0.0E0_realk,CDIAGoccALL)
     call array2_free(CDIAGoccTALL)
  endif
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
  call mem_alloc(ABdecomp,nbasisAux,nbasisAux)
  ABdecompCreate = .TRUE.
  IF(fc)THEN
!     call Build_CalphaMO(MyFragment%mylsitem,master,nbasis,nbasisAux,LUPRI,&
!          & FORCEPRINT,CollaborateWithSlaves,CDIAGvirt%val,nvirt,&
!          & CDIAGoccALL,nocctot,mynum,numnodes,nAtomsAux,Calpha,NBA,&
!          & ABdecomp,ABdecompCreate)
     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 3
     intspec(3) = 'R' !Regular AO basis function on center 4
     intspec(4) = 'C' !Coulomb Operator
     call Build_CalphaMO2(MyFragment%mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
          & FORCEPRINT,CollaborateWithSlaves,CDIAGvirt%val,nvirt,CDIAGoccALL,nocctot,&
          & mynum,numnodes,Calpha,NBA,ABdecomp,ABdecompCreate,intspec,use_bg_buf)
  ELSE
!     call Build_CalphaMO(MyFragment%mylsitem,master,nbasis,nbasisAux,LUPRI,&
!          & FORCEPRINT,CollaborateWithSlaves,CDIAGvirt%val,nvirt,&
!          & CDIAGocc%val,nocc,mynum,numnodes,nAtomsAux,Calpha,NBA,&
!          & ABdecomp,ABdecompCreate)
     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 3
     intspec(3) = 'R' !Regular AO basis function on center 4
     intspec(4) = 'C' !Coulomb Operator
     call Build_CalphaMO2(MyFragment%mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
          & FORCEPRINT,CollaborateWithSlaves,CDIAGvirt%val,nvirt,CDIAGocc%val,nocc,&
          & mynum,numnodes,Calpha,NBA,ABdecomp,ABdecompCreate,intspec,use_bg_buf)
  ENDIF

  CALL LSTIMER('DECRIMP2: CalphaMO',TS2,TE2,LUPRI,FORCEPRINT)
  IF(first_order)THEN
     ABdecompCreate = .FALSE. !do not need to create again
  ELSE
     call array2_free(CDIAGvirt)
     call array2_free(CDIAGocc)
     call mem_dealloc(ABdecomp)
  ENDIF
  !At this point we have the Calpha in the diagonal basis 

! here: wait for async updates on async handle 1
!$acc wait(async_id(1))

  !=====================================================================================
  !  Major Step 5: Generate toccEOS(nvirt,noccEOS,nvirt,noccEOS)
  !=====================================================================================

  dimocc = [nvirt,noccEOS,nvirt,noccEOS]   ! Output order
  IF(NBA.GT.0)THEN
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)     
     !Perform Tiling if tocc(nocc,noccEOS,nvirt,nvirt) does not fit in memory
     IF(use_bg_buf)THEN
        MemInGBCollected = mem_get_bg_buf_free()*8.0E-9_realk
     ELSE
        MemInGBCollected = 0.0E0_realk
        call get_currently_available_memory(MemInGBCollected)
        MemInGBCollected = MemInGBCollected*0.80E0_realk !80%
     ENDIF
     MaxSize = (nocc*noccEOS*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccEOS)*8.0E-9_realk
     PerformTiling = MaxSize.GT.MemInGBCollected
     if(DECinfo%PL>2)then
        WRITE(DECinfo%output,'(A,F10.2,A,F10.2,A)')'DECRIMP2: Perform Tiling  MaxSize=',&
             &MaxSize,' GB > memory available = ',MemInGBCollected,' GB'
     endif
     IF(PerformTiling)THEN 
        MaxVirtSize = MIN(nvirt,FLOOR((MemInGBCollected-&
             & (NBA*nvirt*nocc+nocc*noccEOS)*8.0E-9_realk)/(nocc*noccEOS*nvirt*8.0E-9_realk)))
        if(DECinfo%PL>2)then
           WRITE(DECinfo%output,'(A,I10)')'DECRIMP2: MaxVirtSize =',MaxVirtSize 
        endif
        nTiles =  nvirt/MaxVirtSize 
        IF(nTiles.EQ.0)PerformTiling = .FALSE.
     ENDIF
#if defined(VAR_OPENACC) && defined(VAR_CUDA)
     !In case of GPU usage tocc must also fit on device memory
     call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
     IF(PerformTiling)THEN 
        !check that tilesize determine accoriding to CPU memory is valid for gpu
        IF(DECinfo%PL>2)then
           WRITE(DECinfo%output,'(A,I12)')'DECRIMP2: The CPU requires tiling in step 5  MaxVirtSize=',MaxVirtSize        
        ENDIF
        MaxSize = (noccEOS*noccEOS*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccEOS)*8.0E0_realk+&
             & MaxVirtSize*((nocc+noccEOS)*noccEOS*nvirt)*8.0E0_realk
        IF(Maxsize .GT. free_gpu)THEN
           !reduce MaxVirtSize
           MaxVirtSize = MIN(nvirt,FLOOR( (free_gpu*0.80E0_realk-(noccEOS*noccEOS*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccEOS)*&
                & 8.0E0_realk)/(((nocc+noccEOS)*noccEOS*nvirt)*8.0E0_realk))) 
           IF(DECinfo%PL>2)then
              WRITE(DECinfo%output,'(A,I12)')'DECRIMP2: The GPU requires a smaller tiling in step 5 New MaxVirtSize=',MaxVirtSize        
           ENDIF
           nTiles =  nvirt/MaxVirtSize
           IF(nTiles.EQ.0)PerformTiling = .FALSE.
        ENDIF
     ELSE
        !determine if GPU requires tiling even if CPU does not
        MaxSize = (nocc*noccEOS*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccEOS)*8.0E0_realk  !in BYTES
        PerformTiling = MaxSize.GT.free_gpu*0.80E0_realk
     ENDIF
     IF(PerformTiling)THEN 
        maxsize = (noccEOS*noccEOS*nvirt*nvirt+NBA*nvirt*nocc+noccEOS*nocc+nocc*noccEOS*nvirt+noccEOS*noccEOS*nvirt)*8.0E0_realk
        IF(Maxsize .GT. free_gpu)THEN
           print*,'Calpha requires',NBA*nvirt*nocc*8,'Bytes'
           print*,'U requires',noccEOS*nocc*8,'Bytes'
           print*,'tocc2 requires',noccEOS*noccEOS*nvirt*nvirt*8,'Bytes'
           print*,'tocc which requires at least ',nocc*noccEOS*nvirt*8,'Bytes'
           print*,'tocc2TMP which requires at least',noccEOS*noccEOS*nvirt*8,'Bytes'
           print*,'In total',MaxSize,'Bytes'
           print*,'Free on the GPU: ',free_gpu,'Bytes'
           call lsquit('GPU memory cannot hold required objects')
        ENDIF
        MaxVirtSize = MIN(nvirt,FLOOR( (free_gpu*0.80E0_realk-(noccEOS*noccEOS*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccEOS)*&
             & 8.0E0_realk)/(((nocc+noccEOS)*noccEOS*nvirt)*8.0E0_realk))) 
        nTiles =  nvirt/MaxVirtSize
        IF(nTiles.EQ.0)PerformTiling = .FALSE.
     ENDIF
#endif

#if defined(VAR_OPENACC) && defined(VAR_CUDA)
     if(DECinfo%PL>2)then
        !In case of GPU usage tocc must also fit on device memory
        call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
        WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Memory available on device (step 5)     ',free_gpu*1.0E0_realk,' Bytes'
        IF(PerformTiling)THEN 
           WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: MaxVirtSize',MaxVirtSize
           WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Memory required in Step 5 using tiling  ',&
                & (noccEOS*noccEOS*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccEOS)*8.0E0_realk+MaxVirtSize*((nocc+noccEOS)*noccEOS*nvirt)*8.0E0_realk,' Bytes'
        ELSE
           WRITE(DECinfo%output,'(A,F18.2,A)')'DECRIMP2: Memory required in Step 5 without tiling',&
                & (noccEOS*noccEOS*nvirt*nvirt+NBA*nvirt*nocc+nocc*noccEOS)*8.0E0_realk+nocc*noccEOS*nvirt*nvirt*8.0E0_realk,' Bytes'
        ENDIF
     endif
#endif
     if (DECinfo%RIMP2_tiling)THEN
        PerformTiling = .TRUE. ! enforce tiling
        MaxVirtSize = 1
        nTiles =  nvirt/MaxVirtSize
     ENDIF
     IF(PerformTiling)THEN
        nsize1 = noccEOS*(noccEOS*i8)*nvirt*(nvirt*i8)
        nsize2 = nocc*(noccEOS*i8)*nvirt*(MaxVirtSize*i8)
        nsize3 = noccEOS*noccEOS*(nvirt*MaxVirtSize*i8)
        IF(use_bg_buf)THEN
           call mem_pseudo_alloc(tocc2,nsize1)
           call mem_pseudo_alloc(tocc,nsize2)
           call mem_pseudo_alloc(tocc2TMP,nsize3)
        ELSE
           call mem_alloc(tocc2,nsize1)
           call mem_alloc(tocc,nsize2)
           call mem_alloc(tocc2TMP,nsize3)
        ENDIF
!$acc enter data create(tocc,tocc2,tocc2TMP) copyin(Calpha) 
        DO I=1,nTiles
           offsetV = (I-1)*MaxVirtSize
           call RIMP2_calc_toccA(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,&
                & MaxVirtSize,offsetV)
           !Transform second occupied index (IDIAG,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BDIAG)
           M = noccEOS              !rows of Output Matrix
           N = noccEOS*nvirt*MaxVirtSize  !columns of Output Matrix
           K = nocc                 !summation dimension
#ifdef VAR_OPENACC
!$acc host_data use_device(tocc,UoccEOST,tocc2TMP)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
           call dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2TMP,M)
#elif defined(VAR_CUBLAS)
           stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                & 1.0E0_realk,c_loc(UoccEOST),int(K,kind=4),c_loc(tocc),int(K,kind=4),&
                & 0.0E0_realk,c_loc(tocc2TMP),int(M,kind=4))
#endif
!$acc end host_data
#else
           call dgemm('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2TMP,M)
#endif
           call PlugInTotocc2(tocc2,noccEOS,nvirt,tocc2TMP,MaxVirtSize,offsetV)
        ENDDO
        IF(MOD(nvirt,MaxVirtSize).NE.0)THEN !Remainder
           offsetV = nTiles*MaxVirtSize
           MaxVirtSize = MOD(nvirt,MaxVirtSize)
           call RIMP2_calc_toccA(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,&
                & MaxVirtSize,offsetV)
           !Transform second occupied index (IDIAG,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BDIAG)
           M = noccEOS                    !rows of Output Matrix
           N = noccEOS*nvirt*MaxVirtSize  !columns of Output Matrix
           K = nocc                       !summation dimension
#ifdef VAR_OPENACC
!$acc host_data use_device(tocc,UoccEOST,tocc2TMP)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
           call dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2TMP,M)
#elif defined(VAR_CUBLAS)
           stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                & 1.0E0_realk,c_loc(UoccEOST),int(K,kind=4),c_loc(tocc),int(K,kind=4),&
                & 0.0E0_realk,c_loc(tocc2TMP),int(M,kind=4))
#endif
!$acc end host_data
#else
           call dgemm('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2TMP,M)
#endif
           call PlugInTotocc2(tocc2,noccEOS,nvirt,tocc2TMP,MaxVirtSize,offsetV)
        ENDIF
!$acc exit data delete(tocc,tocc2TMP)
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(tocc2TMP)
           call mem_pseudo_dealloc(tocc)
        ELSE
           call mem_dealloc(tocc2TMP)
           call mem_dealloc(tocc)
        ENDIF
     ELSE
        !Calculate and partial transform to local basis:
        !transform 1 occupied indices (IDIAG,JLOC,ADIAG,BDIAG)
        offsetV=0
        nsize1 = noccEOS*(noccEOS*i8)*nvirt*(nvirt*i8)
        nsize2 = nocc*(noccEOS*i8)*nvirt*(nvirt*i8)
        IF(use_bg_buf)THEN
           call mem_pseudo_alloc(tocc2,nsize1)
           call mem_pseudo_alloc(tocc,nsize2)
        ELSE
           call mem_alloc(tocc2,nsize1)
           call mem_alloc(tocc,nsize2)
        ENDIF
!$acc enter data create(tocc) copyin(Calpha)

        call RIMP2_calc_toccA(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,nvirt,offsetV)

#ifdef VAR_OPENACC
        gpuflops = NBA*nocc*nocc*nvirt*nvirt + nocc*nocc*nvirt*nvirt*noccEOS
        call AddFLOP_FLOPonGPUaccouting(gpuflops)
#endif
        !Transform second occupied index (IDIAG,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BDIAG)
        M = noccEOS              !rows of Output Matrix
        N = noccEOS*nvirt*nvirt  !columns of Output Matrix
        K = nocc                 !summation dimension
!$acc enter data create(tocc2)
#ifdef VAR_OPENACC
!$acc host_data use_device(tocc,UoccEOST,tocc2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
        call dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2,M)
#elif defined(VAR_CUBLAS)
        stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
             & 1.0E0_realk,c_loc(UoccEOST),int(K,kind=4),c_loc(tocc),int(K,kind=4),&
             & 0.0E0_realk,c_loc(tocc2),int(M,kind=4))
#endif
!$acc end host_data
!$acc exit data delete(tocc)
        call addDGEMM_FLOPonGPUaccouting(M,N,K,0.0E0_realk)
#else
        call dgemm('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2,M)
#endif
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(tocc)
        ELSE
           call mem_dealloc(tocc)
        ENDIF
     ENDIF
     !Transform first Virtual index (ILOC,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BLOC)
     M = noccEOS*noccEOS*nvirt  !rows of Output Matrix
     N = nvirt                  !columns of Output Matrix
     K = nvirt                  !summation dimension
     nsize = nvirt*(nvirt*i8)*noccEOS*(i8*noccEOS)
     IF(use_bg_buf)THEN
        call mem_pseudo_alloc(tocc3,nsize)
     ELSE
        call mem_alloc(tocc3,nsize)
     ENDIF
!$acc enter data create(tocc3)
#ifdef VAR_OPENACC
!$acc host_data use_device(tocc2,UvirtT,tocc3)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
     call dgemm_acc('N','N',M,N,K,1.0E0_realk,tocc2,M,UvirtT,K,0.0E0_realk,tocc3,M)
#elif defined(VAR_CUBLAS)
     stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                           & 1.0E0_realk,c_loc(tocc2),int(M,kind=4),c_loc(UvirtT),int(K,kind=4),&
                           & 0.0E0_realk,c_loc(tocc3),int(M,kind=4))
#endif
!$acc end host_data
!$acc exit data delete(tocc2)
     call addDGEMM_FLOPonGPUaccouting(M,N,K,0.0E0_realk)
#else
     call dgemm('N','N',M,N,K,1.0E0_realk,tocc2,M,UvirtT,K,0.0E0_realk,tocc3,M)
#endif     
     !Final virtual transformation and reorder to dimocc
     call tensor_ainit(toccEOS,dimocc,4)
!$acc enter data create(toccEOS%elm1)
     call RIMP2_calc_toccB(nvirt,noccEOS,tocc3,UvirtT,toccEOS%elm1)
!$acc exit data copyout(toccEOS%elm1) async(async_id(2))
!$acc exit data delete(tocc3)
#ifdef VAR_OPENACC
     gpuflops = noccEOS*noccEOS*nvirt*nvirt*nvirt
     call AddFLOP_FLOPonGPUaccouting(gpuflops)
#endif
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(tocc3)     
        call mem_pseudo_dealloc(tocc2)
     ELSE
        call mem_dealloc(tocc3)     
        call mem_dealloc(tocc2)
     ENDIF
     CALL LSTIMER('RIMP2: toccEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     call tensor_ainit(toccEOS,dimocc,4)
     nsize = nvirt*noccEOS*nvirt*noccEOS
     call ls_dzero8(toccEOS%elm1,nsize)
  ENDIF
  CALL LSTIMER('DECRIMP2: tocc          ',TS2,TE2,LUPRI,FORCEPRINT)

#if defined(VAR_OPENACC) && defined(VAR_CUDA)
  !In case of GPU usage tvirt must also fit on device memory
  call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
  if(DECinfo%PL>2)then
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

  !=====================================================================================
  !  Major Step 6: Generate tvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
  !=====================================================================================
  dimvirt = [nvirtEOS,nocc,nvirtEOS,nocc]   ! Output order
  IF(NBA.GT.0)THEN
     !Calculate and partial transform to local basis - transform occupied indices
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
     nsize1 = nocc*nocc*(nvirtEOS*i8)*nvirt
     nsize2 = nocc*nocc*(nvirtEOS*i8)*nvirtEOS
     IF(use_bg_buf)THEN
        call mem_pseudo_alloc(tvirt2,nsize2)
        call mem_pseudo_alloc(tvirt,nsize1) !IDIAG,JDIAG,ALOC,BDIAG        
     ELSE
        call mem_alloc(tvirt2,nsize2)
        call mem_alloc(tvirt,nsize1) !IDIAG,JDIAG,ALOC,BDIAG
     ENDIF
!$acc enter data create(tvirt)
     call RIMP2_calc_tvirtA(nvirt,nocc,nvirtEOS,NBA,Calpha,EVocc,EVvirt,tvirt,UvirtEOST)
!$acc exit data delete(EVocc,EVvirt)
#ifdef VAR_OPENACC
     gpuflops = NBA*nocc*nocc*nvirt*nvirt + nocc*nocc*nvirt*nvirt*nvirtEOS
     call AddFLOP_FLOPonGPUaccouting(gpuflops)
#endif
     call mem_dealloc(EVocc)
     call mem_dealloc(EVvirt)

     !Transform first Virtual index (IDIAG,JDIAG,ALOC,BDIAG) => (IDIAG,JDIAG,ALOC,BLOC)
     M = nocc*nocc*nvirtEOS     !rows of Output Matrix
     N = nvirtEOS               !columns of Output Matrix
     K = nvirt                  !summation dimension
!$acc enter data create(tvirt2)
#ifdef VAR_OPENACC
!$acc host_data use_device(tvirt,UvirtEOST,tvirt2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
     call dgemm_acc('N','N',M,N,K,1.0E0_realk,tvirt,M,UvirtEOST,K,0.0E0_realk,tvirt2,M)
#elif defined(VAR_CUBLAS)
     stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                           & 1.0E0_realk,c_loc(tvirt),int(M,kind=4),c_loc(UvirtEOST),int(K,kind=4),&
                           & 0.0E0_realk,c_loc(tvirt2),int(M,kind=4))
#endif
!$acc end host_data
!$acc exit data delete(tvirt)
     call addDGEMM_FLOPonGPUaccouting(M,N,K,0.0E0_realk)
#else
     call dgemm('N','N',M,N,K,1.0E0_realk,tvirt,M,UvirtEOST,K,0.0E0_realk,tvirt2,M)
#endif
     nsize = nocc*nocc*(nvirtEOS*i8)*nvirtEOS
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(tvirt)
        call mem_pseudo_alloc(tvirt3,nsize)
     ELSE
        call mem_dealloc(tvirt)
        call mem_alloc(tvirt3,nsize)
     ENDIF
     !Transform first occupied index (IDIAG,JDIAG,ALOC,BLOC) => (ILOC,JDIAG,ALOC,BLOC)
     M = nocc                    !rows of Output Matrix
     N = nocc*nvirtEOS*nvirtEOS  !columns of Output Matrix
     K = nocc                    !summation dimension
!$acc enter data create(tvirt3)
#ifdef VAR_OPENACC
!$acc host_data use_device(tvirt2,UoccT,tvirt3)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
     call dgemm_acc('T','N',M,N,K,1.0E0_realk,UoccT,K,tvirt2,M,0.0E0_realk,tvirt3,M)
#elif defined(VAR_CUBLAS)
     stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                           & 1.0E0_realk,c_loc(UoccT),int(K,kind=4),c_loc(tvirt2),int(M,kind=4),&
                           & 0.0E0_realk,c_loc(tvirt3),int(M,kind=4))
#endif
!$acc end host_data
!$acc exit data delete(tvirt2)
     call addDGEMM_FLOPonGPUaccouting(M,N,K,0.0E0_realk)
#else
     call dgemm('T','N',M,N,K,1.0E0_realk,UoccT,K,tvirt2,M,0.0E0_realk,tvirt3,M)
#endif

     !transform last occ index to local basis and reorder 
     call tensor_ainit(tvirtEOS,dimvirt,4)
!$acc enter data create(tvirtEOS%elm1)
     call RIMP2_calc_tvirtB(nvirtEOS,nocc,tvirt3,UoccT,tvirtEOS%elm1)
!$acc exit data copyout(tvirtEOS%elm1) async(async_id(3))
!$acc exit data delete(tvirt3)
#ifdef VAR_OPENACC
     gpuflops = nvirtEOS*nocc*nvirtEOS*nocc*nvirt
     call AddFLOP_FLOPonGPUaccouting(gpuflops)
#endif
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(tvirt3)
        call mem_pseudo_dealloc(tvirt2)
     ELSE
        call mem_dealloc(tvirt3)
        call mem_dealloc(tvirt2)
     ENDIF
     CALL LSTIMER('RIMP2: tvirtEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     call tensor_ainit(tvirtEOS,dimvirt,4)
     call mem_dealloc(EVocc)
     call mem_dealloc(EVvirt)
     nSize = nvirtEOS*nocc*nvirtEOS*nocc
     call ls_dzero8(tvirtEOS%elm1,nSize)
  ENDIF

#if defined(VAR_OPENACC) && defined(VAR_CUDA)
  call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
  if(DECinfo%PL>2)then
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

  !=====================================================================================
  !  Major Step 7: Generate goccEOS(nvirt,noccEOS,nvirt,noccEOS)
  !=====================================================================================

  IF(NBA.GT.0)THEN
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
     nsize = nba*nvirt*noccEOS
     IF(use_bg_buf)THEN
        call mem_pseudo_alloc(Calpha2,nsize)
        call mem_pseudo_alloc(Calpha3,nsize)
     ELSE
        call mem_alloc(Calpha2,nsize)
        call mem_alloc(Calpha3,nsize)
     ENDIF
     ! Transform Calpha(ALPHA,a,i) to local occupied index and local Virt
     ! Transform index delta to local occupied index 
     !(alphaAux;gamma,Jloc) = (alphaAux;gamma,J)*U(J,Jloc)     UoccEOST(iDIAG,iLOC)
     M = nba*nvirt  !rows of Output Matrix
     N = noccEOS          !columns of Output Matrix
     K = nocc             !summation dimension
!$acc enter data create(Calpha2)

#ifdef VAR_OPENACC
!$acc host_data use_device(Calpha,UoccEOST,Calpha2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
     call dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccEOST,K,0.0E0_realk,Calpha2,M)
#elif defined(VAR_CUBLAS)
     stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                           & 1.0E0_realk,c_loc(Calpha),int(M,kind=4),c_loc(UoccEOST),int(K,kind=4),&
                           & 0.0E0_realk,c_loc(Calpha2),int(M,kind=4))
#endif
!$acc end host_data
!$acc exit data delete(UoccEOST) if(.not. first_order)
     call addDGEMM_FLOPonGPUaccouting(M,N,K,0.0E0_realk)
#else
     call dgemm('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccEOST,K,0.0E0_realk,Calpha2,M)
#endif

     IF(.NOT.first_order)call mem_dealloc(UoccEOST)

!$acc enter data create(Calpha3)
     call RIMP2_TransAlpha1(nvirt,noccEOS,nba,UvirtT,Calpha2,Calpha3)
!$acc exit data delete(Calpha2,UvirtT) if(.not. first_order)
!$acc exit data delete(Calpha2) if(first_order)
#ifdef VAR_OPENACC
     gpuflops = NBA*nvirt*nvirt*noccEOS
     call AddFLOP_FLOPonGPUaccouting(gpuflops)
#endif
     IF(.NOT.first_order)call mem_dealloc(UvirtT)
     
     call tensor_ainit(goccEOS,dimocc,4)
#ifdef VAR_OPENACC
!$acc enter data create(goccEOS%elm1)
     call RIMP2_calc_gocc(nvirt,noccEOS,NBA,Calpha3,goccEOS%elm1)
!$acc exit data copyout(goccEOS%elm1) async(async_id(4))
!$acc exit data delete(Calpha3)
#else
     !goccEOS(nvirt,noccEOS,nvirt,noccEOS)
     M = nvirt*noccEOS  !rows of Output Matrix
     N = nvirt*noccEOS  !columns of Output Matrix
     K = NBA            !summation dimension
     call dgemm('T','N',M,N,K,1.0E0_realk,Calpha3,K,Calpha3,K,0.0E0_realk,goccEOS%elm1,M)
#endif
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(Calpha3)
        call mem_pseudo_dealloc(Calpha2)
     ELSE
        call mem_dealloc(Calpha3)
        call mem_dealloc(Calpha2)
     ENDIF
     CALL LSTIMER('RIMP2: goccEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     IF(.NOT.first_order)call mem_dealloc(UoccEOST)
     IF(.NOT.first_order)call mem_dealloc(UvirtT)
     call tensor_ainit(goccEOS,dimocc,4)
     nSize = nvirt*noccEOS*nvirt*noccEOS
     call ls_dzero8(goccEOS%elm1,nsize)
  ENDIF

  !=====================================================================================
  !  Major Step 8: Generate gvirtEOS(nvirtEOS,nocc,nvirtEOS,nocctot)
  !=====================================================================================
#if defined(VAR_OPENACC) && defined(VAR_CUDA)
  call get_dev_mem(total_gpu,free_gpu) !free_gpu is amount of free memory in BYTES     
  if(DECinfo%PL>2)then
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

  dimvirt = [nvirtEOS,nocc,nvirtEOS,nocctot]   ! Output order
  IF(NBA.GT.0)THEN
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
     IF(fc)THEN
        nsize2 = nba*nvirt*nocctot
        nsize3 = nba*nvirt*nocctot
        IF(use_bg_buf)THEN
           call mem_pseudo_alloc(Calpha2,nsize2)
           call mem_pseudo_alloc(Calpha3,nsize3)
        ELSE
           call mem_alloc(Calpha2,nsize2)
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

#ifdef VAR_OPENACC
!$acc host_data use_device(Calpha3,UoccallT,Calpha2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
        call dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha3,M,UoccallT,K,0.0E0_realk,Calpha2,M)
#elif defined(VAR_CUBLAS)
        stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
             & 1.0E0_realk,c_loc(Calpha3),int(M,kind=4),c_loc(UoccallT),int(K,kind=4),&
             & 0.0E0_realk,c_loc(Calpha2),int(M,kind=4))
#endif
!$acc end host_data
!$acc exit data delete(UoccallT,Calpha3) if(fc .and. (.not. first_order))
!$acc exit data delete(Calpha3) if(fc .and. first_order)
        call addDGEMM_FLOPonGPUaccouting(M,N,K,0.0E0_realk)
#else
        call dgemm('N','N',M,N,K,1.0E0_realk,Calpha3,M,UoccallT,K,0.0E0_realk,Calpha2,M)
#endif        
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(Calpha3)
        ELSE
           call mem_dealloc(Calpha3)
        ENDIF
        IF(.NOT.first_order)call mem_dealloc(UoccallT)
     ELSE
        nsize = nba*nvirt*nocc
        IF(use_bg_buf)THEN
           call mem_pseudo_alloc(Calpha2,nsize)
        ELSE
           call mem_alloc(Calpha2,nsize)
        ENDIF
        ! Transform index delta to local occupied index 
        !(alphaAux;gamma,Jloc) = (alphaAux;gamma,J)*U(J,Jloc)     UoccEOST(iDIAG,iLOC)
        M = nba*nvirt        !rows of Output Matrix
        N = nocc             !columns of Output Matrix
        K = nocc             !summation dimension
!$acc enter data create(Calpha2)
#ifdef VAR_OPENACC
!$acc host_data use_device(Calpha,UoccT,Calpha2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
        call dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccT,K,0.0E0_realk,Calpha2,M)
#elif defined(VAR_CUBLAS)
        stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
             & 1.0E0_realk,c_loc(Calpha),int(M,kind=4),c_loc(UoccT),int(K,kind=4),&
             & 0.0E0_realk,c_loc(Calpha2),int(M,kind=4))
#endif
!$acc end host_data
!$acc exit data delete(UoccT,Calpha) if(.not. first_order)
        call addDGEMM_FLOPonGPUaccouting(M,N,K,0.0E0_realk)
#else
        call dgemm('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccT,K,0.0E0_realk,Calpha2,M)
#endif            
     ENDIF
     IF(.NOT.first_order)call mem_dealloc(UoccT)
     nsize = nba*nvirtEOS*nocctot
     IF(use_bg_buf)THEN
        call mem_pseudo_alloc(Calpha3,nsize)
     ELSE
        call mem_alloc(Calpha3,nsize)
     ENDIF
!$acc enter data create(Calpha3)
     call RIMP2_TransAlpha2(nocctot,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha3)
!$acc exit data delete(Calpha2,UvirtEOST) if(.not. first_order)
!$acc exit data delete(Calpha2) if(first_order)
#ifdef VAR_OPENACC
     gpuflops = NBA*nvirtEOS*nocctot*nvirt
     call AddFLOP_FLOPonGPUaccouting(gpuflops)
#endif
     IF(.NOT.first_order)call mem_dealloc(UvirtEOST)

     call tensor_ainit(gvirtEOS,dimvirt,4)
!$acc enter data create(gvirtEOS%elm1)
     call RIMP2_calc_gvirt(nvirtEOS,nocctot,NBA,nocc,Calpha3,gvirtEOS%elm1,offset)
!$acc exit data copyout(gvirtEOS%elm1) async(async_id(5))
!$acc exit data delete(Calpha3)

     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(Calpha3)
        call mem_pseudo_dealloc(Calpha2)
     ELSE
        call mem_dealloc(Calpha3)
        call mem_dealloc(Calpha2)
     ENDIF
     IF(.NOT.first_order)THEN
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(Calpha)
        ELSE
           call mem_dealloc(Calpha)
        ENDIF
     ENDIF

     CALL LSTIMER('RIMP2: gvirtEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     IF(.NOT.first_order)call mem_dealloc(UvirtEOST)
     IF(.NOT.first_order)call mem_dealloc(UoccT)
     call tensor_ainit(gvirtEOS,dimvirt,4)
     nSize = nvirtEOS*nocc*nvirtEOS*nocctot
     call ls_dzero8(gvirtEOS%elm1,nsize)
  ENDIF

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
     nSize = nvirt*noccEOS*nvirt*noccEOS
     call lsmpi_reduction(toccEOS%elm1,nsize,infpar%master,infpar%lg_comm)
     nSize = nvirtEOS*nocc*nvirtEOS*nocc
     call lsmpi_reduction(tvirtEOS%elm1,nSize,infpar%master,infpar%lg_comm)
     nSize = nvirt*noccEOS*nvirt*noccEOS
     call lsmpi_reduction(goccEOS%elm1,nSize,infpar%master,infpar%lg_comm)
     nSize = nvirtEOS*nocc*nvirtEOS*nocctot
     call lsmpi_reduction(gvirtEOS%elm1,nsize,infpar%master,infpar%lg_comm)
     call time_start_phase(PHASE_WORK)
     if (.not. Master) then
        call tensor_free(toccEOS)
        call tensor_free(tvirtEOS)
        call tensor_free(goccEOS)
        call tensor_free(gvirtEOS)
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
        nsize3 = nba*nvirt*noccEOS
        IF(use_bg_buf)THEN
           call mem_pseudo_alloc(Calpha2,nsize2)
           call mem_pseudo_alloc(Calpha3,nsize3)
        ELSE
           call mem_alloc(Calpha2,nsize2)
           call mem_alloc(Calpha3,nsize3)
        ENDIF

        CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
        !(alphaAux;nvirt,JnoccEOS) = (alphaAux;nvirt,J)*U(J,JnoccEOS)
        M = nba*nvirt        !rows of Output Matrix
        N = noccEOS          !columns of Output Matrix
        K = nocc             !summation dimension
!$acc enter data create(Calpha2)
#ifdef VAR_OPENACC
!$acc host_data use_device(Calpha,UoccEOST,Calpha2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
        call dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccEOST,K,0.0E0_realk,Calpha2,M)
#elif defined(VAR_CUBLAS)
        stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                              & 1.0E0_realk,c_loc(Calpha),int(M,kind=4),c_loc(UoccEOST),int(K,kind=4),&
                              & 0.0E0_realk,c_loc(Calpha2),int(M,kind=4))
#endif
!$acc end host_data
#else
        call dgemm('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccEOST,K,0.0E0_realk,Calpha2,M)
#endif
        
        !(alphaAux,nvirtAOS,noccEOS) = (alphaAux;nvirt,noccEOS)*Uvirt(nvirt,nvirtAOS)
!$acc enter data create(Calpha3)
        call RIMP2_TransAlpha2(noccEOS,nvirt,nvirt,nba,UvirtT,Calpha2,Calpha3)
        
        CALL LSTIMER('START ',TS2,TE2,LUPRI)
        IF(DECinfo%frozencore)THEN
!           call Build_CalphaMO(MyFragment%MyLsitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
!                & CollaborateWithSlaves,CDIAGocc%val,nocc,CDIAGoccALL,nocctot,mynum,&
!                & numnodes,nAtomsAux,CalphaOcc,NBA,ABdecomp,ABdecompCreate)
           intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
           intspec(2) = 'R' !Regular AO basis function on center 3
           intspec(3) = 'R' !Regular AO basis function on center 4
           intspec(4) = 'C' !Coulomb Operator
           call Build_CalphaMO2(MyFragment%mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
                & FORCEPRINT,CollaborateWithSlaves,CDIAGocc%val,nocc,CDIAGoccALL,nocctot,&
                & mynum,numnodes,CalphaOcc,NBA,ABdecomp,ABdecompCreate,intspec,use_bg_buf)
           call mem_dealloc(CDIAGoccALL)
        ELSE
!           call Build_CalphaMO(MyFragment%MyLsitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
!                & CollaborateWithSlaves,CDIAGocc%val,nocc,CDIAGocc%val,nocc,mynum,&
!                & numnodes,nAtomsAux,CalphaOcc,NBA,ABdecomp,ABdecompCreate)
           intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
           intspec(2) = 'R' !Regular AO basis function on center 3
           intspec(3) = 'R' !Regular AO basis function on center 4
           intspec(4) = 'C' !Coulomb Operator
           call Build_CalphaMO2(MyFragment%mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
                & FORCEPRINT,CollaborateWithSlaves,CDIAGocc%val,nocc,CDIAGocc%val,nocc,&
                & mynum,numnodes,CalphaOcc,NBA,ABdecomp,ABdecompCreate,intspec,use_bg_buf)
           IF(nocctot.NE.nocc)call lsquit('FC Error RIMP2.',-1)
        ENDIF
        CALL LSTIMER('DECRIMP2: CalphaOO',TS2,TE2,LUPRI,FORCEPRINT)
        call array2_free(CDIAGocc)

        !(alphaAux;nocc,noccAOS=nocctot) = (alphaAux;nocc,nocc)*UoccallT(nocctot,nocctot)
        M = nba*nocc         !rows of Output Matrix
        N = nocctot          !columns of Output Matrix
        K = nocctot          !summation dimension
!$acc enter data copyin(CalphaOcc)
        IF(DECinfo%frozencore)THEN
#ifdef VAR_OPENACC
!$acc host_data use_device(CalphaOcc,UoccallT,Calpha2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
           call dgemm_acc('N','N',M,N,K,1.0E0_realk,CalphaOcc,M,UoccallT,K,0.0E0_realk,Calpha2,M)
#elif defined(VAR_CUBLAS)
           stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                                 & 1.0E0_realk,c_loc(CalphaOcc),int(M,kind=4),c_loc(UoccallT),int(K,kind=4),&
                                 & 0.0E0_realk,c_loc(Calpha2),int(M,kind=4))
#endif
!$acc end host_data
#else
           call dgemm('N','N',M,N,K,1.0E0_realk,CalphaOcc,M,UoccallT,K,0.0E0_realk,Calpha2,M)
#endif
!$acc exit data delete(UoccallT)
           call mem_dealloc(UoccallT)
        ELSE
#ifdef VAR_OPENACC
!$acc host_data use_device(CalphaOcc,UoccT,Calpha2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
           call dgemm_acc('N','N',M,N,K,1.0E0_realk,CalphaOcc,M,UoccT,K,0.0E0_realk,Calpha2,M)
#elif defined(VAR_CUBLAS)
           stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                                 & 1.0E0_realk,c_loc(CalphaOcc),int(M,kind=4),c_loc(UoccT),int(K,kind=4),&
                                 & 0.0E0_realk,c_loc(Calpha2),int(M,kind=4))
#endif
!$acc end host_data
#else
           call dgemm('N','N',M,N,K,1.0E0_realk,CalphaOcc,M,UoccT,K,0.0E0_realk,Calpha2,M)
#endif
        ENDIF
!$acc exit data delete(CalphaOcc)
        nsize = nba*noccEOS*nocctot
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(CalphaOcc)
           call mem_pseudo_alloc(Calpha4,nsize)
        ELSE
           call mem_dealloc(CalphaOcc)
           call mem_alloc(Calpha4,nsize)
        ENDIF
        !(alphaAux,noccEOS,noccAOS=nocctot) = (alphaAux;nocc,noccAOS=nocctot)*UoccEOST(nocc,noccEOS)
!$acc enter data create(Calpha4)
        call RIMP2_TransAlpha2(nocctot,nocc,noccEOS,nba,UoccEOST,Calpha2,Calpha4)
!$acc exit data delete(UoccEOST,Calpha2)
        call mem_dealloc(UoccEOST)
        
        !  djikEOS(nvirtAOS,noccEOS,noccEOS,noccAOS)
        call tensor_ainit(djik,dimvirt,4)
!$acc enter data create(djik%elm1)
        call RIMP2_calc_gen4DimFO(NBA,Calpha3,nvirt,noccEOS,Calpha4,noccEOS,nocctot,djik%elm1)
!$acc exit data delete(Calpha3,Calpha4)
!$acc exit data copyout(djik%elm1) async(async_id(6))
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(Calpha4)
           call mem_pseudo_dealloc(Calpha3)
           call mem_pseudo_dealloc(Calpha2)
        ELSE
           call mem_dealloc(Calpha4)
           call mem_dealloc(Calpha3)
           call mem_dealloc(Calpha2)
        ENDIF
        CALL LSTIMER('RIMP2: djik',TS3,TE3,LUPRI,FORCEPRINT)
     ELSE
        call array2_free(CDIAGocc)
        call mem_dealloc(UoccEOST)
        call tensor_ainit(djik,dimvirt,4)
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
           call mem_pseudo_alloc(Calpha2,nsize)
        ELSE
           call mem_alloc(Calpha2,nsize)
        ENDIF
        !(alphaAux;nvirt,noccAOS) = (alphaAux;nvirt,nocc)*U(nocc,noccAOS)
        M = nba*nvirt        !rows of Output Matrix
        N = nocc             !columns of Output Matrix
        K = nocc             !summation dimension
!$acc enter data create(Calpha2)
#ifdef VAR_OPENACC
!$acc host_data use_device(Calpha,UoccT,Calpha2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
        call dgemm_acc('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccT,K,0.0E0_realk,Calpha2,M)
#elif defined(VAR_CUBLAS)
        stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                              & 1.0E0_realk,c_loc(Calpha),int(M,kind=4),c_loc(UoccT),int(K,kind=4),&
                              & 0.0E0_realk,c_loc(Calpha2),int(M,kind=4))
#endif
!$acc end host_data
#else
        call dgemm('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccT,K,0.0E0_realk,Calpha2,M)
#endif
!$acc exit data delete(UoccT) 
        call mem_dealloc(UoccT)
        
        !(alphaAux,nvirtEOS,noccAOS) = (alphaAux;nvirt,noccAOS)*Uvirt(nvirt,nvirtEOS)
        nsize = nba*nvirtEOS*nocc
        IF(use_bg_buf)THEN
           call mem_pseudo_alloc(Calpha3,nsize)
        ELSE
           call mem_alloc(Calpha3,nsize)
        ENDIF
!$acc enter data create(Calpha3)
        call RIMP2_TransAlpha2(nocc,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha3)   
        !(alphaAux;nvirt,nvirtAOS) = (alphaAux;nvirt,nvirt)*UvirtT(nvirt,nvirt)
        M = nba*nvirt        !rows of Output Matrix
        N = nvirt            !columns of Output Matrix
        K = nocc             !summation dimension

        CALL LSTIMER('START ',TS2,TE2,LUPRI)
        intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
        intspec(2) = 'R' !Regular AO basis function on center 3
        intspec(3) = 'R' !Regular AO basis function on center 4
        intspec(4) = 'C' !Coulomb Operator
        call Build_CalphaMO2(MyFragment%mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
             & FORCEPRINT,CollaborateWithSlaves,CDIAGvirt%val,nvirt,CDIAGvirt%val,nvirt,&
             & mynum,numnodes,CalphaVV,NBA,ABdecomp,ABdecompCreate,intspec,use_bg_buf)
        call mem_dealloc(ABdecomp)
        CALL LSTIMER('DECRIMP2: CalphaVV',TS2,TE2,LUPRI,FORCEPRINT)
        call array2_free(CDIAGvirt)

!$acc enter data copyin(CalphaVV)
#ifdef VAR_OPENACC
!$acc host_data use_device(CalphaVV,UvirtT,Calpha2)
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
        call dgemm_acc('N','N',M,N,K,1.0E0_realk,CalphaVV,M,UvirtT,K,0.0E0_realk,Calpha2,M)
#elif defined(VAR_CUBLAS)
        stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(M,kind=4),int(N,kind=4),int(K,kind=4),&
                              & 1.0E0_realk,c_loc(CalphaVV),int(M,kind=4),c_loc(UvirtT),int(K,kind=4),&
                              & 0.0E0_realk,c_loc(Calpha2),int(M,kind=4))
#endif
!$acc end host_data
#else
        call dgemm('N','N',M,N,K,1.0E0_realk,CalphaVV,M,UvirtT,K,0.0E0_realk,Calpha2,M)
#endif
!$acc exit data delete(CalphaVV,UvirtT,Calpha)
        nsize = nba*nvirtEOS*nvirt
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(CalphaVV)           
           call mem_pseudo_alloc(Calpha4,nsize)
        ELSE
           call mem_dealloc(CalphaVV)
           call mem_alloc(Calpha4,nsize)
        ENDIF
        call mem_dealloc(UvirtT)
        !(alphaAux,nvirtEOS,nvirtAOS) = (alphaAux;nvirt,nvirtAOS)*UvirtEOST(nvirt,nvirtEOS)
!$acc enter data create(Calpha4)
        call RIMP2_TransAlpha2(nvirt,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha4)
!$acc exit data delete(UvirtEOST,Calpha2)
        call mem_dealloc(UvirtEOST)
        
        !generate blad(nvirtEOS,noccAOS,nvirtEOS,nvirtAOS)
        call tensor_ainit(blad,dimvirt,4)
!$acc enter data create(blad%elm1)
        call RIMP2_calc_gen4DimFO(NBA,Calpha3,nvirtEOS,nocc,Calpha4,nvirtEOS,nvirt,blad%elm1)
!$acc exit data delete(Calpha3,Calpha4)
!$acc exit data copyout(blad%elm1) async(async_id(7))
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(Calpha4)
           call mem_pseudo_dealloc(Calpha3)
           call mem_pseudo_dealloc(Calpha2)
           call mem_pseudo_dealloc(Calpha)
        ELSE
           call mem_dealloc(Calpha4)
           call mem_dealloc(Calpha3)
           call mem_dealloc(Calpha2)
           call mem_dealloc(Calpha)
        ENDIF
        CALL LSTIMER('RIMP2: blad',TS3,TE3,LUPRI,FORCEPRINT)
     ELSE
        call array2_free(CDIAGvirt)
        call mem_dealloc(UoccT)
        call mem_dealloc(UvirtT)
        call mem_dealloc(UvirtEOST)
        call tensor_ainit(blad,dimvirt,4)
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

  if(DECinfo%PL>0)THEN
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
#endif

end subroutine RIMP2_integrals_and_amplitudes

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

subroutine Build_CalphaMO(myLSitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
     & CollaborateWithSlaves,Cvirt,nvirt,Cocc,nocc,mynum,numnodes,nAtomsAux,Calpha,&
     & NBA,AlphaBetaDecomp,AlphaBetaDecompCreate,Oper)
  implicit none
  type(lsitem), intent(inout) :: mylsitem
  integer,intent(inout) :: NBA
  integer,intent(in) :: nAtomsAux,nocc,nvirt
  integer,intent(in) :: nbasisAux,LUPRI,nbasis,mynum,numnodes
  logical,intent(in) :: master,FORCEPRINT,CollaborateWithSlaves,AlphaBetaDecompCreate
  real(realk) :: AlphaBetaDecomp(nbasisAux,nbasisAux)
  real(realk),intent(in) :: Cocc(nbasis,nocc),Cvirt(nbasis,nvirt)
  real(realk),pointer :: Calpha(:,:,:)
  integer,optional :: Oper
  !
  integer :: MynbasisAuxMPI
  integer,pointer :: nbasisAuxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
  real(realk),pointer :: AlphaBeta(:,:)
  real(realk),pointer :: TMPAlphaBetaDecomp(:,:),AlphaCD3(:,:,:)
  real(realk),pointer :: AlphaCD5(:,:,:),AlphaCDFull(:,:,:)
  real(realk) :: TS3,TE3,MemInGBCollected,SizeCalpha
  TYPE(MoleculeInfo),pointer      :: molecule1,molecule2,molecule3,molecule4
  integer(kind=long)    :: maxsize,nSize,n8
  integer(kind=ls_mpik) :: node 
  integer :: CurrentWait(2),nAwaitDealloc,iAwaitDealloc,MynAtomsMPI
  integer :: myOriginalRank,OriginalRanknbasisAuxMPI,M,N,K,I,offset,offset2
  integer :: ndimMax,nbasisAuxMPI2(numnodes),MynbasisAuxMPI2,rimp2_nodtot
  integer :: nbuf1,nbuf2,nbuf3,inode,J
  logical :: useAlphaCD5,useAlphaCD6,ChangedDefault,first_order,MessageRecieved
  logical :: PerformReduction,RIMPSubGroupCreated,UseSubGroupCommunicator
  PerformReduction = .TRUE.
  MynbasisAuxMPI2 = 0
  NBA = 0 
  call get_currently_available_memory(MemInGBCollected)
  
  !===========================================================
  !   Determine Scheme to Use (AllReduce, Bcast Method)
  !===========================================================
  CALL LSTIMER('START ',TS3,TE3,LUPRI)

  IF(master)THEN
     IF(DECinfo%RIMP2ForcePDMCalpha)THEN
        PerformReduction = .FALSE.
        IF(numnodes.EQ.1)THEN           
           PerformReduction = .TRUE.
        ENDIF
     ELSE
        !maxsize = max number of floating point elements
        SizeCalpha = (nbasisAux+nbasisAux/numnodes)*nocc*nvirt*8E-9_realk
        IF(SizeCalpha.LT.MemInGBCollected*0.75E0_realk.OR.numnodes.EQ.1)THEN
           !Calpha can fit on all nodes Which means we can do a reduction.
           PerformReduction = .TRUE.
           WRITE(DECinfo%output,'(A,F8.1,A)')'RIMP2: Full (alpha|cd) integral requires ',SizeCalpha,' GB'
           WRITE(DECinfo%output,'(A,F8.1,A)')'RIMP2: Memory available                  ',MemInGBCollected,' GB'
        ELSE
           !Calpha cannot fit so we distribute this - which means more MPI communication.
           PerformReduction = .FALSE.
           !Determine number of nodes to use to construct Calpha. (use same to determine the MPI split)
        ENDIF
     ENDIF
#ifdef VAR_MPI
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     call ls_mpibcast(PerformReduction,infpar%master,infpar%lg_comm)
     call time_start_phase( PHASE_WORK )
#endif
  ELSE
#ifdef VAR_MPI
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     call ls_mpibcast(PerformReduction,infpar%master,infpar%lg_comm)
     call time_start_phase( PHASE_WORK )
#endif
     IF(PerformReduction)THEN
        !the master have enough space 
        !maxsize = max number of floating point elements
        SizeCalpha = nbasisAux*nocc*nvirt*8E-9_realk
        IF(SizeCalpha.GT.MemInGBCollected*0.6E0_realk)THEN
           print*,'WARNING: Master have space for (alpha|cd) but slave do not'
           print*,'WARNING: Full (alpha|cd) integral requires ',SizeCalpha,' GB'
           print*,'WARNING: Memory available                  ',MemInGBCollected,' GB'
        ENDIF
     ENDIF
  ENDIF

  CALL LSTIMER('Calpha1',TS3,TE3,LUPRI)

  !===========================================================
  !   Determine Sizes1: used to calc 3 center integrals
  !===========================================================

  IF(CollaborateWithSlaves)then 
     !all nodes have info about all nodes 
     call mem_alloc(nbasisAuxMPI,numnodes)           !number of Aux basis func assigned to rank
     nbasisAuxMPI = 0 
     call mem_alloc(nAtomsMPI,numnodes)              !atoms assign to rank
     call mem_alloc(startAuxMPI,nAtomsAux,numnodes)  !startindex in full (nbasisAux)
     call mem_alloc(AtomsMPI,nAtomsAux,numnodes)     !identity of atoms in full molecule
     call mem_alloc(nAuxMPI,nAtomsAux,numnodes)      !nauxBasis functions for each of the nAtomsMPI
     IF(DECinfo%AuxAtomicExtent)THEN   
        call getRIbasisMPI(mylsitem%INPUT%AUXMOLECULE,nAtomsAux,numnodes,&
             & nbasisAuxMPI,startAuxMPI,AtomsMPI,nAtomsMPI,nAuxMPI)
     ELSE
        call getRIbasisMPI(mylsitem%SETTING%MOLECULE(1)%p,nAtomsAux,numnodes,&
             & nbasisAuxMPI,startAuxMPI,AtomsMPI,nAtomsMPI,nAuxMPI)
     ENDIF
     MynAtomsMPI = nAtomsMPI(mynum+1)
     MynbasisAuxMPI = nbasisAuxMPI(mynum+1)
     rimp2_nodtot = 0 
     DO I = 1,numnodes
        IF(nbasisAuxMPI(I).GT.0) rimp2_nodtot = rimp2_nodtot + 1
     ENDDO

     call mem_dealloc(AtomsMPI) !not used in this subroutine 
  ELSE
     MynbasisAuxMPI = nbasisAux     
     rimp2_nodtot = numnodes
  ENDIF

  CALL LSTIMER('Calpha2',TS3,TE3,LUPRI)
  
  !=====================================================================================
  ! Master Obtains (alpha|beta) ERI in Auxiliary Basis 
  !=====================================================================================
  IF(AlphaBetaDecompCreate)THEN

     IF(master)THEN
        call mem_alloc(AlphaBeta,nbasisAux,nbasisAux)
        CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
        IF(DECinfo%AuxAtomicExtent)THEN
           molecule1 => mylsitem%SETTING%MOLECULE(1)%p
           molecule2 => mylsitem%SETTING%MOLECULE(2)%p
           molecule3 => mylsitem%SETTING%MOLECULE(3)%p
           molecule4 => mylsitem%SETTING%MOLECULE(4)%p
           mylsitem%SETTING%MOLECULE(1)%p => mylsitem%INPUT%AUXMOLECULE
           mylsitem%SETTING%MOLECULE(2)%p => mylsitem%INPUT%AUXMOLECULE
           mylsitem%SETTING%MOLECULE(3)%p => mylsitem%INPUT%AUXMOLECULE
           mylsitem%SETTING%MOLECULE(4)%p => mylsitem%INPUT%AUXMOLECULE
        ENDIF
        IF(present(Oper))THEN
           call II_get_RI_AlphaBeta_2centerInt(DECinfo%output,DECinfo%output,&
                & AlphaBeta,mylsitem%setting,nbasisAux,Oper)
        ELSE
           call II_get_RI_AlphaBeta_2centerInt(DECinfo%output,DECinfo%output,&
                & AlphaBeta,mylsitem%setting,nbasisAux)
        ENDIF
        IF(DECinfo%AuxAtomicExtent)THEN
           mylsitem%SETTING%MOLECULE(1)%p => molecule1
           mylsitem%SETTING%MOLECULE(2)%p => molecule2
           mylsitem%SETTING%MOLECULE(3)%p => molecule3
           mylsitem%SETTING%MOLECULE(4)%p => molecule4
        ENDIF
        
        CALL LSTIMER('AlphaBeta ',TS3,TE3,LUPRI,FORCEPRINT)
        ! Create the inverse square root AlphaBeta = (alpha|beta)^(-1/2)
        ! Warning the inverse is not unique so in order to make sure all slaves have the same
        ! inverse matrix we calculate it on the master a BCAST to slaves
        call lowdin_diag_S_minus_sqrt(nbasisAux, AlphaBeta,AlphaBetaDecomp, lupri)
        call mem_dealloc(AlphaBeta)
        CALL LSTIMER('AlphaBetamSq ',TS3,TE3,LUPRI,FORCEPRINT)
     ENDIF
#ifdef VAR_MPI
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     call ls_mpibcast(AlphaBetaDecomp,nbasisAux,nbasisAux,infpar%master,infpar%lg_comm)
     call time_start_phase(PHASE_WORK)   
#endif
  ENDIF

  CALL LSTIMER('Calpha3',TS3,TE3,LUPRI)

  !==================================================================
  !   Determine MynbasisAuxMPI2:  Calpha(MynbasisAuxMPI2,nvirt,nocc)
  !==================================================================

  IF(CollaborateWithSlaves)then 
     ndimMax = nbasisAux/numnodes
     do I=1,numnodes
        nbasisAuxMPI2(I) = ndimMax
     enddo
     J=2 !not add to master
     do I=1,MOD(nbasisAux,numnodes)
        nbasisAuxMPI2(J) = nbasisAuxMPI2(J) + 1
        J=J+1
     enddo
     MynbasisAuxMPI2 = nbasisAuxMPI2(mynum+1)
     call mem_alloc(TMPAlphaBetaDecomp,MynbasisAuxMPI2,nbasisAux)
     offset = mynum*ndimMax
     offset2 = numnodes*ndimMax + mynum -1 +1
     IF(MynbasisAuxMPI2.GT.ndimMax)THEN
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(nbasisAux,ndimMax,&
        !$OMP TMPAlphaBetaDecomp,AlphaBetaDecomp,offset,offset2)
        do I=1,nbasisAux
           do J=1,ndimMax
              TMPAlphaBetaDecomp(J,I) = AlphaBetaDecomp(offset+J,I)
           enddo
           TMPAlphaBetaDecomp(ndimMax+1,I) = AlphaBetaDecomp(offset2,I)
        enddo
        !$OMP END PARALLEL DO
     ELSE
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(nbasisAux,ndimMax,&
        !$OMP TMPAlphaBetaDecomp,AlphaBetaDecomp,offset)
        do I=1,nbasisAux
           do J=1,ndimMax
              TMPAlphaBetaDecomp(J,I) = AlphaBetaDecomp(offset+J,I)
           enddo
        enddo
        !$OMP END PARALLEL DO
     ENDIF
     NBA = MynbasisAuxMPI2
  ELSE
     NBA = nbasisAux
  ENDIF

  CALL LSTIMER('Calpha4',TS3,TE3,LUPRI)

  !=====================================================================================
  ! Obtain 3 center RI integrals (alpha,a,i) 
  !=====================================================================================

  IF(MynbasisAuxMPI.GT.0)THEN
     call get_currently_available_memory(MemInGBCollected)
     !maxsize = max number of floating point elements
     !allow to building of 3 center integral to use 60 procent of 
     !currently available memory
     maxsize = 0.60E0_realk*NINT(MemInGBCollected*1.E9_realk)
     !call mem_alloc(AlphaCD3,nbasisAux,nvirt,nocc)
     !It is very annoying but I allocated AlphaCD3 inside 
     !II_get_RI_AlphaCD_3centerInt2 due to memory concerns
     !This Part of the Code is MPI/OpenMP parallel and AlphaCD3 
     !will have the dimensions (MynbasisAuxMPI,nvirt,nocc) 
     !nbasisAuxMPI is nbasisAux divided out on the nodes so roughly 
     !nbasisAuxMPI = nbasisAux/numnodes
     IF(DECinfo%AuxAtomicExtent)THEN
        molecule1 => mylsitem%SETTING%MOLECULE(1)%p
        molecule2 => mylsitem%SETTING%MOLECULE(2)%p
        mylsitem%SETTING%MOLECULE(1)%p => mylsitem%INPUT%AUXMOLECULE
        mylsitem%SETTING%MOLECULE(2)%p => mylsitem%INPUT%AUXMOLECULE
     ENDIF
     IF(present(Oper))THEN
        call II_get_RI_AlphaCD_3centerInt2(DECinfo%output,DECinfo%output,&
             & AlphaCD3,mylsitem%setting,nbasisAux,nbasis,&
             & nvirt,nocc,Cvirt,Cocc,maxsize,mynum,numnodes,Oper)
     ELSE
        call II_get_RI_AlphaCD_3centerInt2(DECinfo%output,DECinfo%output,&
             & AlphaCD3,mylsitem%setting,nbasisAux,nbasis,&
             & nvirt,nocc,Cvirt,Cocc,maxsize,mynum,numnodes)
     ENDIF
     IF(DECinfo%AuxAtomicExtent)THEN
        mylsitem%SETTING%MOLECULE(1)%p => molecule1
        mylsitem%SETTING%MOLECULE(2)%p => molecule2
     ENDIF
  ENDIF

  CALL LSTIMER('Calpha5',TS3,TE3,LUPRI)

  !=====================================================================================
  ! MPI scheme:  PerformReduction  or   a Bcast Routine
  !=====================================================================================

  IF(PerformReduction)THEN

     !=====================================================================================
     ! MPI scheme:  PerformReduction
     !=====================================================================================

     WRITE(DECinfo%output,'(A)')'RIMP2 Calpha Scheme 1: Using Allreduce on (alpha|cd) integral'
     IF(CollaborateWithSlaves)then 
        call mem_alloc(alphaCDFull,nbasisAux,nvirt,nocc)
        n8 = nbasisAux*nocc*nvirt
        call ls_dzero8(alphaCDFull,n8)
        IF(MynbasisAuxMPI.GT.0)THEN
           call PlugInToalphaCDFull(mynum,nAtomsMPI,startAuxMPI,nocc,nvirt,nAuxMPI,&
                & alphaCDFull,alphaCD3,nbasisAux,MynbasisAuxMPI,numnodes,nAtomsAux)
           call mem_dealloc(alphaCD3)
        ENDIF
#ifdef VAR_MPI
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        call lsmpi_allreduce(alphaCDFull,nbasisAux,nvirt,nocc,infpar%lg_comm)
        call time_start_phase( PHASE_WORK )
#endif
        !Calpha = TMPAlphaBetaDecomp(MynbasisAuxMPI,nbasisAux)
        M =  MynbasisAuxMPI2   !rows of Output Matrix
        N =  nvirt*nocc       !columns of Output Matrix
        K =  nbasisAux        !summation dimension
        call mem_alloc(Calpha,MynbasisAuxMPI2,nvirt,nocc)
        call dgemm('N','N',M,N,K,1.0E0_realk,TMPAlphaBetaDecomp,&
             & M,alphaCDFull,K,0.0E0_realk,Calpha,M)
        call mem_dealloc(alphaCDFull)
        call mem_dealloc(TMPAlphaBetaDecomp)
     ELSE
        !Serial version
        M =  MynbasisAuxMPI   !rows of Output Matrix
        N =  nvirt*nocc       !columns of Output Matrix
        K =  nbasisAux        !summation dimension
        call mem_alloc(Calpha,MynbasisAuxMPI,nvirt,nocc)
        call dgemm('N','N',M,N,K,1.0E0_realk,AlphaBetaDecomp,&
             & M,AlphaCD3,K,0.0E0_realk,Calpha,M)
        call mem_dealloc(AlphaCD3)
     ENDIF
  ELSE

     !=====================================================================================
     ! MPI scheme:  Bcast Routine
     !=====================================================================================
     call mem_alloc(Calpha,MynbasisAuxMPI2,nvirt,nocc)
     nsize = MynbasisAuxMPI2*nvirt*nocc
     call ls_dzero8(Calpha,nsize)
     DO inode = 1,rimp2_nodtot
        IF(mynum.EQ.inode-1)THEN
           nbuf1 = nbasisAuxMPI(mynum+1)
           nbuf2 = nvirt
           nbuf3 = nocc
#ifdef VAR_MPI
           call time_start_phase( PHASE_IDLE )
           call lsmpi_barrier(infpar%lg_comm)
           call time_start_phase( PHASE_COMM )
           node = mynum
           call ls_mpibcast(AlphaCD3,nbuf1,nbuf2,nbuf3,node,infpar%lg_comm)
#endif
           call RIMP2_buildOwnCalphaFromAlphaCD(nocc,nvirt,mynum,numnodes,&
                & natomsAux,MynbasisAuxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,&
                & AlphaCD3,Calpha,TMPAlphaBetaDecomp,nbasisAux,&
                & MynbasisAuxMPI2)
           call mem_dealloc(AlphaCD3)
        ELSE
           nbuf1 = nbasisAuxMPI(inode)
           nbuf2 = nvirt
           nbuf3 = nocc
           node = inode-1
           !recieve
           call mem_alloc(AlphaCD5,nbasisAuxMPI(inode),nvirt,nocc)
#ifdef VAR_MPI
           call time_start_phase( PHASE_IDLE )
           call lsmpi_barrier(infpar%lg_comm)
           call time_start_phase( PHASE_COMM )
           call ls_mpibcast(AlphaCD5,nbuf1,nbuf2,nbuf3,node,infpar%lg_comm)
           call time_start_phase( PHASE_WORK )
#endif
           myOriginalRank = inode-1
           OriginalRanknbasisAuxMPI = nbasisAuxMPI(inode)
           call RIMP2_buildCalphaContFromAlphaCD(nocc,nvirt,myOriginalRank,&
                & numnodes,natomsAux,OriginalRanknbasisAuxMPI,&
                & nAtomsMPI,startAuxMPI,nAuxMPI,AlphaCD5,&
                & Calpha,TMPAlphaBetaDecomp,nbasisAux,MynbasisAuxMPI2)
           call mem_dealloc(AlphaCD5)
        ENDIF
     ENDDO
     call mem_dealloc(TMPAlphaBetaDecomp)
  ENDIF
  CALL LSTIMER('Calpha6',TS3,TE3,LUPRI)
  IF(CollaborateWithSlaves)then 
     call mem_dealloc(nbasisAuxMPI)
     call mem_dealloc(startAuxMPI)
     call mem_dealloc(nAtomsMPI)
     call mem_dealloc(nAuxMPI)
  ENDIF
!  call sleep(mynum*5)
!  PRINT*,'MynbasisAuxMPI2',MynbasisAuxMPI2
!  WRITE(6,*)'Final Calph(NBA=',NBA,',nvirt=',nvirt,',nocc=',nocc,')'
!  WRITE(6,*)'Print Subset Final Calph(NBA=',NBA,',1:4)  MYNUM',MYNUM
!  call ls_output(Calpha,1,NBA,1,4,NBA,nvirt*nocc,1,6)

end subroutine Build_CalphaMO

!This should be call my master and slaves
subroutine Build_CalphaMO2(myLSitem,master,nbasis1,nbasis2,nbasisAux,LUPRI,FORCEPRINT,&
     & CollaborateWithSlaves,Cvirt,nvirt,Cocc,nocc,mynum,numnodes,Calpha,&
     & NBA,AlphaBetaDecomp,AlphaBetaDecompCreate,intspec,use_bg_bufInput)
  implicit none
  type(lsitem), intent(inout) :: mylsitem
  integer,intent(inout) :: NBA
  integer,intent(in) :: nocc,nvirt
  integer,intent(in) :: nbasisAux,LUPRI,nbasis1,nbasis2,mynum,numnodes
  logical,intent(in) :: master,FORCEPRINT,CollaborateWithSlaves,AlphaBetaDecompCreate
  real(realk),intent(inout) :: AlphaBetaDecomp(nbasisAux,nbasisAux)
  real(realk),intent(in) :: Cvirt(nbasis1,nvirt),Cocc(nbasis2,nocc)
  real(realk),pointer :: Calpha(:)
  character,intent(in) :: intspec(4) 
  logical,optional :: use_bg_bufInput
  !
  integer,pointer :: IndexToGlobal(:,:)
  real(realk),pointer :: AlphaBeta(:,:)
  real(realk),pointer :: TMPAlphaBetaDecomp(:,:),AlphaCD3(:),NBAR(:,:)
  real(realk),pointer :: AlphaCD5(:),Calpha2(:),CalphaNAF(:)
  real(realk),pointer :: W(:,:),Wprime(:,:),NBARTMP(:,:)
  real(realk) :: TS3,TE3,MemInGBCollected,MemInGBCollected2,TS,TE,SumSV,FullSumSV
  real(realk) :: MemForFullAOINT,MemForFullMOINT,maxsize,MemForPartialMOINT
  real(realk) :: TS4,TE4
  TYPE(MoleculeInfo),pointer      :: molecule1,molecule2,molecule3,molecule4
  integer(kind=long)    :: nSize,n8
  integer(kind=ls_mpik) :: node 
  integer :: MaxNaux,M,N,K,ndimMax1,nbasisAuxMPI2(numnodes),MynbasisAuxMPI2
  integer :: nthreads,PerformReduction,Oper,nAuxMPI(numnodes),MaxnAuxMPI,I
  integer :: J,offset,offset2,inode,nbuf1,nbuf2,nbuf3,dim1,MinAuxBatch,ndimMax2
  integer :: GindexToLocal(nbasisAux),nbasisAuxMPI3(numnodes),NREDLOC,NRED
  real(realk) :: epsilon
  logical :: use_bg_buf
#ifdef VAR_OMP
  integer, external :: OMP_GET_MAX_THREADS
#endif
  CALL LSTIMER('START ',TS,TE,LUPRI,ForcePrint)
  CALL LSTIMER('START ',TS3,TE3,LUPRI,ForcePrint)
  use_bg_buf = .FALSE.
  IF(present(use_bg_bufInput)) use_bg_buf = use_bg_bufInput
  epsilon = DECinfo%NAFthreshold 
  IF(use_bg_buf.AND.DECinfo%NAF)call lsquit('bg_buf and NAF combi not tested',-1)
  IF(DECinfo%PL.GT.2)THEN
   IF(DECinfo%NAF)THEN
    WRITE(DECinfo%output,*)'Use Natural Auxiliary Functions (NAF) with Threshold = ',&
         & epsilon
   ENDIF
  ENDIF
  PerformReduction = -1
  NBA = 0
  MynbasisAuxMPI2 = 0 
  MaxNaux = 0
  print*,'use_bg_buf',use_bg_buf
  IF(use_bg_buf)THEN
     maxsize = mem_get_bg_buf_free()*8.0E-9_realk
  ELSE
     call get_currently_available_memory(MemInGBCollected)
     maxsize = MemInGBCollected*0.65E0_realk
  ENDIF
  call GetOperatorFromCharacter(Oper,intspec(4),mylsitem%Setting)

  IF(DECinfo%AuxAtomicExtent)THEN
     molecule1 => mylsitem%SETTING%MOLECULE(1)%p
     mylsitem%SETTING%MOLECULE(1)%p => mylsitem%INPUT%AUXMOLECULE
  ENDIF
  call determine_maxBatchOrbitalsize(DECinfo%output,&
       & Mylsitem%setting,MinAuxBatch,intspec(1))
  IF(DECinfo%AuxAtomicExtent)THEN
     mylsitem%SETTING%MOLECULE(1)%p => molecule1
  ENDIF

#ifdef VAR_OMP     
  !$OMP PARALLEL SHARED(nthreads)  
  !$OMP MASTER
  nthreads = OMP_GET_MAX_THREADS()
  !$OMP END MASTER
  !$OMP END PARALLEL
#else
  nthreads = 1
#endif
  
  IF(master)THEN
     !Memory requirement to have the full MO integral in memory 
     MemForFullMOINT = (nbasisAux*nvirt*nocc+MinAuxBatch*nthreads*(nbasis1*nocc+nbasis1*nbasis2)+&
          & nbasis1*nvirt+nbasis2*nocc)*8.0E-9_realk
!     print*,'MemForFullMOINT=',MemForFullMOINT,'maxsize',maxsize
!     print*,'MemForFullMOINT.LT.maxsize: ',MemForFullMOINT.LT.maxsize

!     WRITE(DECinfo%output,*)'DECinfo%RIMP2ForcePDMCalpha',DECinfo%RIMP2ForcePDMCalpha
     IF(DECinfo%RIMP2ForcePDMCalpha.OR.MemForFullMOINT.GE.maxsize)THEN
        !Full MO cannot fit in memory       
        PerformReduction = 0
        MaxNaux = 0
!        WRITE(DECinfo%output,*)'DECinfo%RIMP2ForcePDMCalpha: PerformReduction = ',PerformReduction
        IF(mylsitem%setting%scheme%ForceRIMP2memReduced)THEN 
           MaxNaux = MinAuxBatch + 1
        ENDIF
     ELSE
        !Full MO can fit on all nodes Which means we can do a reduction.        
        IF(DECinfo%PL.GT.0)THEN
           WRITE(DECinfo%output,'(A,F8.1,A)')'RIMP2: Full MO (alpha|cd) integral requires ',MemForFullMOINT,' GB'
           IF(.NOT.use_bg_buf)WRITE(DECinfo%output,'(A,F8.1,A)') &
                & 'RIMP2: Memory available                     ',MemInGBCollected,' GB'
           WRITE(DECinfo%output,'(A,F8.1,A)')'RIMP2: Memory available (65%)               ',maxsize,' GB'
           print*,'RIMP2: Full MO (alpha|cd) integral requires ',MemForFullMOINT,' GB'
           IF(.NOT.use_bg_buf)print*,'RIMP2: Memory available                     ',MemInGBCollected,' GB'
           print*,'RIMP2: Memory available (65%)               ',maxsize,' GB'
        ENDIF
        !Memory requirement to have the full AO integral in memory 
        MemForFullAOINT = MAX(nbasisAux*nvirt*nocc+nbasisAux*nbasis1*nocc,&
             & nbasisAux*nbasis1*nocc+nbasisAux*nbasis1*nbasis2)*8.0E-9_realk
        IF(DECinfo%PL.GT.0)THEN
           WRITE(DECinfo%output,'(A,F8.1,A)')'RIMP2: Full AO (alpha|cd) integral requires ',MemForFullAOINT,' GB'           
        ENDIF
        IF(MemForFullAOINT.LT.maxsize)THEN 
           !Full AO can fit on all nodes Which means we can do a reduction.
           MaxNaux = nbasisAux
           PerformReduction = MaxNaux
!           print*,'MemForFullMOINT.LT.maxsize: MaxNaux',MaxNaux
        ELSE !Full MO can fit on all nodes Which means we can do a reduction.
           !cannot have the full AO integral in memory 
           !so we build a subset (MaxNaux,nbasis1,nbasis2) and then 
           !transform to MO basis during the integral evaluation. 
           !The bigger MaxNaux the fewer times we have to do the AO to MO
           !Memory requirement = nbasisAux*nvirt*nocc+MaxNaux*nthreads*
           !  (nbasis1*nocc+nbasis1*nbasis2)+nbasis1*nvirt+nbasis2*nocc 
           !so we choose MaxNaux to be
           MaxNaux = MIN(nbasisAux,FLOOR((MaxSize/8.0E-9_realk-nbasisAux*nvirt*nocc-2*nbasis1*nvirt-2*nbasis2*nocc) &
                & /(nbasis1*nocc+nbasis1*nbasis2*nthreads)))
           PerformReduction = MaxNaux
!           print*,'MemForFullMOINT.LT.maxsize: MaxNaux',MaxNaux

           IF(DECinfo%PL.GT.0)THEN
              WRITE(DECinfo%output,'(A,F8.1,A)')'RIMP2: MaxNaux Determination: Memory available (65%)',maxsize,' GB'           
              WRITE(DECinfo%output,'(A,F8.1,A)')'RIMP2: MaxNaux Determination: Memory available after MO int + CMO',&
                   & (MaxSize-nbasisAux*nvirt*nocc*8.0E-9_realk-2*nbasis1*nvirt*8.0E-9_realk-2*nbasis2*nocc*8.0E-9_realk),' GB'           
              WRITE(DECinfo%output,'(A,I8,A,I8)')'RIMP2: MaxNaux Determination: MaxNaux',MaxNaux,' compared to nAux=',nbasisAux
              WRITE(DECinfo%output,'(A,F8.1,A)')'RIMP2: This MaxNaux Correspond to Memory usage of ',&
                   &   nbasisAux*nvirt*nocc*8.0E-9_realk+2*nbasis1*nvirt*8.0E-9_realk+2*nbasis2*nocc*8.0E-9_realk&
                   & + MaxNaux*nthreads*nbasis1*nbasis2*8.0E-9_realk + MaxNaux*nbasis1*nocc*8.0E-9_realk,' GB'
           ENDIF
        ENDIF
        IF(mylsitem%setting%scheme%ForceRIMP2memReduced)THEN 
           MaxNaux = MinAuxBatch + 1
           PerformReduction = MaxNaux
        ENDIF
     ENDIF
  ENDIF
#ifdef VAR_MPI
  call time_start_phase( PHASE_IDLE )
  call lsmpi_barrier(infpar%lg_comm)
  call time_start_phase( PHASE_COMM )
  call ls_mpibcast(PerformReduction,infpar%master,infpar%lg_comm)
  IF(PerformReduction.NE.0) MaxNaux = PerformReduction
  call time_start_phase( PHASE_WORK )
#endif
  CALL LSTIMER('DF_Calpha:Init ',TS3,TE3,LUPRI,ForcePrint)

  !=====================================================================================
  ! Master Obtains (alpha|beta) ERI in Auxiliary Basis 
  !=====================================================================================

  IF(AlphaBetaDecompCreate)THEN
     IF(master)THEN
!        call mem_alloc(AlphaBeta,nbasisAux,nbasisAux)
        IF(DECinfo%AuxAtomicExtent)THEN
           molecule1 => mylsitem%SETTING%MOLECULE(1)%p
           molecule2 => mylsitem%SETTING%MOLECULE(2)%p
           molecule3 => mylsitem%SETTING%MOLECULE(3)%p
           molecule4 => mylsitem%SETTING%MOLECULE(4)%p
           mylsitem%SETTING%MOLECULE(1)%p => mylsitem%INPUT%AUXMOLECULE
           mylsitem%SETTING%MOLECULE(2)%p => mylsitem%INPUT%AUXMOLECULE
           mylsitem%SETTING%MOLECULE(3)%p => mylsitem%INPUT%AUXMOLECULE
           mylsitem%SETTING%MOLECULE(4)%p => mylsitem%INPUT%AUXMOLECULE
        ENDIF
        CALL LSTIMER('START ',TS4,TE4,LUPRI,ForcePrint)
        call II_get_RI_AlphaBeta_2centerInt(DECinfo%output,DECinfo%output,&
             & AlphaBetaDecomp,mylsitem%setting,nbasisAux,Oper)
        CALL LSTIMER('DF_Calpha:AlphaBeta',TS4,TE4,LUPRI,ForcePrint)
        IF(DECinfo%AuxAtomicExtent)THEN
           mylsitem%SETTING%MOLECULE(1)%p => molecule1
           mylsitem%SETTING%MOLECULE(2)%p => molecule2
           mylsitem%SETTING%MOLECULE(3)%p => molecule3
           mylsitem%SETTING%MOLECULE(4)%p => molecule4
        ENDIF
        ! Create the inverse square root AlphaBeta = (alpha|beta)^(-1/2)
        ! Warning the inverse is not unique so in order to make sure all slaves have the same
        ! inverse matrix we calculate it on the master a BCAST to slaves
!        call lowdin_diag_S_minus_sqrt(nbasisAux, AlphaBeta,AlphaBetaDecomp,lupri)
        call Get_InverseCholeskyFactor(nbasisAux,AlphaBetaDecomp,lupri)
        CALL LSTIMER('DF_Calpha:AlphaBetaDecomp',TS4,TE4,LUPRI,ForcePrint)
!        call mem_dealloc(AlphaBeta)
     ENDIF
#ifdef VAR_MPI
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     call ls_mpibcast(AlphaBetaDecomp,nbasisAux,nbasisAux,infpar%master,infpar%lg_comm)
     call time_start_phase(PHASE_WORK)   
#endif
  ENDIF
  CALL LSTIMER('DF_Calpha:AlphaBeta',TS3,TE3,LUPRI,ForcePrint)

  !==================================================================
  ! Determine:  
  ! nbasisAuxMPI2 is the dimension of the Aux basis for each node
  ! used in Calpha(nbasisAuxMPI2,nvirt,nocc) 
  ! NOT to be confused with nAuxMPI which is the Aux basis for each node
  ! used in the 3 center integral (alpha|nvirt,nocc)
  ! also make TMPAlphaBetaDecomp(MynbasisAuxMPI2,nbasisAux) ready for DGEMM
  !==================================================================
  
  dim1 = nBasisAux

  !Output the Aux dimensions for each node nAuxMPI of the Integral
  IF(DECinfo%AuxAtomicExtent)THEN
     molecule1 => mylsitem%SETTING%MOLECULE(1)%p
     mylsitem%SETTING%MOLECULE(1)%p => mylsitem%INPUT%AUXMOLECULE
  ENDIF
  call II_get_RI_AlphaCD_3CenterIntFullOnAllNNdim(mylsitem%setting,&
       & nAuxMPI,IndexToGlobal,numnodes,MaxnAuxMPI,intspec(1),&
       & GindexToLocal,nbasisAux)
  IF(DECinfo%AuxAtomicExtent)THEN
     mylsitem%SETTING%MOLECULE(1)%p => molecule1
  ENDIF

  IF(CollaborateWithSlaves.OR.DECinfo%RIMP2ForcePDMCalpha)then 
     ndimMax1 = nbasisAux/numnodes
     do I=1,numnodes
        nbasisAuxMPI2(I) = ndimMax1
     enddo
     J=2 !not add to master
     do I=1,MOD(nbasisAux,numnodes)
        nbasisAuxMPI2(J) = nbasisAuxMPI2(J) + 1
        J=J+1
     enddo
     MynbasisAuxMPI2 = nbasisAuxMPI2(mynum+1)
     call mem_alloc(TMPAlphaBetaDecomp,MynbasisAuxMPI2,nbasisAux)
     call buildTMPAlphaBetaDecomp(TMPAlphaBetaDecomp,AlphaBetaDecomp,&
          & MynbasisAuxMPI2,nbasisAux,mynum,ndimMax1,numnodes)
     NBA = MynbasisAuxMPI2
     IF(PerformReduction.EQ.0)THEN
        MemForPartialMOINT = (nAuxMPI(mynum+1)*nvirt*nocc+&
             & MinAuxBatch*nthreads*(nbasis1*nocc+nbasis1*nbasis2)+&
             & nbasis1*nvirt+nbasis2*nocc)*8.0E-9_realk
        IF(MemForPartialMOINT.GE.maxsize)THEN !Error 
           CALL lsquit('Not enough memory in build_calpha bcast schem',-1)
        ELSE
           MaxNaux = MIN(nAuxMPI(mynum+1),&
                &FLOOR((MaxSize/8.0E-9_realk-nAuxMPI(mynum+1)*nvirt*nocc-&
                & nbasis1*nvirt-nbasis2*nocc) &
                & /((nbasis1*nocc+nbasis1*nbasis2)*nthreads)))
        ENDIF
        IF(mylsitem%setting%scheme%ForceRIMP2memReduced)THEN 
           MaxNaux = MinAuxBatch + 1
        ENDIF
        dim1 = nAuxMPI(mynum+1)
     ENDIF
  ELSE
     NBA = nbasisAux
     nAuxMPI = nbasisAux
  ENDIF

  IF(use_bg_buf)THEN
     !allocate now because I need to allocate in order due to push pop mechanisme
     nsize = NBA*nvirt*nocc
     call mem_pseudo_alloc(Calpha,nsize)
  ENDIF

  !=====================================================================================
  ! Obtain 3 center RI integrals (alpha,a,i) 
  !=====================================================================================
  
  IF(DECinfo%AuxAtomicExtent)THEN
     molecule1 => mylsitem%SETTING%MOLECULE(1)%p
     molecule2 => mylsitem%SETTING%MOLECULE(2)%p
     mylsitem%SETTING%MOLECULE(1)%p => mylsitem%INPUT%AUXMOLECULE
     mylsitem%SETTING%MOLECULE(2)%p => mylsitem%INPUT%AUXMOLECULE
  ENDIF
  IF(DECinfo%PL.GT.0)THEN
     IF(dim1.EQ.nbasisAux)THEN
        print*,'NEW (alpha|AI) CODE: Perform Reduction'
     ELSE
        print*,'NEW (alpha|AI) CODE: Perform BCAST Loop using:',dim1
     ENDIF
     IF(MaxNaux.LT.dim1)THEN
        print*,'NEW (alpha|AI) CODE: Internal AotoMO, ',MaxNaux,' out of ',dim1
     ELSE
        print*,'NEW (alpha|AI) CODE: external AotoMO, MaxNaux=',MaxNaux
     ENDIF
     print*,'NEW (alpha|AI) CODE: Number of Threads =',nthreads
  ENDIF
  !Output: AlphaCD3(dim1,nvirt,nocc) 
  call II_get_RI_AlphaCD_3CenterIntFullOnAllNN(DECinfo%output,DECinfo%output,&
       & AlphaCD3,mylsitem%setting,nbasisAux,nbasis1,nbasis2,intspec,MaxNaux,&
       & nvirt,nocc,.TRUE.,Cvirt,Cocc,nthreads,dim1,GindexToLocal,DECinfo%PL,&
       & use_bg_buf)

  CALL LSTIMER('DF_Calpha:3CenterInt',TS3,TE3,LUPRI,ForcePrint)

!  call sleep(mynum*10)
!  print*,'XXX AlphaCD3 dim1=',dim1,'mynum',mynum
!  call ls_output(AlphaCD3,1,dim1,1,4,dim1,nvirt*nocc,1,6)
  
  IF(DECinfo%AuxAtomicExtent)THEN
     mylsitem%SETTING%MOLECULE(1)%p => molecule1
     mylsitem%SETTING%MOLECULE(2)%p => molecule2
  ENDIF
  
#ifdef VAR_MPI
  IF(PerformReduction.NE.0)THEN
     !all nodes have filled in part of alphaCD3(nbasisAux,nvirt,nocc)
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     n8 = nbasisAux*nvirt*nocc
     call lsmpi_allreduce(alphaCD3,n8,infpar%lg_comm)
     call time_start_phase( PHASE_WORK )
  ENDIF
#endif
  IF(PerformReduction.NE.0)THEN
     IF(DECinfo%NAF)THEN
        IF(CollaborateWithSlaves)then 
           call mem_dealloc(TMPAlphaBetaDecomp)
        ENDIF
        call mem_alloc(Wprime,nbasisAux,nbasisAux)
        call RIMP2_buildWprimeFromAlphaCD(AlphaCD3,nbasisAux,nocc,nvirt,Wprime,mynum,numnodes)
#ifdef VAR_MPI
        !Reduction of Wprime
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        call lsmpi_reduction(Wprime,nbasisAux,nbasisAux,infpar%master,infpar%lg_comm)
        call time_start_phase( PHASE_WORK )
#endif
        IF(master)THEN
           call mem_alloc(W,nbasisAux,nbasisAux)
           call NAF_buildW(W,Wprime,AlphaBetaDecomp,nbasisAux)
           call NAF_SVD_W(W,Wprime,NBAR,nbasisAux,epsilon,NRED,SumSV,FullSumSV) 
           IF(DECinfo%PL.GT.2)THEN
              WRITE(DECinfo%output,*)'NAF: Auxiliary functions         = ',nbasisAux
              WRITE(DECinfo%output,*)'NAF: Natural Auxiliary functions = ',NRED
              WRITE(DECinfo%output,*)'NAF: Sum of included eigenvalues = ',SumSV
              WRITE(DECinfo%output,*)'NAF: Sum of neglected eigenvalues= ',FullSumSV-SumSV
              WRITE(DECinfo%output,*)'NAF: Sum of all eigenvalues      = ',FullSumSV
           ENDIF
           call mem_dealloc(W)
        ENDIF
        call mem_dealloc(Wprime)
#ifdef VAR_MPI
        !BCAST OF NBAR
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        call ls_mpibcast(NRED,infpar%master,infpar%lg_comm)
        nbuf1 = NRED
        nbuf2 = nbasisAux
        IF(.NOT.master)call mem_alloc(NBAR,NRED,nbasisAux)
        call ls_mpibcast(NBAR,nbuf1,nbuf2,infpar%master,infpar%lg_comm)
        call time_start_phase( PHASE_WORK )
#endif
        ndimMax2 = NRED/numnodes
        do I=1,numnodes
           nbasisAuxMPI3(I) = ndimMax2
        enddo
        J=2 !not add to master
        do I=1,MOD(NRED,numnodes)
           nbasisAuxMPI3(J) = nbasisAuxMPI3(J) + 1
           J=J+1
        enddo
        NREDLOC = nbasisAuxMPI3(mynum+1)
        call mem_alloc(NBARTMP,NREDLOC,nbasisAux)
        call buildNBARTMP(NBARTMP,NBAR,NREDLOC,NRED,nbasisAux,mynum,ndimMax2,numnodes)
        call mem_dealloc(NBAR)
        call mem_alloc(W,NREDLOC,nbasisAux) !used as TMP
        M =  NREDLOC         !rows of Output Matrix
        N =  nbasisAux       !columns of Output Matrix
        K =  nbasisAux       !summation dimension
        call dgemm('N','N',M,N,K,1.0E0_realk,NBARTMP,M,AlphaBetaDecomp,K,0.0E0_realk,W,M)
        call mem_dealloc(NBARTMP)
        M = NREDLOC          !rows of Output Matrix
        N = nvirt*nocc       !columns of Output Matrix
        K = nbasisAux        !summation dimension
        nsize = NREDLOC*nvirt*nocc
        IF(.NOT.use_bg_buf) call mem_alloc(Calpha,nsize)
        call dgemm('N','N',M,N,K,1.0E0_realk,W,M,AlphaCD3,K,0.0E0_realk,Calpha,M)
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(alphaCD3)
           !This looks weird but Calpha is currently pointing to the first 1:NBA*nvirt*nocc elements
           !of a "permanent" memory array. I only need the first 1:NREDLOC*nvirt*nocc elements
           !so I shrink the array dimension by deassociating (NOT deallocating) and reassociate
           call mem_pseudo_dealloc(Calpha)
           nsize = NREDLOC*nvirt*nocc
           call mem_pseudo_alloc(Calpha,nsize)
        ELSE
           call mem_dealloc(alphaCD3)
        ENDIF
        call mem_dealloc(W)
        NBA = NREDLOC
     ELSE
        IF(CollaborateWithSlaves)then 
           nsize = NBA*nvirt*(nocc*i8)
           IF(.NOT.use_bg_buf)call mem_alloc(Calpha,nsize)
           !Calpha = TMPAlphaBetaDecomp(MynbasisAuxMPI,nbasisAux)
           M =  NBA              !rows of Output Matrix
           N =  nvirt*nocc       !columns of Output Matrix
           K =  nbasisAux        !summation dimension
           call dgemm('N','N',M,N,K,1.0E0_realk,TMPAlphaBetaDecomp,M,alphaCD3,K,0.0E0_realk,Calpha,M)
           call mem_dealloc(TMPAlphaBetaDecomp)
        ELSE !Serial version        
           nsize = nbasisAux*nvirt*nocc
           IF(.NOT.use_bg_buf)call mem_alloc(Calpha,nsize)
           M =  nbasisAux        !rows of Output Matrix
           N =  nvirt*nocc       !columns of Output Matrix
           K =  nbasisAux        !summation dimension
           call dgemm('N','N',M,N,K,1.0E0_realk,AlphaBetaDecomp,M,AlphaCD3,K,0.0E0_realk,Calpha,M)
        ENDIF
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(alphaCD3)
        ELSE
           call mem_dealloc(alphaCD3)
        ENDIF
     ENDIF
     !  print*,'MY RIMP2 INTEGRAL AlphaCD2(1:nA,1:4) NEW VERSION MYNUM',MYNUM
     !  call ls_output(AlphaCD3,1,size(AlphaCD3,1),1,4,size(AlphaCD3,1),nvirt*nocc,1,6)
     CALL LSTIMER('DF_Calpha:Calpha',TS3,TE3,LUPRI,ForcePrint)
  ELSE
     !=====================================================================================
     ! MPI scheme:  Bcast Routine
     !=====================================================================================
!     print*,'NBA',NBA     
     IF(.NOT.use_bg_buf)THEN
        nsize = NBA*nvirt*(nocc*i8)
        call mem_alloc(Calpha,nsize)
     ENDIF
     call ls_dzero8(Calpha,nsize)
     IF(DECinfo%NAF)THEN
        call mem_alloc(Wprime,nbasisAux,nbasisAux)
        nsize = nbasisAux*nbasisAux
        call ls_dzero8(Wprime,nsize)
     ENDIF
     DO inode = 1,numnodes
        nbuf1 = nAuxMPI(inode)
        nbuf2 = nvirt
        nbuf3 = nocc
        IF(mynum.EQ.inode-1)THEN
#ifdef VAR_MPI
           call time_start_phase( PHASE_IDLE )
           call lsmpi_barrier(infpar%lg_comm)
           call time_start_phase( PHASE_COMM )
           node = mynum
           nsize = nbuf1*nbuf2*nbuf3
           call ls_mpibcast(AlphaCD3,nsize,node,infpar%lg_comm)
#endif
           call RIMP2_buildCalphaFromAlphaCD(nocc,nvirt,numnodes,&
                & nAuxMPI,IndexToGlobal,MaxnAuxMPI,AlphaCD3,nAuxMPI(inode),&
                & Calpha,NBA,TMPAlphaBetaDecomp,nbasisAux,inode)
           IF(DECinfo%NAF)THEN
              call RIMP2_buildWprimeFromAlphaCD1(nocc,nvirt,numnodes,&
                   & nAuxMPI,IndexToGlobal,MaxnAuxMPI,AlphaCD3,nAuxMPI(inode),&
                   & Wprime,nbasisAux,inode)              
           ENDIF
        ELSE
           node = inode-1
           !recieve
           nsize = nAuxMPI(inode)*nvirt*(nocc*i8)
           IF(use_bg_buf)THEN
              call mem_pseudo_alloc(AlphaCD5,nsize)
           ELSE
              call mem_alloc(AlphaCD5,nsize)
           ENDIF
#ifdef VAR_MPI
           call time_start_phase( PHASE_IDLE )
           call lsmpi_barrier(infpar%lg_comm)
           call time_start_phase( PHASE_COMM )
           nsize = nbuf1*nbuf2*nbuf3
           call ls_mpibcast(AlphaCD5,nsize,node,infpar%lg_comm)
           call time_start_phase( PHASE_WORK )
#endif
           call RIMP2_buildCalphaFromAlphaCD(nocc,nvirt,numnodes,&
                & nAuxMPI,IndexToGlobal,MaxnAuxMPI,AlphaCD5,nAuxMPI(inode),&
                & Calpha,NBA,TMPAlphaBetaDecomp,nbasisAux,inode)
           IF(DECinfo%NAF)THEN
              call RIMP2_buildWprimeFromAlphaCD2(nocc,nvirt,numnodes,&
                   & nAuxMPI,IndexToGlobal,MaxnAuxMPI,AlphaCD5,nAuxMPI(inode),&
                   & Wprime,nbasisAux,inode,AlphaCD3,nAuxMPI(mynum+1),mynum)
           ENDIF
           IF(use_bg_buf)THEN
              call mem_pseudo_dealloc(AlphaCD5)
           ELSE
              call mem_dealloc(AlphaCD5)
           ENDIF
        ENDIF
     ENDDO
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(AlphaCD3)
     ELSE
        call mem_dealloc(AlphaCD3)
     ENDIF
     CALL LSTIMER('DF_Calpha:Calpha',TS3,TE3,LUPRI,ForcePrint)
     IF(DECinfo%NAF)THEN
#ifdef VAR_MPI
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        call lsmpi_reduction(Wprime,nbasisAux,nbasisAux,infpar%master,infpar%lg_comm)
        call time_start_phase( PHASE_WORK )
#endif
        IF(master)THEN
           call mem_alloc(W,nbasisAux,nbasisAux)
           call NAF_buildW(W,Wprime,AlphaBetaDecomp,nbasisAux)           
           call NAF_SVD_W(W,Wprime,NBAR,nbasisAux,epsilon,NRED,SumSV,FullSumSV) 
           IF(DECinfo%PL.GT.2)THEN
              WRITE(DECinfo%output,*)'NAF: Auxiliary functions         = ',nbasisAux
              WRITE(DECinfo%output,*)'NAF: Natural Auxiliary functions = ',NRED
              WRITE(DECinfo%output,*)'NAF: Sum of included eigenvalues = ',SumSV
              WRITE(DECinfo%output,*)'NAF: Sum of neglected eigenvalues= ',FullSumSV-SumSV
              WRITE(DECinfo%output,*)'NAF: Sum of all eigenvalues      = ',FullSumSV
           ENDIF
           call mem_dealloc(W)
        ENDIF
        call mem_dealloc(Wprime)
#ifdef VAR_MPI
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        call ls_mpibcast(NRED,infpar%master,infpar%lg_comm)
        nbuf1 = NRED
        nbuf2 = nbasisAux
        IF(.NOT.master)call mem_alloc(NBAR,NRED,nbasisAux)
        call ls_mpibcast(NBAR,nbuf1,nbuf2,infpar%master,infpar%lg_comm)
        call time_start_phase( PHASE_WORK )
#endif
        ndimMax2 = NRED/numnodes
        do I=1,numnodes
           nbasisAuxMPI3(I) = ndimMax2
        enddo
        J=2 !not add to master
        do I=1,MOD(NRED,numnodes)
           nbasisAuxMPI3(J) = nbasisAuxMPI3(J) + 1
           J=J+1
        enddo
        NREDLOC = nbasisAuxMPI3(mynum+1)
        nsize = NREDLOC*nvirt*nocc
        IF(.NOT.use_bg_buf)call mem_alloc(CalphaNAF,nsize)
        call ls_dzero8(CalphaNAF,nsize)
        call mem_alloc(NBARTMP,NREDLOC,nbasisAux)
        call buildNBARTMP(NBARTMP,NBAR,NREDLOC,NRED,nbasisAux,mynum,ndimMax2,numnodes)
        call mem_dealloc(NBAR)
        DO inode = 1,numnodes
           nbuf1 = nbasisAuxMPI2(inode) !dim1 of Calpha()
           nbuf2 = nvirt
           nbuf3 = nocc
           IF(mynum.EQ.inode-1)THEN
#ifdef VAR_MPI
              call time_start_phase( PHASE_IDLE )
              call lsmpi_barrier(infpar%lg_comm)
              call time_start_phase( PHASE_COMM )
              node = mynum
              nsize = nbuf1*nbuf2*nbuf3
              call ls_mpibcast(Calpha,nsize,node,infpar%lg_comm)
#endif
              offset = (inode-1)*ndimMax1
              offset2 = numnodes*ndimMax1 + inode - 1 
              call NAF_buildCalphaNAF(Calpha,nbasisAuxMPI2(inode),nvirt,nocc,&
                   & CalphaNAF,NREDLOC,NBARTMP,nbasisAux,offset,&
                   & offset2,ndimMax1,numnodes)
           ELSE
              node = inode-1
              !recieve
              nsize = nbasisAuxMPI2(inode)*nvirt*nocc
              IF(use_bg_buf)THEN
                 call mem_pseudo_alloc(Calpha2,nsize)
              ELSE
                 call mem_alloc(Calpha2,nsize)
              ENDIF
#ifdef VAR_MPI
              call time_start_phase( PHASE_IDLE )
              call lsmpi_barrier(infpar%lg_comm)
              call time_start_phase( PHASE_COMM )
              nsize = nbuf1*nbuf2*nbuf3
              call ls_mpibcast(Calpha2,nsize,node,infpar%lg_comm)
              call time_start_phase( PHASE_WORK )
#endif
              
              offset = (inode-1)*ndimMax1
              offset2 = numnodes*ndimMax1 + inode - 1 
              call NAF_buildCalphaNAF(Calpha2,nbasisAuxMPI2(inode),nvirt,nocc,&
                   & CalphaNAF,NREDLOC,NBARTMP,nbasisAux,offset,&
                   & offset2,ndimMax1,numnodes)
              IF(use_bg_buf)THEN
                 call mem_pseudo_dealloc(Calpha2)
              ELSE
                 call mem_dealloc(Calpha2)
              ENDIF
           ENDIF
        ENDDO
        call mem_dealloc(NBARTMP)
        NBA = NREDLOC
        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(Calpha)
           Calpha => CalphaNAF
           call lsquit('clearly not working',-1)
        ELSE
           call mem_dealloc(Calpha)
           Calpha => CalphaNAF           
        ENDIF
        CALL LSTIMER('DF_Calpha:NAF',TS3,TE3,LUPRI,ForcePrint)
     ENDIF
     !allocated inside II_get_RI_AlphaCD_3CenterIntFullOnAllNNdim
     IF(CollaborateWithSlaves.OR.DECinfo%RIMP2ForcePDMCalpha)THEN
        call mem_dealloc(TMPAlphaBetaDecomp)        
     ENDIF
  ENDIF
  call mem_dealloc(IndexToGlobal) 

!  call sleep(mynum*5)
!  PRINT*,'Build_CalphaMO2:Nba',Nba
!  WRITE(6,*)'Build_CalphaMO2:Final Calph(NBA=',NBA,',nvirt=',nvirt,',nocc=',nocc,')'
!  WRITE(6,*)'Build_CalphaMO2:Print Subset Final Calph(NBA=',NBA,',1:4)  MYNUM',MYNUM
!  call ls_output(Calpha,1,NBA,1,4,NBA,nvirt*nocc,1,6)
!
  CALL LSTIMER('Build_CalphaMO2',TS,TE,LUPRI,ForcePrint)

end subroutine Build_CalphaMO2

!computes the Cholesky factorization of the inverse of 
!a real symmetric positive definite matrix A = U^T * U 
subroutine Get_InverseCholeskyFactor(n,A,lupri)
  implicit none
  integer, intent(in)          :: n,lupri
  real(realk), intent(inout)   :: A(n,n)
  !
  integer                      :: i,j,np,k,info
  real(realk) :: TS,TE
  call LSTIMER('START ',TS,TE,lupri,.TRUE.)
  call DPOTRF('U', N, A, N, INFO ) !U=Cholesky factor
  IF(INFO.ne. 0) THEN
     print *, 'DPOTRF NR 1 Failed in Get_InverseCholeskyFactor',INFO
     call lsquit('DPOTRF NR 1 Failed in Get_InverseCholeskyFactor',-1)
  ENDIF
  call DPOTRI('U', N, A, N, INFO ) !U=inverse of a original U
  IF(INFO.ne. 0) THEN
     print *, 'DPOTRI Failed in Get_InverseCholeskyFactor',INFO
     call lsquit('DPOTRI Failed in Get_InverseCholeskyFactor',-1)
  ENDIF
  call DPOTRF('U', N, A, N, INFO ) !U=Cholesky factor of inverse of a original U
  IF(INFO.ne. 0) THEN
     print *, 'DPOTRF NR 2 Failed in Get_InverseCholeskyFactor',INFO
     call lsquit('DPOTRF NR 2 Failed in Get_InverseCholeskyFactor',-1)
  ENDIF
  do J=1,N
     do I=J+1,N
        A(I,J) = 0.0E0_realk
     enddo
  enddo
  call LSTIMER('Get_CholeskyFactor',TS,TE,lupri,.TRUE.)
end subroutine Get_InverseCholeskyFactor

subroutine NAF_buildCalphaNAF(Calpha,nAux2,nvirt,nocc,&
     & CalphaNAF,NREDLOC,NBARTMP,nbasisAux,&
     & offset,offset2,ndimMax,numnodes)
  implicit none 
  integer,intent(in) :: nAux2,nvirt,nocc,offset,offset2
  integer,intent(in) :: ndimMax,numnodes,nbasisAux,NREDLOC
  real(realk),intent(in) :: Calpha(nAux2,nvirt*nocc)
  real(realk),intent(in) :: NBARTMP(NREDLOC,nbasisAux)
  real(realk),intent(inout) :: CalphaNAF(NREDLOC,nvirt*nocc)
  !
  integer :: J,IREDLOC,AUX
  real(realk) :: TMP
  IF(nAux2.GT.ndimMax)THEN
     !$OMP PARALLEL DO DEFAULT(none) PRIVATE(IREDLOC,J,AUX,&
     !$OMP TMP) SHARED(ndimMax,nvirt,nocc,NREDLOC,NBARTMP,&
     !$OMP Calpha,offset,offset2,nAux2,CalphaNAF)
     DO J=1,nvirt*nocc
        DO IREDLOC=1,NREDLOC
           TMP = 0.0E0_realk
           DO AUX=1,ndimMax
              TMP = TMP + NBARTMP(IREDLOC,offset+AUX)*Calpha(AUX,J)
           ENDDO
           TMP = TMP + NBARTMP(IREDLOC,offset2)*Calpha(nAux2,J)
           CalphaNAF(IREDLOC,J) = CalphaNAF(IREDLOC,J) + TMP
        ENDDO
     ENDDO
     !$OMP END PARALLEL DO
  ELSE
     !$OMP PARALLEL DO DEFAULT(none) PRIVATE(IREDLOC,J,AUX,&
     !$OMP TMP) SHARED(ndimMax,nvirt,nocc,NREDLOC,NBARTMP,&
     !$OMP Calpha,offset,CalphaNAF)
     DO J=1,nvirt*nocc
        DO IREDLOC=1,NREDLOC
           TMP = 0.0E0_realk
           DO AUX=1,ndimMax
              TMP = TMP + NBARTMP(IREDLOC,offset+AUX)*Calpha(AUX,J)
           ENDDO
           CalphaNAF(IREDLOC,J) = CalphaNAF(IREDLOC,J) + TMP
        ENDDO
     ENDDO
     !$OMP END PARALLEL DO
  ENDIF
end subroutine NAF_buildCalphaNAF

subroutine buildTMPAlphaBetaDecomp(TMPAlphaBetaDecomp,AlphaBetaDecomp,&
     & MynbasisAuxMPI2,nbasisAux,mynum,ndimMax1,numnodes)
implicit none
integer,intent(in) :: MynbasisAuxMPI2,nbasisAux,mynum,ndimMax1,numnodes
real(realk),intent(in) :: AlphaBetaDecomp(nbasisAux,nbasisAux)
real(realk),intent(inout) :: TMPAlphaBetaDecomp(MynbasisAuxMPI2,nbasisAux)
!local variables
integer :: offset,offset2,I,J
offset = mynum*ndimMax1
offset2 = numnodes*ndimMax1 + mynum -1 +1
IF(MynbasisAuxMPI2.GT.ndimMax1)THEN
   !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(nbasisAux,ndimMax1,&
   !$OMP TMPAlphaBetaDecomp,AlphaBetaDecomp,offset,offset2)
   do I=1,nbasisAux
      do J=1,ndimMax1
         TMPAlphaBetaDecomp(J,I) = AlphaBetaDecomp(offset+J,I)
      enddo
      TMPAlphaBetaDecomp(ndimMax1+1,I) = AlphaBetaDecomp(offset2,I)
   enddo
   !$OMP END PARALLEL DO
ELSE
   !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(nbasisAux,ndimMax1,&
   !$OMP TMPAlphaBetaDecomp,AlphaBetaDecomp,offset)
   do I=1,nbasisAux
      do J=1,ndimMax1
         TMPAlphaBetaDecomp(J,I) = AlphaBetaDecomp(offset+J,I)
      enddo
   enddo
   !$OMP END PARALLEL DO
ENDIF
end subroutine buildTMPAlphaBetaDecomp

subroutine buildNBARTMP(NBARTMP,NBAR,NREDLOC,NRED,nbasisAux,mynum,ndimMax,numnodes)
  implicit none
  integer,intent(in) :: NREDLOC,nbasisAux,mynum,ndimMax,NRED,numnodes
  real(realk),intent(in) :: NBAR(NRED,nbasisAux)
  real(realk),intent(inout) :: NBARTMP(NREDLOC,nbasisAux)
  !local variables
  integer :: offset,offset2,I,J
  offset = mynum*ndimMax
  offset2 = numnodes*ndimMax + mynum -1 +1
  IF(NREDLOC.GT.ndimMax)THEN
     !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(nbasisAux,ndimMax,&
     !$OMP NBARTMP,NBAR,offset,offset2)
     do I=1,nbasisAux
        do J=1,ndimMax
           NBARTMP(J,I) = NBAR(offset+J,I)
        enddo
        NBARTMP(ndimMax+1,I) = NBAR(offset2,I)
     enddo
     !$OMP END PARALLEL DO
  ELSE
     !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(nbasisAux,ndimMax,&
     !$OMP NBARTMP,NBAR,offset)
     do I=1,nbasisAux
        do J=1,ndimMax
           NBARTMP(J,I) = NBAR(offset+J,I)
        enddo
     enddo
     !$OMP END PARALLEL DO
  ENDIF
end subroutine buildNBARTMP

subroutine NAF_SVD_W(W,TMP,NBAR,N,epsilon,nred,SumSV,FullSumSV)
  implicit none
  integer,intent(in)        :: N
  real(realk),intent(in)    :: W(N,N)
  real(realk),intent(inout) :: TMP(N,N),SumSV,FullSumSV
  real(realk),pointer       :: NBAR(:,:)
  real(realk),intent(in)    :: epsilon
  integer,intent(inout)     :: nred
  !local variables
  integer                :: lwork,INFO,I,K,J,infdiag,liwork
  real(realk), pointer   :: work(:),U(:,:),VT(:,:)
  integer,pointer        :: IPVT(:)
  real(realk)            :: RCOND, dummy(2),maxSV,SVm1,idummy(2),SV(N)
!
  real(realk),pointer :: W2(:,:)
  logical :: doSVD
  infdiag = 0
!  doSVD = .FALSE. ! do diagonalization instead (faster)

  !Perform a SVD  decomposition 
  ! W = U * SIGMA * transpose(V)
  ! where SIGMA is an M-by-N matrix which is zero except for 
  ! its min(m,n) diagonal elements
  ! for efficient storage the SIGMA non zero elements are stored in 
  ! SV (non singular values)
  !only the first min(m,n) columns of U (the left singular
  !vectors) are returned in the array U;
  !call mem_alloc(U,n,N)  Use TMP
  !only the first min(m,n) rows of V**T (the transposed right singular
  !vectors) are returned in the array V;
!  call mem_alloc(U,N,N)
  !S(n,m) = U(n,N) SV(N) VT(N,m)
!!$  IF(doSVD)THEN
!!$     call mem_alloc(VT,N,N)
!!$     lwork = -1      !workspace query
!!$     call dgesvd('S','S',N,N,W,N,SV,TMP,N,VT,N,dummy,lwork,INFO)
!!$     lwork = dummy(1)
!!$     call mem_alloc(work,lwork)
!!$     call dgesvd('S','S',N,N,W,N,SV,TMP,N,VT,N,work,lwork,INFO)
!!$     call mem_dealloc(VT)
!!$     !content of W destroyed
!!$     IF(INFO.NE.0)THEN
!!$        print*,'dgesvd in NAF_SVD_W failed  INFO=',INFO
!!$     ENDIF
!!$     call mem_dealloc(work)
!!$  ELSE
     call my_dsyev('V', 'U', N, W, SV)
!!$  ENDIF

  nred = 0 
  SumSV = 0.0E0_realk
  FullSumSV = 0.0E0_realk
  DO I=1,N
     IF(DECinfo%PL.GT.2) print*,'W SV(I)',SV(I),'SV(I).GT.epsilon',SV(I).GT.epsilon
     IF(SV(I).GT.epsilon)THEN
        nred = nred + 1
        SumSV = SumSV + SV(I)
     ENDIF
     FullSumSV = FullSumSV + SV(I)
  ENDDO
  call mem_alloc(NBAR,nred,N)
  nred = 0
!!$  IF(doSVD)THEN
!!$     DO I=1,N
!!$        IF(ABS(SV(I)).GT.epsilon)THEN
!!$           nred = nred + 1
!!$           DO J=1,N
!!$              !NBAR(nred,J) = U(n,nmin)
!!$              NBAR(nred,J) = TMP(J,I)
!!$           ENDDO
!!$        ENDIF
!!$     ENDDO
!!$  ELSE
     DO I=1,N
        IF(SV(I).GT.epsilon)THEN
           nred = nred + 1
           DO J=1,N
              NBAR(nred,J) = W(J,I)
           ENDDO
        ENDIF
     ENDDO
!!$  ENDIF
end subroutine NAF_SVD_W

subroutine NAF_buildW(W,Wprime,AlphaBetaDecomp,N) 
  implicit none
  integer,intent(in) :: N
  real(realk),intent(in)    :: Wprime(N,N)
  real(realk),intent(in)    :: AlphaBetaDecomp(N,N) !overlap^(-1/2)
  real(realk),intent(inout) :: W(N,N)
  !local variables
  real(realk),pointer :: TMP(:,:)
  !W(A,B) = Decomp(A,C)*Wprime(C,D)*Decomp(D,B)
  call mem_alloc(TMP,N,N)
  call dgemm('N','T',N,N,N,1.0E0_realk,Wprime,N,AlphaBetaDecomp,N,0.0E0_realk,TMP,N)
  call dgemm('N','N',N,N,N,1.0E0_realk,AlphaBetaDecomp,N,TMP,N,0.0E0_realk,W,N)
  call mem_dealloc(TMP)
end subroutine NAF_buildW

subroutine RIMP2_buildWprimeFromAlphaCD(AlphaCDl,nbasisAux,nocc,nvirt,&
     & Wprime,mynum,numnodes)
implicit none
integer,intent(in) :: nocc,nvirt,nbasisAux,mynum,numnodes
real(realk),intent(in) :: AlphaCDl(nbasisAux,nvirt*nocc)
real(realk),intent(inout) :: Wprime(nbasisAux,nbasisAux)
!
integer :: IB,B,A,StartB
real(realk) :: TMP
IB=1+mynum !1,2,3
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(B,A,TMP) &
!$OMP SHARED(nocc,nvirt,nbasisAux,AlphaCDl,Wprime,IB)
DO B = 1,nbasisAux
   TMP = AlphaCDl(B,IB)
   DO A = 1,nbasisAux
      Wprime(A,B) = AlphaCDl(A,IB)*TMP
   enddo
ENDDO
!$OMP END PARALLEL DO
StartB = 1+numnodes+mynum !4,5,6
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
!$OMP PRIVATE(IB,B,A,TMP) &
!$OMP SHARED(nocc,nvirt,nbasisAux,AlphaCDl,Wprime,StartB,numnodes)
do IB = StartB,nocc*nvirt,numnodes
   DO B = 1,nbasisAux
      TMP = AlphaCDl(B,IB)
      DO A = 1,nbasisAux
         Wprime(A,B) = Wprime(A,B) + AlphaCDl(A,IB)*TMP
      enddo
   ENDDO
ENDDO
!$OMP END PARALLEL DO
end subroutine RIMP2_buildWprimeFromAlphaCD

subroutine RIMP2_buildWprimeFromAlphaCD1(nocc,nvirt,numnodes,&
     & nAuxMPI,IndexToGlobal3,MaxnAuxMPI,AlphaCDk,nAux2,&
     & Wprime,nbasisAux,integralnum)
implicit none
integer,intent(in) :: nocc,nvirt,numnodes,nbasisAux,MaxnAuxMPI
integer,intent(in) :: nAuxMPI(numnodes),integralnum,nAux2
integer,intent(in) :: IndexToGlobal3(MaxnAuxMPI,numnodes)
real(realk),intent(in) :: AlphaCDk(nAux2,nvirt*nocc)
real(realk),intent(inout) :: Wprime(nbasisAux,nbasisAux)
!
integer :: IB,B,A,BETA,ALPHA
real(realk) :: TMP
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
!$OMP PRIVATE(B,A,BETA,ALPHA,TMP,IB) &
!$OMP SHARED(nocc,nvirt,nAux2,IndexToGlobal3,&
!$OMP        AlphaCDk,Wprime,integralnum)
DO B = 1,nAux2
   DO A = 1,nAux2
      BETA = IndexToGlobal3(B,integralnum)   
      ALPHA = IndexToGlobal3(A,integralnum)   
      TMP = 0.0E0_realk
      do IB = 1,nocc*nvirt
         TMP = TMP + AlphaCDk(A,IB)*AlphaCDk(B,IB)
      enddo
      Wprime(ALPHA,BETA) = TMP
   ENDDO
ENDDO
!$OMP END PARALLEL DO
end subroutine RIMP2_buildWprimeFromAlphaCD1

subroutine RIMP2_buildWprimeFromAlphaCD2(nocc,nvirt,numnodes,&
     & nAuxMPI,IndexToGlobal4,MaxnAuxMPI,AlphaCDk,nAux2,&
     & Wprime,nbasisAux,integralnum,AlphaCDj,nAux3,mynum)
implicit none
integer,intent(in) :: nocc,nvirt,numnodes,nbasisAux,MaxnAuxMPI
integer,intent(in) :: nAuxMPI(numnodes),integralnum,nAux2,mynum
integer,intent(in) :: IndexToGlobal4(MaxnAuxMPI,numnodes),nAux3
real(realk),intent(in) :: AlphaCDk(nAux2,nvirt*nocc)
real(realk),intent(in) :: AlphaCDj(nAux3,nvirt*nocc)
real(realk),intent(inout) :: Wprime(nbasisAux,nbasisAux)
!
integer :: IB,B,A,BETA,ALPHA
real(realk) :: TMP
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
!$OMP PRIVATE(B,A,BETA,ALPHA,TMP,IB) &
!$OMP SHARED(nocc,nvirt,nAux2,nAux3,IndexToGlobal4,&
!$OMP        AlphaCDj,AlphaCDk,Wprime,integralnum,mynum)
DO B = 1,nAux2
   DO A = 1,nAux3
      BETA = IndexToGlobal4(B,integralnum)   
      ALPHA = IndexToGlobal4(A,mynum+1)
      TMP = 0.0E0_realk
      do IB = 1,nocc*nvirt
         TMP = TMP + AlphaCDj(A,IB)*AlphaCDk(B,IB)
      enddo
      Wprime(ALPHA,BETA) = TMP
   ENDDO
ENDDO
!$OMP END PARALLEL DO
end subroutine RIMP2_buildWprimeFromAlphaCD2

subroutine RIMP2_buildCalphaFromAlphaCD(nocc,nvirt,numnodes,&
     & nAuxMPI,IndexToGlobal2,MaxnAuxMPI,AlphaCDi,nAux2,Calpha,NBA,TMPAlphaBetaDecomp,&
     & nbasisAux,integralnum)
  implicit none
  integer,intent(in) :: nocc,nvirt,numnodes,nbasisAux,MaxnAuxMPI
  integer,intent(in) :: nAuxMPI(numnodes),integralnum,nAux2,NBA
  integer,intent(in) :: IndexToGlobal2(MaxnAuxMPI,numnodes)
  real(realk),intent(in) :: AlphaCDi(nAux2,nvirt*nocc)
  real(realk),intent(in) :: TMPAlphaBetaDecomp(NBA,nbasisAux)
  real(realk),intent(inout) :: Calpha(NBA,nvirt*nocc)
  !
  integer :: IB,B,BETA,ALPHA
  real(realk) :: TMP
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
  !$OMP PRIVATE(IB,B,BETA,ALPHA,TMP) &
  !$OMP SHARED(nocc,nvirt,nAux2,IndexToGlobal2,NBA,&
  !$OMP        AlphaCDi,Calpha,TMPAlphaBetaDecomp,integralnum)
  do IB = 1,nocc*nvirt
     DO B = 1,nAux2
        TMP = AlphaCDi(B,IB)
        BETA = IndexToGlobal2(B,integralnum)
        do ALPHA = 1,NBA
           Calpha(ALPHA,IB) = Calpha(ALPHA,IB) + TMPAlphaBetaDecomp(ALPHA,BETA)*TMP
        enddo
     ENDDO
  enddo
!$OMP END PARALLEL DO
end subroutine RIMP2_buildCalphaFromAlphaCD

subroutine PlugInToalphaCDFull(mynum,nAtomsMPI,startAuxMPI,nocc,&
     & nvirt,nAuxMPI,alphaCDFull,alphaCD3,nbasisAux,MynbasisAuxMPI,&
     & numnodes,nAtomsAux)
  implicit none
  integer,intent(in) :: mynum,numnodes,nocc,nvirt,nAtomsAux,MynbasisAuxMPI
  integer,intent(in) :: nAtomsMPI(numnodes),startAuxMPI(nAtomsAux,numnodes)
  integer,intent(in) :: nAuxMPI(nAtomsAux,numnodes),nbasisAux
  real(realk),intent(in) :: alphaCD3(MynbasisAuxMPI,nvirt*nocc)
  real(realk),intent(inout) :: alphaCDFull(nbasisAux,nvirt*nocc)
  !
  integer :: iatomB,startB,IA,BETA,startB2,nAuxLoc
  !$OMP PARALLEL DEFAULT(none) PRIVATE(iatomB,startB,IA,&
  !$OMP BETA,startB2,nAuxLoc) SHARED(mynum,nAtomsMPI,startAuxMPI,&
  !$OMP nocc,nvirt,nAuxMPI,alphaCDFull,alphaCD3)
  startB2 = 0
  DO iAtomB=1,nAtomsMPI(mynum+1)
     StartB = startAuxMPI(iAtomB,mynum+1)
     nAuxLoc = nAuxMPI(iAtomB,mynum+1)
     !$OMP DO
     do IA = 1,nocc*nvirt
        do BETA = 1,nAuxLoc
           alphaCDFull(startB + BETA,IA) = alphaCD3(startB2 + BETA,IA)
        enddo
     enddo
     !$OMP END DO
     startB2 = startB2 + nAuxMPI(iAtomB,mynum+1)
  ENDDO
  !$OMP END PARALLEL
end subroutine PlugInToalphaCDFull

!alphaCD(NBA,nvirt,nocc) is in the diagonal basis 
subroutine RIMP2_calc_toccA(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,nvirt2,offset2)
  implicit none
  integer,intent(in) :: nvirt,nocc,noccEOS,NBA,nvirt2,offset2
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(in) :: EVocc(nocc),EVvirt(nvirt),UoccEOST(nocc,noccEOS)
  real(realk),intent(inout) :: tocc(nocc,noccEOS,nvirt,nvirt2)
  !
  integer :: BDIAG,ADIAG,IDIAG,JDIAG,ALPHAAUX,ILOC,JLOC
  real(realk) :: gmocont,deltaEPS,TMP
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
  real(realk) :: gmocont,deltaEPS,TMP,tvirtTMP(nvirt)
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
  real(realk) :: TMP
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
  real(realk) :: TMP
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
  real(realk) :: TMP
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

subroutine RIMP2_buildCalphaContFromAlphaCD(nocc,nvirt,myOriginalRank,numnodes,natoms,&
     & OriginalRanknbasisAuxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,AlphaCD5,&
     & Calpha,TMPAlphaBeta_inv,nbasisAux,MynbasisAuxMPI2)
  implicit none
  integer,intent(in) :: nocc,nvirt,myOriginalRank,numnodes,natoms,MynbasisAuxMPI2
  integer,intent(in) :: OriginalRanknbasisAuxMPI,nbasisAux
  integer,intent(in) :: nAtomsMPI(numnodes)
  integer,intent(in) :: startAuxMPI(nAtoms,numnodes)
  integer,intent(in) :: nAuxMPI(nAtoms,numnodes)
  real(realk),intent(in) :: AlphaCD5(OriginalRanknbasisAuxMPI,nvirt*nocc)
  real(realk),intent(inout) :: Calpha(MynbasisAuxMPI2,nvirt*nocc)
  real(realk),intent(in) :: TMPAlphaBeta_inv(MynbasisAuxMPI2,nbasisAux)
  !
  integer :: IB,startB2,iAtomB,StartB,BETA,ALPHA
  real(realk) :: TMP
  !Step 2: Obtain part of Calpha from this contribution
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(IB,startB2,iAtomB,startB,TMP,ALPHA,BETA) &
!$OMP SHARED(nocc,nvirt,myOriginalRank,numnodes,natoms,MynbasisAuxMPI2,&
!$OMP        OriginalRanknbasisAuxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,&
!$OMP        AlphaCD5,Calpha,TMPAlphaBeta_inv,nbasisAux)
  do IB = 1,nocc*nvirt
     startB2 = 0
     DO iAtomB=1,nAtomsMPI(myOriginalRank+1)
        StartB = startAuxMPI(iAtomB,myOriginalRank+1)
        do BETA = 1,nAuxMPI(iAtomB,myOriginalRank+1)
           TMP = AlphaCD5(startB2 + BETA,IB)
           do ALPHA = 1,MynbasisAuxMPI2
              Calpha(ALPHA,IB) = Calpha(ALPHA,IB) + TMPAlphaBeta_inv(ALPHA,startB + BETA)*TMP
           enddo
        enddo
        startB2 = startB2 + nAuxMPI(iAtomB,myOriginalRank+1)
     ENDDO
  enddo
!$OMP END PARALLEL DO
end subroutine RIMP2_buildCalphaContFromAlphaCD

subroutine  RIMP2_buildOwnCalphaFromAlphaCD(nocc,nvirt,mynum,numnodes,natoms,&
     & MynbasisAuxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,AlphaCD3,Calpha,&
     & TMPAlphaBeta_inv,nbasisAux,MynbasisAuxMPI2)
  implicit none
  integer,intent(in) :: nocc,nvirt,mynum,numnodes,natoms
  integer,intent(in) :: MynbasisAuxMPI,nbasisAux,MynbasisAuxMPI2
  integer,intent(in) :: nAtomsMPI(numnodes)
  integer,intent(in) :: startAuxMPI(nAtoms,numnodes)
  integer,intent(in) :: nAuxMPI(nAtoms,numnodes)
  real(realk),intent(in) :: AlphaCD3(MynbasisAuxMPI,nvirt*nocc)
  real(realk),intent(inout) :: Calpha(MynbasisAuxMPI2,nvirt*nocc)
  real(realk),intent(in) :: TMPAlphaBeta_inv(MynbasisAuxMPI2,nbasisAux)
  !
  integer :: IB,startB2,iAtomB,StartB,BETA,ALPHA
  real(realk) :: TMP

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(IB,startB2,iAtomB,StartB,BETA,ALPHA,TMP) &
!$OMP SHARED(nocc,nvirt,nAtomsMPI,startAuxMPI,mynum,MynbasisAuxMPI2,&
!$OMP        nAuxMPI,AlphaCD3,Calpha,TMPAlphaBeta_inv,MynbasisAuxMPI) 
  do IB = 1,nocc*nvirt
     startB2 = 0
     DO iAtomB=1,nAtomsMPI(mynum+1)
        StartB = startAuxMPI(iAtomB,mynum+1)
        do BETA = 1,nAuxMPI(iAtomB,mynum+1)
           TMP = AlphaCD3(startB2 + BETA,IB)
           do ALPHA = 1,MynbasisAuxMPI2
              Calpha(ALPHA,IB) = Calpha(ALPHA,IB) + &
                   & TMPAlphaBeta_inv(ALPHA,startB + BETA)*TMP
           enddo
        enddo
        startB2 = startB2 + nAuxMPI(iAtomB,mynum+1)
     ENDDO
  enddo
!$OMP END PARALLEL DO
end subroutine RIMP2_buildOwnCalphaFromAlphaCD

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
