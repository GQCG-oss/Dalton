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
  real(realk),pointer :: AlphaBeta(:,:),AlphaBeta_minus_sqrt(:,:)
  real(realk),pointer :: AlphaCD3(:,:,:),AlphaCD5(:,:,:),AlphaCD6(:,:,:)
  real(realk),pointer :: Calpha(:,:,:),Calpha2(:,:,:),Calpha3(:,:,:),CalphaVV(:,:,:)
  real(realk),pointer :: UoccEOS(:,:),UvirtEOS(:,:)
  real(realk),pointer :: gao(:,:,:,:),gmo(:,:,:,:),gocc(:,:,:,:),gocc2(:,:,:,:)
  real(realk),pointer :: tocc(:,:,:,:),UoccEOST(:,:),UvirtT(:,:),tocc3(:,:,:,:)
  real(realk),pointer :: toccTMP(:,:),TMPAlphaBeta_minus_sqrt(:,:),tocc2(:,:,:,:)
  real(realk),pointer :: tvirtTMP(:,:),tvirt(:,:,:,:),UoccT(:,:),UvirtEOST(:,:)
  real(realk),pointer :: tvirt2(:,:,:,:),tvirt3(:,:,:,:),Calpha4(:,:,:)
  real(realk),pointer :: UoccallT(:,:),CalphaOcc(:,:,:),tocc2TMP(:,:,:,:)
  real(realk) :: deltaEPS,goccAIBJ,goccBIAJ,Gtmp,Ttmp,Eocc,TMP,Etmp,twmpi2
  real(realk) :: gmocont,Gtmp1,Gtmp2,Eocc2,TMP1,flops,tmpidiff,EnergyMPI(2)
  real(realk) :: tcpu, twall,tcpu1,twall1,tcpu2,twall2,tcmpi1,tcmpi2,twmpi1
  real(realk) :: Evirt,Evirt2,dummy(2),MemInGBCollected
  integer(kind=long) :: maxsize
  Integer :: iAtomA,nBastLocA,startRegA,endRegA,nAuxA,startAuxA,endAuxA,lupri
  integer :: MynAtomsMPI,startA2,StartA,B,I,startB2,iAtomB,StartB,node,myOriginalRank
  Integer :: OriginalRanknbasisAuxMPI,NBA,dimocc(4),dimvirt(4)
  real(realk) :: time_i,time_c,time_w
  real(realk),pointer :: OccContribsFull(:),VirtContribsFull(:),Calpha_debug(:,:,:)
  real(realk),pointer :: occ_tmp(:),virt_tmp(:)
  integer,pointer :: IPVT(:)
  integer,pointer :: nbasisAuxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
  TYPE(MOLECULARORBITALINFO) :: orbitalInfo
  real(realk), pointer   :: work1(:),Etmp2222(:)
  real(realk)            :: RCOND
  integer(kind=ls_mpik)  :: COUNT,TAG,IERR,request,Receiver,sender,J,COUNT2
  real(realk) :: TS,TE,TS2,TE2,TS3,TE3
  real(realk) :: tcpu_start,twall_start, tcpu_end,twall_end,MemEstimate,memstep2
  integer ::CurrentWait(2),nAwaitDealloc,iAwaitDealloc,oldAORegular,oldAOdfAux
  integer :: MaxVirtSize,nTiles,offsetV
  logical :: useAlphaCD5,useAlphaCD6,ChangedDefault,first_order,PerformTiling
  integer(kind=ls_mpik)  :: request5,request6
  real(realk) :: phase_cntrs(nphases)
  integer(kind=long) :: nSize
  TYPE(MoleculeInfo),pointer      :: molecule1,molecule2,molecule3,molecule4
  ! cublas stuff
  type(c_ptr) :: cublas_handle
  integer*4 :: stat
  !> async handles
  integer :: num_ids
#ifdef VAR_OPENACC
  integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
  integer(kind=acc_device_kind) :: acc_device_type
#ifdef VAR_PGF90
  integer*4, external :: acc_set_cuda_stream
#endif
#else
  integer, pointer, dimension(:) :: async_id
#endif
  type(c_ptr) :: tocc_dev, tocc2_dev, tocc3_dev
#ifdef VAR_MPI
  INTEGER(kind=ls_mpik) :: HSTATUS
  CHARACTER*(MPI_MAX_PROCESSOR_NAME) ::  HNAME
  TAG = 1319
#endif  

  ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
  ! handle 1: UoccEOST, EVocc, and EVvirt
  ! handle 2: UvirtT
  ! handle 3: UvirtEOST
  ! handle 4: UoccT
  num_ids = 4
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
  nvirt = MyFragment%nunoccAOS     ! virtual AOS
  noccEOS = MyFragment%noccEOS     ! occupied EOS
  nvirtEOS = MyFragment%nunoccEOS  ! virtual EOS
  nocctot = MyFragment%nocctot     ! total occ: core+valence (identical to nocc without frozen core)
  ncore = MyFragment%ncore         ! number of core orbitals
  ! For frozen core energy calculation, we never need core orbitals
  ! (but we do if first order integrals are required)
  if(DECinfo%frozencore .and. (.not. first_order)) nocctot = nocc
  ! In general, for frozen core AND first order integrals, special care must be taken
  ! No frozen core OR frozen core calculation for just energy uses the same
  ! code from now on because the frozen core approximation is "built into" the fragment,
  if(DECinfo%frozencore .and. first_order) then
     fc = .true.
  else
     fc =.false.
  end if
  if(first_order) then
     if(DECinfo%PL>0) write(DECinfo%output,*) 'Calculating RIMP2 integrals (both energy and density) and RIMP2 amplitudes...'
  else
     if(DECinfo%PL>0) write(DECinfo%output,*) 'Calculating RIMP2 integrals (only energy) and RIMP2 amplitudes...'
  end if

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
  ! For frozen core energy calculation, we never need core orbitals
  ! (but we do if first order integrals are required)
  if(DECinfo%frozencore) nocctot = nocc

  if(master.AND.DECinfo%PL>0)THEN
     MemInGBCollected = 0.0E0_realk
     call get_currently_available_memory(MemInGBCollected)
     WRITE(DECinfo%output,'(1X,A)')'RIMP2MEM: RIMP2_integrals_and_amplitudes: Internal memory bookkeeping'
     WRITE(DECinfo%output,'(1X,A)')'RIMP2MEM: Memory Statistics at the beginning of the subroutine'
     write(DECinfo%output,'(1X,a,g12.4)') 'RIMP2MEM: Total memory:    ', DECinfo%memory
     WRITE(DECinfo%output,'(1X,a,g12.4)') 'RIMP2MEM: Memory Available:', MemInGBCollected

     !building Calpha 
#ifdef VAR_MPI
     MemEstimate = 2*nbasisAux/infpar%lg_nodtot*nvirt*nocc
     MemStep2 = nbasisAux/infpar%lg_nodtot*nvirt*nocc
#else
     MemEstimate = nbasisAux*nvirt*nocc
     MemStep2 = nbasisAux*nvirt*nocc
#endif
     !building toccEOS
     MemEstimate=MAX(MemEstimate,noccEOS*nocc*nvirt*nvirt+MemStep2)     
     MemEstimate=MAX(MemEstimate,noccEOS*noccEOS*nvirt*nvirt &
          & +nocc*nocc*nvirtEOS*nvirt+MemStep2)!building tvirtEOS + storage of toccEOS     
     MemEstimate=MAX(MemEstimate,2.0E0_realk*noccEOS*noccEOS*nvirt*nvirt&
          & +2.0E0_realk*nocc*nocc*nvirtEOS*nvirtEOS)  !storage of (toccEOS,tvirtEOS,goccEOS,gvirtEOS)
     IF(first_order) then !additional storage
        MemEstimate = MemEstimate + nvirt*noccEOS*noccEOS*nocctot + nvirtEOS*nocc*nvirtEOS*nvirt
     ENDIF
     MemEstimate = MemEstimate*mem_realsize
     WRITE(DECinfo%output,'(1X,a,g12.4)') 'RIMP2MEM: Estimated memory usage',MemEstimate
     if(DECinfo%PL>0)THEN
        call stats_globalmem(DECinfo%output)
     endif
  endif
  CALL LSTIMER('DECRIMP2: MEMcheck ',TS2,TE2,LUPRI,FORCEPRINT)

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
!$acc enter data copyin(EVocc,EVvirt) async(async_id(1))

  ! Extract occupied EOS indices from rows of Uocc
  call array2_extract_EOS(Uocc,MyFragment,'O','R',tmparray2)

  !make UoccEOS(noccEOS,nocc)
  call mem_alloc(UoccEOST,nocc,noccEOS) 
  M = noccEOS   !row of Input Matrix
  N = nocc      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,tmparray2%val,0.0E0_realk,UoccEOST)
!$acc enter data copyin(UoccEOST) async(async_id(1))
  call array2_free(tmparray2)

  call mem_alloc(UvirtT,nvirt,nvirt) 
  M = nvirt      !row of Input Matrix
  N = nvirt      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,Uvirt%val,0.0E0_realk,UvirtT)
!$acc enter data copyin(UvirtT) async(async_id(2))

  ! Extract virtual EOS indices from rows of Uvirt
  call array2_extract_EOS(Uvirt,MyFragment,'V','R',tmparray2)
  call array2_free(Uvirt)
  call mem_alloc(UvirtEOST,nvirt,nvirtEOS)
  M = nvirtEOS   !row of Input Matrix
  N = nvirt      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,tmparray2%val,0.0E0_realk,UvirtEOST)
!$acc enter data copyin(UvirtEOST) async(async_id(3))
  call array2_free(tmparray2)

  call mem_alloc(UoccT,nocc,nocc) 
  M = nocc      !row of Input Matrix
  N = nocc      !columns of Input Matrix
  call mat_transpose(M,N,1.0E0_realk,Uocc%val,0.0E0_realk,UoccT)
!$acc enter data copyin(UoccT) async(async_id(4))
  call array2_free(Uocc)

  if(fc) then
     call mem_alloc(UoccallT,nocctot,nocctot) 
     M = nocctot      !row of Input Matrix
     N = nocctot      !columns of Input Matrix
     call mat_transpose(M,N,1.0E0_realk,Uoccall%val,0.0E0_realk,UoccallT)
     call array2_free(Uoccall)
  endif
  CALL LSTIMER('DECRIMP2: TransMats ',TS2,TE2,LUPRI,FORCEPRINT)

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
  call Build_CalphaMO(MyFragment%mylsitem,master,nbasis,nbasisAux,LUPRI,&
       & FORCEPRINT,CollaborateWithSlaves,CDIAGocc%val,nocc,&
       & CDIAGvirt%val,nvirt,mynum,numnodes,nAtomsAux,Calpha,NBA)
  CALL LSTIMER('DECRIMP2: CalphaMO',TS2,TE2,LUPRI,FORCEPRINT)
  IF(.NOT.first_order)call array2_free(CDIAGvirt)
  IF(.NOT.first_order)call array2_free(CDIAGocc)
  !At this point we have the Calpha in the diagonal basis 

  !=====================================================================================
  !  Major Step 5: Generate toccEOS(nvirt,noccEOS,nvirt,noccEOS)
  !=====================================================================================

  dimocc = [nvirt,noccEOS,nvirt,noccEOS]   ! Output order
  IF(NBA.GT.0)THEN
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
     !Perform tiling if the tocc(nocc,noccEOS,nvirt,nvirt) cannot fit on the device
     !Janus will set this variable correctly. For now I set it true if it is a debug
     !run and false for release run. In this way all the code is tested
     PerformTiling = .FALSE.
     if (DECinfo%RIMP2_tiling) PerformTiling = .TRUE.
     MaxVirtSize = nvirt/2           !should be determined by Janus is some way
     nTiles =  nvirt/MaxVirtSize 
     IF(nTiles.EQ.0)PerformTiling = .FALSE.
     IF(PerformTiling)THEN
        call mem_alloc(tocc2,noccEOS,noccEOS,nvirt,nvirt)
        call mem_alloc(tocc,nocc,noccEOS,nvirt,MaxVirtSize)
        call mem_alloc(tocc2TMP,noccEOS,noccEOS,nvirt,MaxVirtSize)
!$acc enter data create(tocc,tocc2,tocc2TMP) copyin(Calpha) if(PerformTiling)
! here: wait for UoccEOST, EVocc, and EVvirt on async handle 1
!$acc wait(async_id(1))
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
!$acc exit data delete(tocc,tocc2TMP) if(PerformTiling)
        call mem_dealloc(tocc)
        call mem_dealloc(tocc2TMP)
     ELSE
        !Calculate and partial transform to local basis:
        !transform 1 occupied indices (IDIAG,JLOC,ADIAG,BDIAG)
        offsetV=0
        call mem_alloc(tocc,nocc,noccEOS,nvirt,nvirt)
!$acc enter data create(tocc) copyin(Calpha) if(.not. PerformTiling)

! here: wait for UoccEOST, EVocc, and EVvirt on async handle 1
!$acc wait(async_id(1))
        call RIMP2_calc_toccA(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST,nvirt,offsetV)
        !Transform second occupied index (IDIAG,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BDIAG)
        M = noccEOS              !rows of Output Matrix
        N = noccEOS*nvirt*nvirt  !columns of Output Matrix
        K = nocc                 !summation dimension
        call mem_alloc(tocc2,noccEOS,noccEOS,nvirt,nvirt)
!$acc enter data create(tocc2) if(.not. PerformTiling)
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
!$acc exit data delete(tocc) if(.not. PerformTiling)
#else
        call dgemm('T','N',M,N,K,1.0E0_realk,UoccEOST,K,tocc,K,0.0E0_realk,tocc2,M)
#endif
        call mem_dealloc(tocc)
     ENDIF
     !Transform first Virtual index (ILOC,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BLOC)
     M = noccEOS*noccEOS*nvirt  !rows of Output Matrix
     N = nvirt                  !columns of Output Matrix
     K = nvirt                  !summation dimension
     call mem_alloc(tocc3,nvirt,nvirt,noccEOS,noccEOS)
!$acc enter data create(tocc3)
#ifdef VAR_OPENACC
! here: wait for UvirtT on async handle 2
!$acc wait(async_id(2))
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
#else
     call dgemm('N','N',M,N,K,1.0E0_realk,tocc2,M,UvirtT,K,0.0E0_realk,tocc3,M)
#endif
     call mem_dealloc(tocc2)

     !Final virtual transformation and reorder to dimocc
     call tensor_ainit(toccEOS,dimocc,4)
!$acc enter data create(toccEOS%elm1)
     call RIMP2_calc_toccB(nvirt,noccEOS,tocc3,UvirtT,toccEOS%elm1)
!$acc exit data copyout(toccEOS%elm1) delete(tocc3)
     call mem_dealloc(tocc3)     
     CALL LSTIMER('RIMP2: toccEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     call tensor_ainit(toccEOS,dimocc,4)
     nsize = nvirt*noccEOS*nvirt*noccEOS
     call ls_dzero8(toccEOS%elm1,nsize)
  ENDIF
  CALL LSTIMER('DECRIMP2: tocc          ',TS2,TE2,LUPRI,FORCEPRINT)

#ifdef VAR_MPI
  IF(CollaborateWithSlaves) then  
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )       
     nSize = nvirt*noccEOS*nvirt*noccEOS
     call lsmpi_reduction(toccEOS%elm1,nsize,infpar%master,infpar%lg_comm)
     call time_start_phase(PHASE_WORK)   
     IF(.NOT.Master )call tensor_free(toccEOS)
  ENDIF
#endif

  !=====================================================================================
  !  Major Step 6: Generate tvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
  !=====================================================================================
  dimvirt = [nvirtEOS,nocc,nvirtEOS,nocc]   ! Output order
  IF(NBA.GT.0)THEN
     !Calculate and partial transform to local basis - transform occupied indices
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
     call mem_alloc(tvirt,nocc,nocc,nvirtEOS,nvirt) !IDIAG,JDIAG,ALOC,BDIAG
!$acc enter data create(tvirt)
! here: wait for UvirtEOST on async handle 3
!$acc wait(async_id(3))
     call RIMP2_calc_tvirtA(nvirt,nocc,nvirtEOS,NBA,Calpha,EVocc,EVvirt,tvirt,UvirtEOST)
!$acc exit data delete(EVocc,EVvirt)
     call mem_dealloc(EVocc)
     call mem_dealloc(EVvirt)

     !Transform first Virtual index (IDIAG,JDIAG,ALOC,BDIAG) => (IDIAG,JDIAG,ALOC,BLOC)
     M = nocc*nocc*nvirtEOS     !rows of Output Matrix
     N = nvirtEOS               !columns of Output Matrix
     K = nvirt                  !summation dimension
     call mem_alloc(tvirt2,nocc,nocc,nvirtEOS,nvirtEOS)
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
#else
     call dgemm('N','N',M,N,K,1.0E0_realk,tvirt,M,UvirtEOST,K,0.0E0_realk,tvirt2,M)
#endif
     call mem_dealloc(tvirt)

     !Transform first occupied index (IDIAG,JDIAG,ALOC,BLOC) => (ILOC,JDIAG,ALOC,BLOC)
     M = nocc                    !rows of Output Matrix
     N = nocc*nvirtEOS*nvirtEOS  !columns of Output Matrix
     K = nocc                    !summation dimension
     call mem_alloc(tvirt3,nocc,nocc,nvirtEOS,nvirtEOS)
!$acc enter data create(tvirt3)
#ifdef VAR_OPENACC
! here: wait for UoccT on async handle 4
!$acc wait(async_id(4))
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
#else
     call dgemm('T','N',M,N,K,1.0E0_realk,UoccT,K,tvirt2,M,0.0E0_realk,tvirt3,M)
#endif
     call mem_dealloc(tvirt2)

     !transform last occ index to local basis and reorder 
     call tensor_ainit(tvirtEOS,dimvirt,4)
!$acc enter data create(tvirtEOS%elm1)
     call RIMP2_calc_tvirtB(nvirtEOS,nocc,tvirt3,UoccT,tvirtEOS%elm1)
!$acc exit data delete(tvirt3) copyout(tvirtEOS%elm1)
     call mem_dealloc(tvirt3)
     CALL LSTIMER('RIMP2: tvirtEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     call tensor_ainit(tvirtEOS,dimvirt,4)
     call mem_dealloc(EVocc)
     call mem_dealloc(EVvirt)
     nSize = nvirtEOS*nocc*nvirtEOS*nocc
     call ls_dzero8(tvirtEOS%elm1,nSize)
  ENDIF

#ifdef VAR_MPI
  IF(CollaborateWithSlaves) then         
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     nSize = nvirtEOS*nocc*nvirtEOS*nocc
     call lsmpi_reduction(tvirtEOS%elm1,nSize,infpar%master,infpar%lg_comm)
     call time_start_phase(PHASE_WORK)   
     IF(.NOT.Master )call tensor_free(tvirtEOS)
  ENDIF
#endif

  !=====================================================================================
  !  Major Step 7: Generate goccEOS(nvirt,noccEOS,nvirt,noccEOS)
  !=====================================================================================

  IF(NBA.GT.0)THEN
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
     ! Transform Calpha(ALPHA,a,i) to local occupied index and local Virt
     ! Transform index delta to local occupied index 
     !(alphaAux;gamma,Jloc) = (alphaAux;gamma,J)*U(J,Jloc)     UoccEOST(iDIAG,iLOC)
     M = nba*nvirt  !rows of Output Matrix
     N = noccEOS          !columns of Output Matrix
     K = nocc             !summation dimension
     call mem_alloc(Calpha2,nba,nvirt,noccEOS)
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
#else
     call dgemm('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccEOST,K,0.0E0_realk,Calpha2,M)
#endif

     IF(.NOT.first_order)call mem_dealloc(UoccEOST)



     call mem_alloc(Calpha3,nba,nvirt,noccEOS)
!$acc enter data create(Calpha3)
     call RIMP2_TransAlpha1(nvirt,noccEOS,nba,UvirtT,Calpha2,Calpha3)
!$acc exit data delete(Calpha2,UvirtT) if(.not. first_order)
!$acc exit data delete(Calpha2) if(first_order)
     call mem_dealloc(Calpha2)
     IF(.NOT.first_order)call mem_dealloc(UvirtT)
     
     call tensor_ainit(goccEOS,dimocc,4)
!$acc enter data create(goccEOS%elm1) 
     call RIMP2_calc_gocc(nvirt,noccEOS,NBA,Calpha3,goccEOS%elm1)
!$acc exit data delete(Calpha3) copyout(goccEOS%elm1)
     call mem_dealloc(Calpha3)
     CALL LSTIMER('RIMP2: goccEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     IF(.NOT.first_order)call mem_dealloc(UoccEOST)
     IF(.NOT.first_order)call mem_dealloc(UvirtT)
     call tensor_ainit(goccEOS,dimocc,4)
     nSize = nvirt*noccEOS*nvirt*noccEOS
     call ls_dzero8(goccEOS%elm1,nsize)
  ENDIF
#ifdef VAR_MPI
  IF(CollaborateWithSlaves) then         
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     nSize = nvirt*noccEOS*nvirt*noccEOS
     call lsmpi_reduction(goccEOS%elm1,nSize,infpar%master,infpar%lg_comm)
     call time_start_phase(PHASE_WORK)   
     IF(.NOT.Master )call tensor_free(goccEOS)
  ENDIF
#endif

  !=====================================================================================
  !  Major Step 8: Generate gvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
  !=====================================================================================

  IF(NBA.GT.0)THEN
     CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
     ! Transform index delta to local occupied index 
     !(alphaAux;gamma,Jloc) = (alphaAux;gamma,J)*U(J,Jloc)     UoccEOST(iDIAG,iLOC)
     M = nba*nvirt  !rows of Output Matrix
     N = nocc             !columns of Output Matrix
     K = nocc             !summation dimension
     call mem_alloc(Calpha2,nba,nvirt,nocc)
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
#else
     call dgemm('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccT,K,0.0E0_realk,Calpha2,M)
#endif

     IF(.NOT.first_order)call mem_dealloc(UoccT)
     IF(.NOT.first_order)call mem_dealloc(Calpha)
     call mem_alloc(Calpha3,nba,nvirtEOS,nocc)

!$acc enter data create(Calpha3)
     call RIMP2_TransAlpha2(nocc,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha3)
!$acc exit data delete(Calpha2,UvirtEOST)
     IF(.NOT.first_order)call mem_dealloc(UvirtEOST)
     call mem_dealloc(Calpha2)

     call tensor_ainit(gvirtEOS,dimvirt,4)
!$acc enter data create(gvirtEOS%elm1)
     call RIMP2_calc_gvirt(nvirtEOS,nocc,NBA,Calpha3,gvirtEOS%elm1)
!$acc exit data delete(Calpha3) copyout(gvirtEOS%elm1)
     call mem_dealloc(Calpha3)
     CALL LSTIMER('RIMP2: gvirtEOS',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     IF(.NOT.first_order)call mem_dealloc(UvirtEOST)
     IF(.NOT.first_order)call mem_dealloc(UoccT)
     call tensor_ainit(gvirtEOS,dimvirt,4)
     nSize = nvirtEOS*nocc*nvirtEOS*nocc
     call ls_dzero8(gvirtEOS%elm1,nsize)
  ENDIF
#ifdef VAR_MPI
  IF(CollaborateWithSlaves) then         
     call time_start_phase( PHASE_IDLE )
     call lsmpi_barrier(infpar%lg_comm)
     call time_start_phase( PHASE_COMM )
     nSize = nvirtEOS*nocc*nvirtEOS*nocc
     call lsmpi_reduction(gvirtEOS%elm1,nsize,infpar%master,infpar%lg_comm)
     call time_start_phase( PHASE_WORK )   
     IF(.NOT.Master )call tensor_free(gvirtEOS)
  ENDIF
#endif

  IF(first_order)THEN

     !=====================================================================================
     !  first_order prop: Generate djik(nvirtAOS,noccEOS,noccEOS,noccAOS=nocctot)
     !=====================================================================================
     dimvirt = [nvirt,noccEOS,noccEOS,nocctot]   ! Output order
     IF(NBA.GT.0)THEN
        CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
        !(alphaAux;nvirt,JnoccEOS) = (alphaAux;nvirt,J)*U(J,JnoccEOS)
        M = nba*nvirt        !rows of Output Matrix
        N = noccEOS          !columns of Output Matrix
        K = nocc             !summation dimension
        call mem_alloc(Calpha2,nba,nvirt,noccEOS)
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
        call mem_alloc(Calpha3,nba,nvirt,noccEOS)
!$acc enter data create(Calpha3)
        call RIMP2_TransAlpha2(noccEOS,nvirt,nvirt,nba,UvirtT,Calpha2,Calpha3)
!$acc exit data delete(Calpha2)
        call mem_dealloc(Calpha2)
        
        CALL LSTIMER('START ',TS2,TE2,LUPRI)
        IF(DECinfo%frozencore)THEN
           call Build_CalphaMO(MyFragment%MyLsitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
                & CollaborateWithSlaves,CDIAGocc%val,nocc,CDIAGoccTALL%val,nocctot,mynum,&
                & numnodes,nAtomsAux,CalphaOcc,NBA)
           call array2_free(CDIAGoccTALL)
        ELSE
           call Build_CalphaMO(MyFragment%MyLsitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
                & CollaborateWithSlaves,CDIAGocc%val,nocc,CDIAGocc%val,nocc,mynum,&
                & numnodes,nAtomsAux,CalphaOcc,NBA)
           IF(nocctot.NE.nocc)call lsquit('FC Error RIMP2.',-1)
        ENDIF
        CALL LSTIMER('DECRIMP2: CalphaOO',TS2,TE2,LUPRI,FORCEPRINT)
        call array2_free(CDIAGocc)

        !(alphaAux;nocc,noccAOS=nocctot) = (alphaAux;nocc,nocc)*UoccallT(nocctot,nocctot)
        M = nba*nocc         !rows of Output Matrix
        N = nocctot          !columns of Output Matrix
        K = nocctot          !summation dimension
        call mem_alloc(Calpha2,nba,nocctot,nocctot)
        IF(DECinfo%frozencore)THEN
           call dgemm('N','N',M,N,K,1.0E0_realk,CalphaOcc,M,UoccallT,K,0.0E0_realk,Calpha2,M)
           call mem_dealloc(UoccallT)
        ELSE
           call dgemm('N','N',M,N,K,1.0E0_realk,CalphaOcc,M,UoccT,K,0.0E0_realk,Calpha2,M)
        ENDIF
        call mem_dealloc(CalphaOcc)
        !(alphaAux,noccEOS,noccAOS=nocctot) = (alphaAux;nocc,noccAOS=nocctot)*UoccEOST(nocc,noccEOS)
        call mem_alloc(Calpha4,nba,noccEOS,nocctot)
        call RIMP2_TransAlpha2(nocctot,nocc,noccEOS,nba,UoccEOST,Calpha2,Calpha4)
!$acc exit data delete(UoccEOST,Calpha2)
        call mem_dealloc(UoccEOST)
        call mem_dealloc(Calpha2)
        
        !  djikEOS(nvirtAOS,noccEOS,noccEOS,noccAOS)
        call tensor_ainit(djik,dimvirt,4)
        call RIMP2_calc_gen4DimFO(NBA,Calpha3,nvirt,noccEOS,Calpha4,noccEOS,nocctot,djik%elm1)
        call mem_dealloc(Calpha3)
        call mem_dealloc(Calpha4)
        CALL LSTIMER('RIMP2: djik',TS3,TE3,LUPRI,FORCEPRINT)
     ELSE
        call array2_free(CDIAGocc)
        call mem_dealloc(UoccEOST)
        call tensor_ainit(djik,dimvirt,4)
        nSize = nvirt*noccEOS*noccEOS*nocctot
        call ls_dzero8(djik%elm1,nsize)
     ENDIF
#ifdef VAR_MPI
     IF(CollaborateWithSlaves) then         
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        nSize = nvirt*noccEOS*noccEOS*nocctot
        call lsmpi_reduction(djik%elm1,nsize,infpar%master,infpar%lg_comm)
        call time_start_phase( PHASE_WORK )   
        IF(.NOT.Master )call tensor_free(djik)
     ENDIF
#endif

     !=====================================================================================
     !  first_order prop: Generate blad(nvirtEOS,noccAOS,nvirtEOS,nvirtAOS)
     !=====================================================================================  

     dimvirt = [nvirtEOS,nocc,nvirtEOS,nvirt]   ! Output order
     IF(NBA.GT.0)THEN
        CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
        CALL LSTIMER('START ',TS2,TE2,LUPRI)
        call Build_CalphaMO(MyFragment%MyLsitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
             & CollaborateWithSlaves,CDIAGvirt%val,nvirt,CDIAGvirt%val,nvirt,mynum,&
             & numnodes,nAtomsAux,CalphaVV,NBA)
        CALL LSTIMER('DECRIMP2: CalphaVV',TS2,TE2,LUPRI,FORCEPRINT)
        call mem_dealloc(CDIAGvirt%val)

        !(alphaAux;nvirt,noccAOS) = (alphaAux;nvirt,nocc)*U(nocc,noccAOS)
        M = nba*nvirt        !rows of Output Matrix
        N = nocc             !columns of Output Matrix
        K = nocc             !summation dimension
        call mem_alloc(Calpha2,nba,nvirt,nocc)
        call dgemm('N','N',M,N,K,1.0E0_realk,Calpha,M,UoccT,K,0.0E0_realk,Calpha2,M)
        call mem_dealloc(UoccT)
        
        !(alphaAux,nvirtEOS,noccAOS) = (alphaAux;nvirt,noccAOS)*Uvirt(nvirt,nvirtEOS)
        call mem_alloc(Calpha3,nba,nvirtEOS,nocc)
        call RIMP2_TransAlpha2(nocc,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha3)
        call mem_dealloc(Calpha2)
   
        !(alphaAux;nvirt,nvirtAOS) = (alphaAux;nvirt,nvirt)*UvirtT(nvirt,nvirt)
        M = nba*nvirt        !rows of Output Matrix
        N = nvirt            !columns of Output Matrix
        K = nocc             !summation dimension
        call mem_alloc(Calpha2,nba,nvirt,nvirt)
        call dgemm('N','N',M,N,K,1.0E0_realk,CalphaVV,M,UvirtT,K,0.0E0_realk,Calpha2,M)
        call mem_dealloc(CalphaVV)
        call mem_dealloc(UvirtT)
        call mem_dealloc(Calpha)
        
        !(alphaAux,nvirtEOS,nvirtAOS) = (alphaAux;nvirt,nvirtAOS)*UvirtEOST(nvirt,nvirtEOS)
        call mem_alloc(Calpha4,nba,nvirtEOS,nvirt)
        call RIMP2_TransAlpha2(nvirt,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha4)
        call mem_dealloc(UvirtEOST)
        call mem_dealloc(Calpha2)
        
        !generate blad(nvirtEOS,noccAOS,nvirtEOS,nvirtAOS)
        call tensor_ainit(blad,dimvirt,4)
        call RIMP2_calc_gen4DimFO(NBA,Calpha3,nvirtEOS,nocc,Calpha4,nvirtEOS,nvirt,blad%elm1)
        call mem_dealloc(Calpha3)
        call mem_dealloc(Calpha4)
        CALL LSTIMER('RIMP2: blad',TS3,TE3,LUPRI,FORCEPRINT)
     ELSE
        call mem_dealloc(CDIAGvirt%val)
        call mem_dealloc(UoccT)
        call mem_dealloc(UvirtT)
        call mem_dealloc(UvirtEOST)
        call tensor_ainit(djik,dimvirt,4)
        nSize = nvirtEOS*nocc*nvirtEOS*nvirt
        call ls_dzero8(djik%elm1,nsize)
     ENDIF
#ifdef VAR_MPI
     IF(CollaborateWithSlaves) then         
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        nSize = nvirtEOS*nocc*nvirtEOS*nvirt
        call lsmpi_reduction(blad%elm1,nsize,infpar%master,infpar%lg_comm)
        call time_start_phase( PHASE_WORK )   
        IF(.NOT.Master )call tensor_free(blad)
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
     call end_flop_counter(flops)
  end if

#ifdef VAR_MPI
  ! If slaves were not invoked
  ! then we of course skip the reduction.
  MPIcollect: if(wakeslave) then
     ! FLOP counting
     if(master) then
        flops=0.0E0_realk  ! we want to count only flops from slaves (these were set above)
        ! Total time for all slaves (not local master itself)
        MyFragment%slavetime_work(MODEL_RIMP2)=0.0E0_realk
     end if
     if(master)MyFragment%flops_slaves=flops !save flops for local slaves (not local master)
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
     FLUSH(DECinfo%output)
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
     & CollaborateWithSlaves,Cocc,nocc,Cvirt,nvirt,mynum,&
     & numnodes,nAtomsAux,Calpha,NBA)
  implicit none
  type(lsitem), intent(inout) :: mylsitem
  integer,intent(inout) :: NBA
  integer,intent(in) :: nAtomsAux,nocc,nvirt
  integer,intent(in) :: nbasisAux,LUPRI,nbasis,mynum,numnodes
  logical,intent(in) :: master,FORCEPRINT,CollaborateWithSlaves
  real(realk),intent(in) :: Cocc(nbasis,nocc),Cvirt(nbasis,nvirt)
  real(realk),pointer :: Calpha(:,:,:)
  !
  integer :: MynbasisAuxMPI
  integer,pointer :: nbasisAuxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
  real(realk),pointer :: AlphaBeta(:,:),AlphaBeta_minus_sqrt(:,:)
  real(realk),pointer :: TMPAlphaBeta_minus_sqrt(:,:),AlphaCD3(:,:,:)
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
     call II_get_RI_AlphaBeta_2centerInt(DECinfo%output,DECinfo%output,&
          & AlphaBeta,mylsitem%setting,nbasisAux)
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
     call mem_alloc(AlphaBeta_minus_sqrt,nbasisAux,nbasisAux)
     call lowdin_diag_S_minus_sqrt(nbasisAux, AlphaBeta,AlphaBeta_minus_sqrt, lupri)
     call mem_dealloc(AlphaBeta)
     CALL LSTIMER('AlphaBetamSq ',TS3,TE3,LUPRI,FORCEPRINT)
  ELSE
     call mem_alloc(AlphaBeta_minus_sqrt,nbasisAux,nbasisAux)
  ENDIF
#ifdef VAR_MPI
  call time_start_phase( PHASE_IDLE )
  call lsmpi_barrier(infpar%lg_comm)
  call time_start_phase( PHASE_COMM )
  call ls_mpibcast(AlphaBeta_minus_sqrt,nbasisAux,nbasisAux,infpar%master,infpar%lg_comm)
  call time_start_phase(PHASE_WORK)   
#endif

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
     call mem_alloc(TMPAlphaBeta_minus_sqrt,MynbasisAuxMPI2,nbasisAux)
     offset = mynum*ndimMax
     offset2 = numnodes*ndimMax + mynum -1 +1
     IF(MynbasisAuxMPI2.GT.ndimMax)THEN
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(nbasisAux,ndimMax,&
        !$OMP TMPAlphaBeta_minus_sqrt,AlphaBeta_minus_sqrt,offset,offset2)
        do I=1,nbasisAux
           do J=1,ndimMax
              TMPAlphaBeta_minus_sqrt(J,I) = AlphaBeta_minus_sqrt(offset+J,I)
           enddo
           TMPAlphaBeta_minus_sqrt(ndimMax+1,I) = AlphaBeta_minus_sqrt(offset2,I)
        enddo
        !$OMP END PARALLEL DO
     ELSE
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(nbasisAux,ndimMax,&
        !$OMP TMPAlphaBeta_minus_sqrt,AlphaBeta_minus_sqrt,offset)
        do I=1,nbasisAux
           do J=1,ndimMax
              TMPAlphaBeta_minus_sqrt(J,I) = AlphaBeta_minus_sqrt(offset+J,I)
           enddo
        enddo
        !$OMP END PARALLEL DO
     ENDIF
     call mem_dealloc(AlphaBeta_minus_sqrt)
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
     maxsize = NINT(MemInGBCollected*1.E9_realk)
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
     call II_get_RI_AlphaCD_3centerInt2(DECinfo%output,DECinfo%output,&
          & AlphaCD3,mylsitem%setting,nbasisAux,nbasis,&
          & nvirt,nocc,Cvirt,Cocc,maxsize,mynum,numnodes)
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
        !Calpha = TMPAlphaBeta_minus_sqrt(MynbasisAuxMPI,nbasisAux)
        M =  MynbasisAuxMPI2   !rows of Output Matrix
        N =  nvirt*nocc       !columns of Output Matrix
        K =  nbasisAux        !summation dimension
        call mem_alloc(Calpha,MynbasisAuxMPI2,nvirt,nocc)
        call dgemm('N','N',M,N,K,1.0E0_realk,TMPAlphaBeta_minus_sqrt,&
             & M,alphaCDFull,K,0.0E0_realk,Calpha,M)
        call mem_dealloc(alphaCDFull)
        call mem_dealloc(TMPAlphaBeta_minus_sqrt)
     ELSE
        !Serial version
        M =  MynbasisAuxMPI   !rows of Output Matrix
        N =  nvirt*nocc       !columns of Output Matrix
        K =  nbasisAux        !summation dimension
        call mem_alloc(Calpha,MynbasisAuxMPI,nvirt,nocc)
        call dgemm('N','N',M,N,K,1.0E0_realk,AlphaBeta_minus_sqrt,&
             & M,AlphaCD3,K,0.0E0_realk,Calpha,M)
        call mem_dealloc(AlphaCD3)
        call mem_dealloc(AlphaBeta_minus_sqrt)
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
                & AlphaCD3,Calpha,TMPAlphaBeta_minus_sqrt,nbasisAux,&
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
                & Calpha,TMPAlphaBeta_minus_sqrt,nbasisAux,MynbasisAuxMPI2)
           call mem_dealloc(AlphaCD5)
        ENDIF
     ENDDO
     call mem_dealloc(TMPAlphaBeta_minus_sqrt)
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
  !acc present(tocc,Calpha,UoccEOST,EVocc,EVvirt)
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

subroutine RIMP2_calc_gvirt(nvirtEOS,nocc,NBA,Calpha3,gvirtEOS)
  implicit none
  integer,intent(in) :: nvirtEOS,nocc,NBA
  real(realk),intent(in) :: Calpha3(NBA,nvirtEOS,nocc)
  real(realk),intent(inout) :: gvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
  !local variables
  integer :: BLOC,JLOC,ILOC,ALOC,ALPHAAUX
  real(realk) :: TMP
#ifdef VAR_OPENACC
  !$ACC PARALLEL LOOP COLLAPSE(4) &
  !$ACC PRIVATE(BLOC,JLOC,ILOC,ALOC,ALPHAAUX,TMP) &
  !$ACC COPYIN(nvirtEOS,nocc,NBA) &
  !$acc present(Calpha3,gvirtEOS)
#else
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(BLOC,JLOC,ILOC,ALOC,ALPHAAUX,TMP) &
  !$OMP SHARED(nvirtEOS,nocc,NBA,Calpha3,gvirtEOS)
#endif
  do jLOC=1,nocc
     do bLOC=1,nvirtEOS
        do iLOC=1,nocc
           do aLOC=1,nvirtEOS
              TMP = 0.0E0_realk
#ifdef VAR_OPENACC
              !$ACC loop seq
#endif 
              do ALPHAAUX = 1,nba
                 tmp = tmp + Calpha3(alphaAUX,ALOC,ILOC)*Calpha3(alphaAUX,BLOC,JLOC) 
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
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(A,B,C,D,ALPHAAUX,TMP) &
  !$OMP SHARED(n1,n2,n3,n4,NBA,Calpha3,Calpha4,djik)
  do d=1,n4
     do c=1,n3
        do b=1,n2
           do a=1,n1
              TMP = 0.0E0_realk
              do ALPHAAUX = 1,nba
                 tmp = tmp + Calpha3(alphaAUX,A,B)*Calpha4(alphaAUX,C,D)
              enddo
              djik(A,B,C,D) = tmp
           enddo
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
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
