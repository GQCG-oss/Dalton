!> @file
!> Resolution of the identity (density-fitting) tools
!> \ author Thomas Kjaergaard

module ri_util_module

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
  use iso_c_binding
#ifdef VAR_OPENACC
  use openacc
#endif
#if defined(VAR_CUDA) || defined(VAR_OPENACC)
  use gpu_interfaces
#endif
  ! DEC DEPENDENCIES (within deccc directory)
  ! *****************************************
  use dec_fragment_utils

contains
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
  character,intent(in) :: intspec(5) 
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
  integer :: MetricOper
  real(realk) :: epsilon
  logical :: use_bg_buf,noOMP
#ifdef VAR_OMP
  integer, external :: OMP_GET_MAX_THREADS
#endif
  CALL LSTIMER('START ',TS,TE,LUPRI,ForcePrint)
  CALL LSTIMER('START ',TS3,TE3,LUPRI,ForcePrint)
  noOMP = mylsitem%setting%scheme%noOMP
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
  IF(use_bg_buf)THEN
     maxsize = mem_get_bg_buf_free()*8.0E-9_realk
  ELSE
     call get_currently_available_memory(MemInGBCollected)
     maxsize = MemInGBCollected*0.65E0_realk
  ENDIF
  call GetOperatorFromCharacter(MetricOper,intspec(5),mylsitem%Setting)

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
        IF(DECinfo%RIMP2_lowdin)THEN
           call mem_alloc(AlphaBeta,nbasisAux,nbasisAux)
        ENDIF
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
        IF(DECinfo%RIMP2_lowdin)THEN
           call II_get_RI_AlphaBeta_2centerInt(DECinfo%output,DECinfo%output,&
                & AlphaBeta,mylsitem%setting,nbasisAux,MetricOper)
        ELSE
           call II_get_RI_AlphaBeta_2centerInt(DECinfo%output,DECinfo%output,&
                & AlphaBetaDecomp,mylsitem%setting,nbasisAux,MetricOper)
        ENDIF
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
        IF(DECinfo%RIMP2_lowdin)THEN
           call lowdin_diag_S_minus_sqrt(nbasisAux,AlphaBeta,AlphaBetaDecomp,lupri)
           call mem_dealloc(AlphaBeta)
        ELSE
           call Get_InverseCholeskyFactor(nbasisAux,AlphaBetaDecomp,lupri)
        ENDIF
        CALL LSTIMER('DF_Calpha:AlphaBetaDecomp',TS4,TE4,LUPRI,ForcePrint)
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
  call GetOperatorFromCharacter(Oper,intspec(4),mylsitem%Setting)

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
           print*,'MemForPartialMOINT',MemForPartialMOINT
           print*,'maxsize',maxsize
           print*,'==================================='
           print*,'MinAuxBatch',MinAuxBatch
           print*,'nAuxMPI(mynum+1)',nAuxMPI(mynum+1)
           print*,'nvirt*nocc',nvirt*nocc,'nvirt,nocc',nvirt,nocc
           print*,'nthreads',nthreads
           print*,'nbasis1',nbasis1
           print*,'nbasis2',nbasis2
           IF(nthreads.GT.1)THEN
              !use Critical instead ! 
              !deactivating OpenMP 
              mylsitem%setting%scheme%noOMP = .TRUE.
              nthreads = 1              
           ELSE
              CALL lsquit('Not enough memory in build_calpha bcast schem',-1)
           ENDIF
        ENDIF

        MaxNaux = MIN(nAuxMPI(mynum+1),&
             &FLOOR((MaxSize/8.0E-9_realk-nAuxMPI(mynum+1)*nvirt*nocc-&
             & nbasis1*nvirt-nbasis2*nocc) &
             & /((nbasis1*nocc+nbasis1*nbasis2)*nthreads)))
        
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
       & AlphaCD3,mylsitem%setting,nbasisAux,nbasis1,nbasis2,intspec(1:4),MaxNaux,&
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
  mylsitem%setting%scheme%noOMP = noOMP
!  call sleep(mynum*5)
!  PRINT*,'Build_CalphaMO2:Nba',Nba
!  WRITE(6,*)'Build_CalphaMO2:Final Calph(NBA=',NBA,',nvirt=',nvirt,',nocc=',nocc,')'
!  WRITE(6,*)'Build_CalphaMO2:Print Subset Final Calph(NBA=',NBA,',1:4)  MYNUM',MYNUM
!  call ls_output(Calpha,1,NBA,1,4,NBA,nvirt*nocc,1,6)
!
  CALL LSTIMER('Build_CalphaMO2',TS,TE,LUPRI,ForcePrint)

end subroutine Build_CalphaMO2

!This should be call my master and slaves
subroutine Build_RobustERImatU(myLSitem,master,nbasis1,nbasis2,nbasisAux,&
     & LUPRI,FORCEPRINT,&
     & CollaborateWithSlaves,Cvirt,nvirt,Cocc,nocc,mynum,numnodes,&
     & AlphaBetaDecomp,intspec,Umat)
  implicit none
  type(lsitem), intent(inout) :: mylsitem
  integer,intent(in) :: nocc,nvirt
  integer,intent(in) :: nbasisAux,LUPRI,nbasis1,nbasis2,mynum,numnodes
  logical,intent(in) :: master,FORCEPRINT,CollaborateWithSlaves
  real(realk),intent(in) :: AlphaBetaDecomp(nbasisAux,nbasisAux)
  real(realk),intent(inout) :: Umat(nbasisAux,nbasisAux)
  real(realk),intent(in) :: Cvirt(nbasis1,nvirt),Cocc(nbasis2,nocc)
  character,intent(in) :: intspec
  !
  real(realk),pointer :: AlphaBeta(:,:),TMP(:,:)
  real(realk) :: TS3,TE3,MemInGBCollected,MemInGBCollected2,TS,TE
  real(realk) :: MemForFullAOINT,MemForFullMOINT,maxsize,MemForPartialMOINT
  real(realk) :: TS4,TE4
  TYPE(MoleculeInfo),pointer      :: molecule1,molecule2,molecule3,molecule4
  integer(kind=long)    :: nSize,n8
  integer(kind=ls_mpik) :: node 
  integer :: M,N,K,MetricOper
  integer :: J,offset
  CALL LSTIMER('START ',TS,TE,LUPRI,ForcePrint)
  print*,'Build_RobustERImatU intspec=',intspec
  call GetOperatorFromCharacter(MetricOper,intspec,mylsitem%Setting)

  !=====================================================================================
  ! Master Obtains (alpha|beta) ERI in Auxiliary Basis 
  !=====================================================================================
  IF(master)THEN
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
          & Umat,mylsitem%setting,nbasisAux,MetricOper)
     CALL LSTIMER('DF_Calpha:AlphaBeta',TS4,TE4,LUPRI,ForcePrint)
     IF(DECinfo%AuxAtomicExtent)THEN
        mylsitem%SETTING%MOLECULE(1)%p => molecule1
        mylsitem%SETTING%MOLECULE(2)%p => molecule2
        mylsitem%SETTING%MOLECULE(3)%p => molecule3
        mylsitem%SETTING%MOLECULE(4)%p => molecule4
     ENDIF
     N = nbasisAux
     call mem_alloc(TMP,N,N)
     call dgemm('N','T',N,N,N,1.0E0_realk,Umat,N,AlphaBetaDecomp,N,0.0E0_realk,TMP,N)
     call dgemm('N','N',N,N,N,1.0E0_realk,AlphaBetaDecomp,N,TMP,N,0.0E0_realk,Umat,N)
     call mem_dealloc(TMP)
  ENDIF
#ifdef VAR_MPI
  call time_start_phase( PHASE_IDLE )
  call lsmpi_barrier(infpar%lg_comm)
  call time_start_phase( PHASE_COMM )
  call ls_mpibcast(Umat,nbasisAux,nbasisAux,infpar%master,infpar%lg_comm)
  call time_start_phase(PHASE_WORK)   
#endif
  CALL LSTIMER('Build_RobustERImatU',TS,TE,LUPRI,ForcePrint)

end subroutine Build_RobustERImatU

!computes the Cholesky factorization of the inverse of 
!a real symmetric positive definite matrix A = U^T * U 
subroutine Get_InverseCholeskyFactor(n,A,lupri)
  implicit none
  integer, intent(in)          :: n,lupri
  real(realk), intent(inout)   :: A(n,n)
  !
  integer                      :: i,j,np,k,info
  real(realk) :: TS,TE
  call LSTIMER('START ',TS,TE,lupri)
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
  call LSTIMER('Get_CholeskyFactor',TS,TE,lupri)
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

end module ri_util_module

