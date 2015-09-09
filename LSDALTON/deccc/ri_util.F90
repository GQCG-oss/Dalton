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
  use lapackMod
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
  use background_buffer_module
#ifdef VAR_OPENACC
  use openacc
#endif
#if defined(VAR_CUDA) || defined(VAR_OPENACC)
  use gpu_interfaces
#endif  
  ! DEC DEPENDENCIES (within deccc directory)
  ! *****************************************
  use dec_fragment_utils
  private
  public :: Build_CalphaMO2,BuildnAuxMPIUsedRI,BuildnAuxMPIUsedRIinfo,&
       & Build_RobustERImatU, Build_RIMP2grad, BuilDUmatTmpRIF12

contains
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
  integer(kind=long)    :: nSize,n8,nbasisAux8,SaveMaxMemoryUsage
  integer(kind=ls_mpik) :: node 
  integer :: MaxNaux,M,N,K,ndimMax1,nbasisAuxMPI2(numnodes),MynbasisAuxMPI2
  integer :: nthreads,PerformReduction,Oper,nAuxMPI(numnodes),MaxnAuxMPI,I
  integer :: J,offset,offset2,inode,nbuf1,nbuf2,nbuf3,dim1,MinAuxBatch,ndimMax2
  integer :: GindexToLocal(nbasisAux),nbasisAuxMPI3(numnodes),NREDLOC,NRED
  integer :: MetricOper,AuxDimUsedInAOcode
  real(realk) :: epsilon,MaxMemoryUsageGB,MemoryEstimateGB
  logical :: noOMP,use_bg_buf,noOMPsave
#ifdef VAR_OMP
  integer, external :: OMP_GET_MAX_THREADS
#endif
  CALL LSTIMER('START ',TS,TE,LUPRI,ForcePrint)
  CALL LSTIMER('START ',TS3,TE3,LUPRI,ForcePrint)
  noOMP = mylsitem%setting%scheme%noOMP  
  nbasisAux8 = nbasisAux   
  epsilon = DECinfo%NAFthreshold 
  use_bg_buf = .FALSE.
  IF(present(use_bg_bufInput)) use_bg_buf = use_bg_bufInput
  IF(use_bg_buf.AND.DECinfo%NAF)call lsquit('bg_buf and NAF combi not tested',-1)
  IF(DECinfo%PL.GT.2)THEN
   IF(DECinfo%NAF)THEN
    WRITE(DECinfo%output,*)'Use Natural Auxiliary Functions (NAF) with Threshold = ',&
         & epsilon
   ENDIF
  ENDIF
  PerformReduction = -1 !Initialization.
  !PerformReduction=0 means No Reduction
  !PerformReduction=X means Perform Reduction

  IF(DECinfo%MemDebugPrint)THEN
     WRITE(DECinfo%output,*)'RIUTILmeminfo: Memory Before RIMP2 '
     IF(use_bg_buf)THEN
        WRITE(DECinfo%output,*)'RIUTILmeminfo: Memory in Background Buffer=',mem_get_bg_buf_free()*8.0E-9_realk,' GB'
        call printBGinfo()     
     ELSE
        MemInGBCollected = 0.0E0_realk
        call get_currently_available_memory(MemInGBCollected)  
        WRITE(DECinfo%output,*)'RIUTILmeminfo: Memory available=',MemInGBCollected,' GB'
     ENDIF
  ENDIF

  NBA = 0
  MynbasisAuxMPI2 = 0 
  IF(use_bg_buf)THEN
     maxsize = mem_get_bg_buf_free()*8.0E-9_realk*0.90E0_realk
     IF(DECinfo%MemDebugPrint)THEN
        call printBGinfo()
     ENDIF     
  ELSE
     call get_currently_available_memory(MemInGBCollected)
     maxsize = MemInGBCollected*0.75E0_realk
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
     IF(DECinfo%RIMP2ForcePDMCalpha)THEN
        !Force use of Bcast scheme - used when full MO cannot fit in memory       
        PerformReduction = 0 !Perform Bcast scheme
        IF(DECinfo%MemDebugPrint)print*,'RIUTIL: Force Bcast scheme'
     ELSE
        !Memory accounting for Reduction of full MO integral alphaCD(nbasisAux*nvirt*nocc): 
        !allocate TMPAlphaBetaDecomp(MynbasisAuxMPI2,nbasisAux) NOT allocated in BG buffer
        !BG  : allocate Calpha(nbasisAux*nvirt*nocc) 
        !noBG: allocate Calpha(MynbasisAux*nvirt*nocc) 
        !allocate alphaCD(nbasisAux*nvirt*nocc)
        !allocate AO to MO tmp (MinAuxBatch*nthreads*(nbasis1*nocc+nbasis1*nbasis2))
        !Assuming internal AO to MO transform using MaxNaux = MinAuxBatch
        IF(use_bg_buf)THEN
           !only looking at the memory allocated in BG buffer 
           MemForFullMOINT = (2*nbasisAux*nvirt*nocc+MinAuxBatch*nthreads*&
                & (nbasis1*nocc+nbasis1*nbasis2))*8.0E-9_realk &
                & + (nbasisAux/numnodes+1)*nocc*nvirt*8.0E-9_realk        
        ELSE
           MemForFullMOINT = (2*nbasisAux*nvirt*nocc+nbasisAux*nvirt*nocc/(numnodes+2)&
                & +MinAuxBatch*nthreads*(nbasis1*nocc+nbasis1*nbasis2)&
                & +nbasis1*nvirt+nbasis2*nocc+nbasisAux*nbasisAux)*8.0E-9_realk 
        ENDIF
        IF(DECinfo%MemDebugPrint)THEN
           print*,'RIUTIL: Full MO (alpha|cd) integral requires ',MemForFullMOINT,' GB'
           print*,'RIUTIL: Full MEM1: ',2*nbasisAux*nvirt*nocc*8.0E-9_realk ,' GB'
           print*,'RIUTIL: Full MEM2: ',MinAuxBatch*nthreads*(nbasis1*nocc+nbasis1*nbasis2)*8.0E-9_realk ,' GB'
           print*,'RIUTIL: Full MEM3: ',(nbasisAux/numnodes+1)*nocc*nvirt*8.0E-9_realk,' GB'
           IF(use_bg_buf)THEN
              print*,'RIUTIL: Mem used : ',buf_realk%offset*8.0E-9_realk,' GB'
              print*,'RIUTIL: Resulting in Memory Estimate of  ',&
                   & (buf_realk%offset+MemForFullMOINT/8.0E-9_realk)*8.0E-9_realk,' GB'
              print*,'RIUTIL: Smaller than buffer size of      ',buf_realk%nmax*8.0E-9_realk,' GB'
           ENDIF
           print*,'RIUTIL: Memory available (75%/90% bg_buffer)',maxsize,' GB'
           IF(MemForFullMOINT.GE.maxsize)THEN
              print*,'RIUTIL: Full MO cannot fit in memory, we cannot do a simple reduction'
           ELSE
              print*,'RIUTIL: Full MO can fit in memory, we can do a simple reduction'
           ENDIF
        ENDIF
        IF(MemForFullMOINT.GE.maxsize)THEN
           !Full MO cannot fit in memory       
           print*,'MemForFullMOINT',MemForFullMOINT
           print*,'maxsize',maxsize
           print*,'MemForFullMOINT.GE.maxsize',MemForFullMOINT.GE.maxsize
           PerformReduction = 0 !Perform Bcast scheme
           IF(numnodes.LT.2)THEN
              call lsquit('Not enough memory to store Calpha(NBA,nvirt,nocc) increase memory or numnodes',-1)
           ENDiF
        ELSE
           PerformReduction = 1
        ENDIF
     ENDIF
  ENDIF
#ifdef VAR_MPI
  call time_start_phase( PHASE_IDLE )
  call lsmpi_barrier(infpar%lg_comm)
  call time_start_phase( PHASE_COMM )
  call ls_mpibcast(PerformReduction,infpar%master,infpar%lg_comm)  
  call time_start_phase( PHASE_WORK )
#endif
  CALL LSTIMER('DF_Calpha:Init ',TS3,TE3,LUPRI,ForcePrint)
  !=====================================================================================
  ! Master Obtains (alpha|beta) ERI in Auxiliary Basis 
  !=====================================================================================

  IF(AlphaBetaDecompCreate)THEN
     IF(master)THEN
        IF(DECinfo%RIMP2_lowdin)THEN
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: alloc AlphaBeta(',nbasisAux8*nbasisAux8,')'
              call mem_pseudo_alloc(AlphaBeta,nbasisAux8,nbasisAux8)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc AlphaBeta(',nbasisAux*nbasisAux,')'
              call mem_alloc(AlphaBeta,nbasisAux,nbasisAux,'RIUTIL:AlphaBeta')
           ENDIF
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
           IF(use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc AlphaBeta(',size(AlphaBeta),')'
              call mem_pseudo_dealloc(AlphaBeta)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc AlphaBeta(',size(AlphaBeta),')'
              call mem_dealloc(AlphaBeta)
           ENDIF
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
  noOMP = .FALSE.
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
     call BuildnAuxMPIUsedRI(nbasisAux,numnodes,nbasisAuxMPI2)
     ndimMax1 = nbasisAux/numnodes
     MynbasisAuxMPI2 = nbasisAuxMPI2(mynum+1)
     IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
     IF(DECinfo%MemDebugPrint)print*,'STD: alloc TMPAlphaBetaDecomp(',MynbasisAuxMPI2*nbasisAux,')'
     call mem_alloc(TMPAlphaBetaDecomp,MynbasisAuxMPI2,nbasisAux,'RIUTIL:TMPAlphaBetaDecomp')
     call buildTMPAlphaBetaDecomp(TMPAlphaBetaDecomp,AlphaBetaDecomp,&
          & MynbasisAuxMPI2,nbasisAux,mynum,ndimMax1,numnodes)
     NBA = MynbasisAuxMPI2

     IF(use_bg_buf)THEN
        !allocate now because I need to allocate in order due to push pop mechanisme
        nsize = NBA*nvirt*nocc     
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha(',nsize,')'
        call mem_pseudo_alloc(Calpha,nsize)
     ENDIF
     IF(PerformReduction.EQ.0)THEN 
        dim1 = nAuxMPI(mynum+1)    !DO NOT PERFORM REDUCTION each node have part of AlphaCD(nAuxMPI(mynum+1),nvirt,nocc)
        AuxDimUsedInAOcode = MaxnAuxMPI
     ELSE
        dim1 = nBasisAux           !Reduction all nodes have full AlphaCD allocated (partial filled in)
        AuxDimUsedInAOcode = dim1
     ENDIF

     !Determine internal or external AO to MO transform 
     IF(mylsitem%setting%scheme%ForceRIMP2memReduced)THEN 
        !Force Internal AO to MO transform 
        MaxNaux = MinAuxBatch + 1
     ELSE
        call DetermineMaxNauxRI(use_bg_buf,noOMP,dim1,AuxDimUsedInAOcode,nbasis1,nbasis2,&
             & nocc,nvirt,MinAuxBatch,numnodes,nAuxMPI,mynum,MaxNaux,nthreads)
        IF(noOMP)THEN
           noOMPsave= mylsitem%setting%scheme%noOMP
           mylsitem%setting%scheme%noOMP = .TRUE.
        ENDIF
     ENDIF
  ELSE
     NBA = nbasisAux
     AuxDimUsedInAOcode = nbasisAux
     nAuxMPI = nbasisAux
     !Serial Code
     IF(use_bg_buf)THEN
        !allocate now because I need to allocate in order due to push pop mechanisme
        nsize = NBA*nvirt*nocc     
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha(',nsize,')'
        call mem_pseudo_alloc(Calpha,nsize)
     ENDIF
     !determine Internal or external AO to MO transform. 
     IF(mylsitem%setting%scheme%ForceRIMP2memReduced)THEN 
        !Force Internal AO to MO transform 
        MaxNaux = MinAuxBatch + 1
     ELSE
        call DetermineMaxNauxRI(use_bg_buf,noOMP,dim1,AuxDimUsedInAOcode,nbasis1,nbasis2,&
             & nocc,nvirt,MinAuxBatch,numnodes,nAuxMPI,mynum,MaxNaux,nthreads)
        IF(noOMP)THEN
           noOMPsave= mylsitem%setting%scheme%noOMP
           mylsitem%setting%scheme%noOMP = .TRUE.
        ENDIF
     ENDIF
  ENDIF

  IF(DECinfo%RIMP2_deactivateOpenMP)THEN
     noOMP = .TRUE.
     noOMPsave= mylsitem%setting%scheme%noOMP
     mylsitem%setting%scheme%noOMP = .TRUE.
     nthreads = 1                           
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
  IF(DECinfo%PL.GT.0.OR.DECinfo%MemDebugPrint)THEN
     IF(dim1.EQ.nbasisAux)THEN
        print*,'NEW (alpha|AI) CODE: Perform Reduction'
     ELSE
        print*,'NEW (alpha|AI) CODE: Perform BCAST Loop using:',dim1
     ENDIF
     IF(MaxNaux.LT.dim1)THEN
        print*,'NEW (alpha|AI) CODE: Internal AotoMO, ',MaxNaux,' out of ',dim1
     ELSE
        print*,'NEW (alpha|AI) CODE: external AotoMO'
     ENDIF
     print*,'NEW (alpha|AI) CODE: Number of Threads =',nthreads
  ENDIF

  !Memory Estimate

  IF(DECinfo%MemDebugPrint)THEN
     IF(dim1.EQ.MaxNaux)THEN 
        !external AotoMO
        IF(use_bg_buf)THEN
           !Calpha(AuxDimUsedInAOcode*nvirt*nocc) have been allocated on BG buffer and is part of offset !
           MemoryEstimateGB = (AuxDimUsedInAOcode*nbasis1*nbasis2+AuxDimUsedInAOcode*nbasis1*nocc)*8.0E-9_realk 
           MemoryEstimateGB = MemoryEstimateGB + buf_realk%offset*8.0E-9_realk
           SaveMaxMemoryUsage = buf_realk%max_usage
           buf_realk%max_usage = buf_realk%offset
        ELSE 
           MemoryEstimateGB = (MAX(AuxDimUsedInAOcode*nbasis1*nbasis2+AuxDimUsedInAOcode*nbasis1*nocc,&
                & AuxDimUsedInAOcode*nbasis1*nocc+AuxDimUsedInAOcode*nvirt*nocc)+1)*8.0E-9_realk 
           MemoryEstimateGB = MemoryEstimateGB + mem_allocated_global*1.0E-9_realk
           SaveMaxMemoryUsage = max_mem_used_global
           max_mem_used_global = mem_allocated_global
        ENDIF
        WRITE(DECinfo%output,*)'RIUTILinfo: EXTERNAL MEM1=',AuxDimUsedInAOcode*nbasis1*nbasis2*8.0E-9_realk,'GB'
        WRITE(DECinfo%output,*)'RIUTILinfo: EXTERNAL MEM2=',AuxDimUsedInAOcode*nbasis1*nocc*8.0E-9_realk,'GB'
        WRITE(DECinfo%output,*)'RIUTILinfo: EXTERNAL MEM3=',AuxDimUsedInAOcode*nbasis1*nocc*8.0E-9_realk,'GB'
        IF(use_bg_buf)THEN
           WRITE(DECinfo%output,*)'RIUTILinfo: MEMORY IN USE=',buf_realk%offset*8.0E-9_realk,'GB'
        ELSE
           WRITE(DECinfo%output,*)'RIUTILinfo: EXTERNAL MEM4=',AuxDimUsedInAOcode*nvirt*nocc*8.0E-9_realk,'GB'
           WRITE(DECinfo%output,*)'RIUTILinfo: MEMORY IN USE=',mem_allocated_global*1.0E-9_realk,'GB'
        ENDIF
        WRITE(DECinfo%output,*)'RIUTILinfo: external AotoMO MemoryEstimateGB=',MemoryEstimateGB
        IF(use_bg_buf)THEN
           WRITE(DECinfo%output,*)'RIUTILinfo: buffer size                      ',mem_get_bg_buf_n()*8.0E-9_realk
!           WRITE(DECinfo%output,*)'RIUTILinfo: memory available                 ',mem_get_bg_buf_free()*8.0E-9_realk*0.90E0_realk
        ELSE
           WRITE(DECinfo%output,*)'RIUTILinfo: memory available                 ',maxsize
        ENDIF
     ELSE
        !internal AOtoMO
        IF(use_bg_buf)THEN
           MemoryEstimateGB=(dim1*nocc*nvirt+MaxNaux*nbasis1*nbasis2*nthreads+MaxNaux*nbasis1*nocc)*8.0E-9_realk &
                & + buf_realk%offset*8.0E-9_realk
           SaveMaxMemoryUsage = buf_realk%max_usage
           buf_realk%max_usage = buf_realk%offset
        ELSE
           MemoryEstimateGB=(dim1*nocc*nvirt+MaxNaux*nbasis1*nbasis2*nthreads+MaxNaux*nbasis1*nocc+&
                & nbasis1*nvirt+nbasis2*nocc)*8.0E-9_realk + mem_allocated_global*1.0E-9_realk 
           SaveMaxMemoryUsage = max_mem_used_global
           max_mem_used_global = mem_allocated_global
        ENDIF
        WRITE(DECinfo%output,*)'RIUTILinfo: INTERNAL MEM1=',dim1*nocc*nvirt*8.0E-9_realk,'GB',dim1*nocc*nvirt
        WRITE(DECinfo%output,*)'RIUTILinfo: INTERNAL MEM2=',MaxNaux*nbasis1*nbasis2*nthreads*8.0E-9_realk,'GB'
        WRITE(DECinfo%output,*)'RIUTILinfo: INTERNAL MEM3=',MaxNaux*nbasis1*nocc*8.0E-9_realk,'GB'
        IF(.NOT.use_bg_buf)THEN
           WRITE(DECinfo%output,*)'RIUTILinfo: INTERNAL MEM4=',nbasis1*nvirt*8.0E-9_realk,'GB'
           WRITE(DECinfo%output,*)'RIUTILinfo: INTERNAL MEM5=',nbasis2*nocc*8.0E-9_realk,'GB'
        ENDIF
        IF(use_bg_buf)THEN
           WRITE(DECinfo%output,*)'RIUTILinfo: MEMORY IN USE=',buf_realk%offset*8.0E-9_realk,'GB'
        ELSE
           WRITE(DECinfo%output,*)'RIUTILinfo: MEMORY IN USE=',mem_allocated_global*1.0E-9_realk,'GB'
        ENDIF
        WRITE(DECinfo%output,*)'RIUTILinfo: internal AotoMO MemoryEstimateGB=',MemoryEstimateGB
        IF(use_bg_buf)THEN
           WRITE(DECinfo%output,*)'RIUTILinfo: buffer size                      ',mem_get_bg_buf_n()*8.0E-9_realk
        ENDIF
     ENDIF
     
  ENDIF

  !Output: AlphaCD3(dim1,nvirt,nocc) 
  call II_get_RI_AlphaCD_3CenterIntFullOnAllNN(DECinfo%output,DECinfo%output,&
       & AlphaCD3,mylsitem%setting,nbasisAux,nbasis1,nbasis2,intspec(1:4),MaxNaux,&
       & nvirt,nocc,.TRUE.,Cvirt,Cocc,nthreads,dim1,GindexToLocal,DECinfo%PL,&
       & DECinfo%MemDebugPrint,use_bg_buf)
  
  IF(DECinfo%MemDebugPrint)THEN
     !-2 to avoid round off isssues. 
     IF(use_bg_buf)THEN
        MaxMemoryUsageGB = (buf_realk%max_usage-2)*8.0E-9_realk
        buf_realk%max_usage = MAX(SaveMaxMemoryUsage,buf_realk%max_usage)
     ELSE
        MaxMemoryUsageGB = (max_mem_used_global-2)*1.0E-9_realk
        max_mem_used_global = MAX(max_mem_used_global,SaveMaxMemoryUsage)
     ENDIF
     WRITE(DECinfo%output,*)'RIUTILinfo: MaxMemoryUsage=',MaxMemoryUsageGB
     WRITE(DECinfo%output,*)'RIUTILinfo: MemoryEstimateGB=',MemoryEstimateGB

     IF(use_bg_buf)THEN
        IF(MaxMemoryUsageGB.GT.MemoryEstimateGB)THEN
           print*,'RIUTILinfo: Error MaxMemoryUsage.GT.MemoryEstimateGB'
           print*,'RIUTILinfo: MaxMemoryUsage=',MaxMemoryUsageGB
           print*,'RIUTILinfo: MemoryEstimateGB=',MemoryEstimateGB
           call lsquit('RIUTILinfo: MaxMemoryUsage.GT.MemoryEstimateGB',-1)
        ENDIF
     ELSE
        !MemoryEstimateGB + 25% of available memory should be less that use memory        
        IF(MaxMemoryUsageGB.GT.MemInGBCollected)THEN
           print*,'RIUTILinfo: MaxMemoryUsage=',MaxMemoryUsageGB
           print*,'RIUTILinfo: MemoryEstimateGB=',MemoryEstimateGB
           print*,'RIUTILinfo: MemInGBCollected=',MemInGBCollected
           call lsquit('RIUTILinfo: MaxMemoryUsage.GT.MemoryEstimateGB',-1)
        ENDIF
     ENDIF
  ENDIF

  IF(noOMP)THEN
     mylsitem%setting%scheme%noOMP = noOMPsave  !restore default values
  ENDIF
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
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc TMPAlphaBetaDecomp(',size(TMPAlphaBetaDecomp),')'
           call mem_dealloc(TMPAlphaBetaDecomp)
        ENDIF
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: alloc Wprime(',nbasisAux*nbasisAux,')'
        call mem_alloc(Wprime,nbasisAux,nbasisAux,'RIUTIL:Wprime')
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
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc W(',nbasisAux*nbasisAux,')'
           call mem_alloc(W,nbasisAux,nbasisAux,'RIUTIL:W')
           call NAF_buildW(W,Wprime,AlphaBetaDecomp,nbasisAux)
           call NAF_SVD_W(W,Wprime,NBAR,nbasisAux,epsilon,NRED,SumSV,FullSumSV) 
           IF(DECinfo%PL.GT.2)THEN
              WRITE(DECinfo%output,*)'NAF: Auxiliary functions         = ',nbasisAux
              WRITE(DECinfo%output,*)'NAF: Natural Auxiliary functions = ',NRED
              WRITE(DECinfo%output,*)'NAF: Sum of included eigenvalues = ',SumSV
              WRITE(DECinfo%output,*)'NAF: Sum of neglected eigenvalues= ',FullSumSV-SumSV
              WRITE(DECinfo%output,*)'NAF: Sum of all eigenvalues      = ',FullSumSV
           ENDIF
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Wprime(',size(W),')'
           call mem_dealloc(W)
        ENDIF
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: alloc Wprime(',size(Wprime),')'
        call mem_dealloc(Wprime)
#ifdef VAR_MPI
        !BCAST OF NBAR
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        call ls_mpibcast(NRED,infpar%master,infpar%lg_comm)
        nbuf1 = NRED
        nbuf2 = nbasisAux
        IF(.NOT.master)THEN
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc NBAR(',NRED*nbasisAux,')'
           call mem_alloc(NBAR,NRED,nbasisAux,'RIUTIL:NBAR')
        ENDIF
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
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: alloc NBARTMP(',NREDLOC*nbasisAux,')'
        call mem_alloc(NBARTMP,NREDLOC,nbasisAux,'RIUTIL:NBARTMP')
        call buildNBARTMP(NBARTMP,NBAR,NREDLOC,NRED,nbasisAux,mynum,ndimMax2,numnodes)
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc NBAR(',size(NBAR),')'
        call mem_dealloc(NBAR)
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: alloc W(',NREDLOC*nbasisAux,')'
        call mem_alloc(W,NREDLOC,nbasisAux,'RIUTIL:W2') !used as TMP
        M =  NREDLOC         !rows of Output Matrix
        N =  nbasisAux       !columns of Output Matrix
        K =  nbasisAux       !summation dimension
        call dgemm('N','N',M,N,K,1.0E0_realk,NBARTMP,M,AlphaBetaDecomp,K,0.0E0_realk,W,M)
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc NBARTMP(',size(NBARTMP),')'
        call mem_dealloc(NBARTMP)
        M = NREDLOC          !rows of Output Matrix
        N = nvirt*nocc       !columns of Output Matrix
        K = nbasisAux        !summation dimension
        nsize = NREDLOC*nvirt*nocc
        IF(.NOT.use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha(',nsize,')'
           call mem_alloc(Calpha,nsize,'RIUTIL:Calpha')
        ENDIF
        call dgemm('N','N',M,N,K,1.0E0_realk,W,M,AlphaCD3,K,0.0E0_realk,Calpha,M)
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc alphaCD3(',size(alphaCD3),')'
           call mem_pseudo_dealloc(alphaCD3)
           !This looks weird but Calpha is currently pointing to the first 1:NBA*nvirt*nocc elements
           !of a "permanent" memory array. I only need the first 1:NREDLOC*nvirt*nocc elements
           !so I shrink the array dimension by deassociating (NOT deallocating) and reassociate
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha(',size(Calpha),')'
           call mem_pseudo_dealloc(Calpha)
           nsize = NREDLOC*nvirt*nocc
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha(',nsize,')'
           call mem_pseudo_alloc(Calpha,nsize)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc alphaCD3(',size(alphaCD3),')'
           call mem_dealloc(alphaCD3)
        ENDIF
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc W(',size(W),')'
        call mem_dealloc(W)
        NBA = NREDLOC
     ELSE
        IF(CollaborateWithSlaves)then 
           nsize = NBA*nvirt*(nocc*i8)
           IF(.NOT.use_bg_buf)THEN 
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha(',nsize,')'
              call mem_alloc(Calpha,nsize,'RIUTIL:Calpha2')
           ENDIF
           !Calpha = TMPAlphaBetaDecomp(MynbasisAuxMPI,nbasisAux)
           M =  NBA              !rows of Output Matrix
           N =  nvirt*nocc       !columns of Output Matrix
           K =  nbasisAux        !summation dimension
           call dgemm('N','N',M,N,K,1.0E0_realk,TMPAlphaBetaDecomp,M,alphaCD3,K,0.0E0_realk,Calpha,M)
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc TMPAlphaBetaDecomp(',size(TMPAlphaBetaDecomp),')'
           call mem_dealloc(TMPAlphaBetaDecomp)
        ELSE !Serial version        
           nsize = nbasisAux*nvirt*nocc
           IF(.NOT.use_bg_buf)THEN
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha(',nsize,')'
              call mem_alloc(Calpha,nsize,'RIUTIL:Calpha3')
           ENDIF
           M =  nbasisAux        !rows of Output Matrix
           N =  nvirt*nocc       !columns of Output Matrix
           K =  nbasisAux        !summation dimension
           call dgemm('N','N',M,N,K,1.0E0_realk,AlphaBetaDecomp,M,AlphaCD3,K,0.0E0_realk,Calpha,M)
        ENDIF
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc alphaCD3(',size(alphaCD3),')'
           call mem_pseudo_dealloc(alphaCD3)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc alphaCD3(',size(alphaCD3),')'
           call mem_dealloc(alphaCD3)
        ENDIF
     ENDIF
     !  print*,'MY RIMP2 INTEGRAL AlphaCD2(1:nA,1:4) NEW VERSION MYNUM',MYNUM
     !  call ls_output(AlphaCD3,1,size(AlphaCD3,1),1,4,size(AlphaCD3,1),nvirt*nocc,1,6)
  ELSE
     !=====================================================================================
     ! MPI scheme:  Bcast Routine
     !=====================================================================================
!     print*,'NBA',NBA     
     IF(.NOT.use_bg_buf)THEN
        nsize = NBA*nvirt*(nocc*i8)
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha(',nsize,')'
        call mem_alloc(Calpha,nsize,'RIUTIL:Calpha4')
     ENDIF
     call ls_dzero8(Calpha,nsize)
     IF(DECinfo%NAF)THEN
        nsize = nbasisAux*nbasisAux
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: alloc Wprime(',nsize,')'
        call mem_alloc(Wprime,nbasisAux,nbasisAux,'RIUTIL:Wprime2')
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
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc alphaCD5(',size(alphaCD5),')'
              call mem_pseudo_alloc(AlphaCD5,nsize)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: alloc AlphaCD5(',nsize,')'
              call mem_alloc(AlphaCD5,nsize,'RIUTIL:AlphaCD5')
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
              IF(DECinfo%MemDebugPrint)call printBGinfo()
              IF(DECinfo%MemDebugPrint)print*,'BG: dealloc alphaCD5(',size(alphaCD5),')'
              call mem_pseudo_dealloc(AlphaCD5)
           ELSE
              IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
              IF(DECinfo%MemDebugPrint)print*,'STD: dealloc alphaCD5(',size(alphaCD5),')'
              call mem_dealloc(AlphaCD5)
           ENDIF
        ENDIF
     ENDDO
     IF(use_bg_buf)THEN
        IF(DECinfo%MemDebugPrint)call printBGinfo()
        IF(DECinfo%MemDebugPrint)print*,'BG: dealloc alphaCD3(',size(alphaCD3),')'
        call mem_pseudo_dealloc(AlphaCD3)
     ELSE
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc alphaCD3(',size(alphaCD3),')'
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
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc W(',nbasisAux*nbasisAux,')'
           call mem_alloc(W,nbasisAux,nbasisAux,'RIUTIL:WC')
           call NAF_buildW(W,Wprime,AlphaBetaDecomp,nbasisAux)           
           call NAF_SVD_W(W,Wprime,NBAR,nbasisAux,epsilon,NRED,SumSV,FullSumSV) 
           IF(DECinfo%PL.GT.2)THEN
              WRITE(DECinfo%output,*)'NAF: Auxiliary functions         = ',nbasisAux
              WRITE(DECinfo%output,*)'NAF: Natural Auxiliary functions = ',NRED
              WRITE(DECinfo%output,*)'NAF: Sum of included eigenvalues = ',SumSV
              WRITE(DECinfo%output,*)'NAF: Sum of neglected eigenvalues= ',FullSumSV-SumSV
              WRITE(DECinfo%output,*)'NAF: Sum of all eigenvalues      = ',FullSumSV
           ENDIF
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc W(',size(W),')'
           call mem_dealloc(W)
        ENDIF
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Wprime(',size(Wprime),')'
        call mem_dealloc(Wprime)
#ifdef VAR_MPI
        call time_start_phase( PHASE_IDLE )
        call lsmpi_barrier(infpar%lg_comm)
        call time_start_phase( PHASE_COMM )
        call ls_mpibcast(NRED,infpar%master,infpar%lg_comm)
        nbuf1 = NRED
        nbuf2 = nbasisAux
        IF(.NOT.master)THEN
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc NBAR(',NRED*nbasisAux,')'
           call mem_alloc(NBAR,NRED,nbasisAux,'RIUTIL:NBARC')
        ENDIF
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
        IF(.NOT.use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: alloc CalphaNAF(',nsize,')'
           call mem_alloc(CalphaNAF,nsize,'RIUTIL:CalphaNAF')
        ENDIF
        call ls_dzero8(CalphaNAF,nsize)
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: alloc NBARTMP(',NREDLOC*nbasisAux,')'
        call mem_alloc(NBARTMP,NREDLOC,nbasisAux,'RIUTIL:NBARTMPD')
        call buildNBARTMP(NBARTMP,NBAR,NREDLOC,NRED,nbasisAux,mynum,ndimMax2,numnodes)
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc NBAR(',size(NBAR),')'
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
                 IF(DECinfo%MemDebugPrint)call printBGinfo()
                 IF(DECinfo%MemDebugPrint)print*,'BG: alloc Calpha2(',size(Calpha2),')'
                 call mem_pseudo_alloc(Calpha2,nsize)
              ELSE
                 IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
                 IF(DECinfo%MemDebugPrint)print*,'STD: alloc Calpha2(',nsize,')'
                 call mem_alloc(Calpha2,nsize,'RIUTIL:Calpha2D')
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
                 IF(DECinfo%MemDebugPrint)call printBGinfo()
                 IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha2(',size(Calpha2),')'
                 call mem_pseudo_dealloc(Calpha2)
              ELSE
                 IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
                 IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha2(',size(Calpha2),')'
                 call mem_dealloc(Calpha2)
              ENDIF
           ENDIF
        ENDDO
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc NBARTMP(',size(NBARTMP),')'
        call mem_dealloc(NBARTMP)
        NBA = NREDLOC
        IF(use_bg_buf)THEN
           IF(DECinfo%MemDebugPrint)call printBGinfo()
           IF(DECinfo%MemDebugPrint)print*,'BG: dealloc Calpha(',size(Calpha),')'
           call mem_pseudo_dealloc(Calpha)
           Calpha => CalphaNAF
           call lsquit('clearly not working',-1)
        ELSE
           IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
           IF(DECinfo%MemDebugPrint)print*,'STD: dealloc Calpha(',size(Calpha),')'
           call mem_dealloc(Calpha)
           Calpha => CalphaNAF           
        ENDIF
        CALL LSTIMER('DF_Calpha:NAF',TS3,TE3,LUPRI,ForcePrint)
     ENDIF
     !allocated inside II_get_RI_AlphaCD_3CenterIntFullOnAllNNdim
     IF(CollaborateWithSlaves.OR.DECinfo%RIMP2ForcePDMCalpha)THEN
        IF(DECinfo%MemDebugPrint)call stats_globalmem(6)
        IF(DECinfo%MemDebugPrint)print*,'STD: dealloc TMPAlphaBetaDecomp(',size(TMPAlphaBetaDecomp),')'
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

subroutine BuildnAuxMPIUsedRI(nbasisAux,numnodes,nbasisAuxMPI2)
  implicit none
  integer,intent(in) :: nbasisAux,numnodes
  integer,intent(inout) :: nbasisAuxMPI2(numnodes)
  !local variables
  integer :: ndimMax1,I,J
  ndimMax1 = nbasisAux/numnodes
  do I=1,numnodes
     nbasisAuxMPI2(I) = ndimMax1
  enddo
  J=2 !not add to master
  do I=1,MOD(nbasisAux,numnodes)
     nbasisAuxMPI2(J) = nbasisAuxMPI2(J) + 1
     J=J+1
  enddo
end subroutine BuildnAuxMPIUsedRI

subroutine BuildnAuxMPIUsedRIinfo(nbasisAux,numnodes,mynum,AuxMPIstartMy,iAuxMPIextraMy)
  implicit none
  integer,intent(in) :: nbasisAux,numnodes,mynum
  integer,intent(inout) :: AuxMPIstartMy
  integer,intent(inout) :: iAuxMPIextraMy
  !local variables 
  integer :: ndimMax1
  ndimMax1 = nbasisAux/numnodes
  AuxMPIstartMy = mynum*ndimMax1
  iAuxMPIextraMy = 0
  IF(MOD(nbasisAux,numnodes).GT.0)THEN
     IF(mynum.GT.0)THEN
        iAuxMPIextraMy = numnodes*ndimMax1 + mynum
     ENDIF
     IF(iAuxMPIextraMy.GT.nbasisAux)iAuxMPIextraMy=0
  ENDIF
end subroutine BuildnAuxMPIUsedRIinfo

subroutine DetermineMaxNauxRI(use_bg_buf,noOMP,dim1,AuxDimUsedInAOcode,nbasis1,nbasis2,&
     & nocc,nvirt,MinAuxBatch,numnodes,nAuxMPI,mynum,MaxNaux,nthreads)
  implicit none       
  logical,intent(in) :: use_bg_buf
  logical,intent(inout) :: noOMP
  integer,intent(in)    :: dim1,AuxDimUsedInAOcode,nbasis1,nbasis2,nocc,nvirt,MinAuxBatch
  integer,intent(in)    :: numnodes,mynum
  integer,intent(in)    :: nAuxMPI(numnodes)
  integer,intent(inout) :: MaxNaux,nthreads
  !
  real(realk) :: maxsize,MemInGBCollected,ExternalAOtoMO,InternalAOtoMO 
  !Determine memory available 
  IF(use_bg_buf)THEN
     maxsize = mem_get_bg_buf_free()*8.0E-9_realk*0.90E0_realk
     IF(DECinfo%MemDebugPrint)THEN
        print*,'DetermineMaxNauxRI: mem_get_bg_buf_free=',mem_get_bg_buf_free()
        print*,'DetermineMaxNauxRI: maxsize',maxsize,' GB'
        call printBGinfo()
     ENDIF
  ELSE
     call get_currently_available_memory(MemInGBCollected)
     maxsize = MemInGBCollected*0.75E0_realk
     IF(DECinfo%MemDebugPrint)print*,'BCAST Scheme: maxsize',maxsize,' GB'
  ENDIF

  !Memory estimate for the external AO to MO transformation
  IF(use_bg_buf)THEN
     !FullAlphaCD(dim1*nbasis1*nbasis2) + TmpAlphaCD(dim1*nbasis1*nocc) - No AlphaCD(dim1*nvirt*nocc) 
     !because FullAlphaCD(dim1*nbasis1*nbasis2) is used and then the size is shrunk
     ExternalAOtoMO = (AuxDimUsedInAOcode*nbasis1*nbasis2+AuxDimUsedInAOcode*nbasis1*nocc)*8.0E-9_realk 
  ELSE 
     ExternalAOtoMO = (MAX(AuxDimUsedInAOcode*nbasis1*nbasis2+AuxDimUsedInAOcode*nbasis1*nocc,&
          & AuxDimUsedInAOcode*nbasis1*nocc+AuxDimUsedInAOcode*nvirt*nocc)+1)*8.0E-9_realk 
  ENDIF
  
  IF(DECinfo%MemDebugPrint)THEN
     WRITE(DECinfo%output,*)'RIUTILinfo: External AO to MO requires',ExternalAOtoMO,' GB' 
     IF(use_bg_buf)THEN
        WRITE(DECinfo%output,*)'RIUTILinfo: external AotoMO MemoryEstimateGB=',&
             & (buf_realk%offset+ExternalAOtoMO/8.0E-9_realk)*8.0E-9_realk,' GB'
        WRITE(DECinfo%output,*)'RIUTILinfo: buffer size                      ',mem_get_bg_buf_n()*8.0E-9_realk
        WRITE(DECinfo%output,*)'RIUTILinfo: 90% buffer size                  ',mem_get_bg_buf_n()*8.0E-9_realk*0.9E0_realk
     ENDIF
     WRITE(DECinfo%output,*)'RIUTILinfo: maxsize                   ',maxsize,' GB' 
  ENDIF

  IF(ExternalAotoMO.LE.maxsize)THEN 
     IF(DECinfo%PL.GT.0.OR.DECinfo%MemDebugPrint)THEN
        print*,'External AO to MO Chosen'
     ENDIF
     !External AO to MO chosen
     MaxNaux = dim1
  ELSE
     !Internal AO to MO chosen. We determine the parameter MaxNaux
     !which denotes the block (MaxNaux,nbasis1,nbasis2) that should be calculated 
     !before it is transformed into (MaxNaux,nMO1,nMO2) and pluged into AlphaCD(dim1,nvirt,nocc)
     
     !Assume minimal MaxNaux (MaxNaux=MinAuxBatch) - If this is still too much memory 
     !We need to deactivate OpenMP!
     !   FullAlphaCD+                        W0
     IF(use_bg_buf)THEN
        InternalAOtoMO = (dim1*nocc*nvirt+MinAuxBatch*nbasis1*nbasis2*nthreads+MinAuxBatch*nbasis1*nocc)*8.0E-9_realk
     ELSE
        InternalAOtoMO = (dim1*nocc*nvirt+MinAuxBatch*nbasis1*nbasis2*nthreads+MinAuxBatch*nbasis1*nocc+&
             & nbasis1*nvirt+nbasis2*nocc)*8.0E-9_realk 
     ENDIF
     IF(DECinfo%PL.GT.0.OR.DECinfo%MemDebugPrint)THEN
        print*,'Minimal MaxNaux = ',MinAuxBatch
        print*,'InternalAOtoMO',InternalAOtoMO
        IF(use_bg_buf)THEN
           print*,'A Resulting in Memory Estimate of  ',&
                & (buf_realk%offset+InternalAOtoMO/8.0E-9_realk)*8.0E-9_realk,' GB'           
           print*,'90% buffer size of                 ',buf_realk%nmax*8.0E-9_realk*0.9E0_realk,' GB'
        ENDIF
        print*,'maxsize                 ',maxsize
     ENDIF
     IF(InternalAOtoMO.GT.maxsize)THEN        
        IF(nthreads.GT.1)THEN
           IF(DECinfo%MemDebugPrint)print*,'Deactivating OpenMP due to memory requirements'
           !use Critical instead ! 
           !Deactivating OpenMP 
           noOMP = .TRUE.
           nthreads = 1                           
        ELSE
           CALL lsquit('Not enough memory in build_calpha X1',-1)
        ENDIF
        IF(use_bg_buf)THEN
           InternalAOtoMO = (dim1*nocc*nvirt+MinAuxBatch*nbasis1*nbasis2*nthreads+MinAuxBatch*nbasis1*nocc)*8.0E-9_realk
        ELSE
           InternalAOtoMO = (dim1*nocc*nvirt+MinAuxBatch*nbasis1*nbasis2*nthreads+MinAuxBatch*nbasis1*nocc+&
                & nbasis1*nvirt+nbasis2*nocc)*8.0E-9_realk 
        ENDIF
        print*,'Minimal MaxNaux = ',MinAuxBatch,' and No OpenMP '
        print*,'InternalAOtoMO',InternalAOtoMO
        IF(use_bg_buf)THEN
           print*,'B Resulting in Memory Estimate of  ',&
                & (buf_realk%offset+InternalAOtoMO/8.0E-9_realk)*8.0E-9_realk,' GB'           
           print*,'90% buffer size of                 ',buf_realk%nmax*8.0E-9_realk*0.9E0_realk,' GB'
        ENDIF
        print*,'maxsize                 ',maxsize
        IF(InternalAOtoMO.GT.maxsize)THEN
           CALL lsquit('Not enough memory in build_calpha X2, even with OpenMP deactivated',-1)
        ELSE           
           !find Optimal MaxNaux (with nthreads=1)
           IF(use_bg_buf)THEN
              MaxNaux = FLOOR((MaxSize/8.0E-9_realk-dim1*nvirt*nocc) &
                   & /(nbasis1*nocc+nbasis1*nbasis2*nthreads))
              IF(MaxNaux.GT.dim1)MaxNaux = dim1 
              IF(MaxNaux.LT.1)call lsquit('Weird error should never happen',-1)
              InternalAOtoMO = (dim1*nocc*nvirt+MaxNaux*nbasis1*nbasis2*nthreads+MaxNaux*nbasis1*nocc)*8.0E-9_realk

              IF(DECinfo%PL.GT.0.OR.DECinfo%MemDebugPrint)THEN
                 print*,'MaxNaux',MaxNaux
                 print*,'InternalAOtoMO',InternalAOtoMO
                 print*,'C Resulting in Memory Estimate of  ',&
                      & (buf_realk%offset+InternalAOtoMO/8.0E-9_realk)*8.0E-9_realk,' GB'           
                 print*,'90% buffer size of                 ',buf_realk%nmax*8.0E-9_realk*0.9E0_realk,' GB'
                 print*,'maxsize                 ',maxsize
              ENDIF
           ELSE
              MaxNaux = FLOOR((MaxSize/8.0E-9_realk-dim1*nvirt*nocc-nbasis1*nvirt-nbasis2*nocc) &
                   & /(nbasis1*nocc+nbasis1*nbasis2*nthreads))
              IF(MaxNaux.GT.dim1)MaxNaux = dim1 
              IF(MaxNaux.LT.1)call lsquit('Weird error should never happen',-1)
              InternalAOtoMO = (dim1*nocc*nvirt+MaxNaux*nbasis1*nbasis2*nthreads+MaxNaux*nbasis1*nocc &
                   & +nbasis1*nvirt+nbasis2*nocc)*8.0E-9_realk               
              IF(DECinfo%PL.GT.0.OR.DECinfo%MemDebugPrint)THEN
                 print*,'InternalAOtoMO',InternalAOtoMO,' GB'
                 print*,'maxsize       ',maxsize,' GB'
              ENDIF
           ENDIF           
           IF(InternalAOtoMO.GT.maxsize)THEN
              CALL lsquit('Miscalc of MaxNaux',-1)
           ENDIF
        ENDIF
     ELSE
        IF(DECinfo%PL.GT.0.OR.DECinfo%MemDebugPrint)THEN
           print*,'minimal MaxNaux=',MinAuxBatch,' OK -> Find optimal MaxNaux:'
        ENDIF
        !minimal MaxNaux OK -> Find optimal MaxNaux
        !MinAuxBatch = (maxsize/8.0E-9_realk-dim1*nocc*nvirt)/(nbasis1*nbasis2*nthreads+nbasis1*nocc)
        IF(use_bg_buf)THEN
           MaxNaux = MIN(dim1,&
                & FLOOR((MaxSize/8.0E-9_realk-dim1*nvirt*nocc) &
                & /(nbasis1*nocc+nbasis1*nbasis2*nthreads)))

           IF(DECinfo%PL.GT.0.OR.DECinfo%MemDebugPrint)THEN
              print*,'MaxNaux',MaxNaux
              InternalAOtoMO = (dim1*nocc*nvirt+MaxNaux*nbasis1*nbasis2*nthreads+MaxNaux*nbasis1*nocc)*8.0E-9_realk
              print*,'InternalAOtoMO',InternalAOtoMO,' GB'
              print*,'DetermineMaxNauxRI(elements)A:',dim1*nocc*nvirt
              print*,'DetermineMaxNauxRI(elements)B:',MaxNaux*nbasis1*nbasis2*nthreads
              print*,'DetermineMaxNauxRI(elements)C:',MaxNaux*nbasis1*nocc              
              print*,'D Resulting in Memory Estimate of  ',&
                   & (buf_realk%offset+InternalAOtoMO/8.0E-9_realk)*8.0E-9_realk,' GB'           
              print*,'90% buffer size of                 ',buf_realk%nmax*8.0E-9_realk*0.9E0_realk,' GB'
              print*,'maxsize                 ',maxsize,' GB'
           ENDIF
        ELSE
           MaxNaux = MIN(dim1,&
                &FLOOR((MaxSize/8.0E-9_realk-dim1*nvirt*nocc-nbasis1*nvirt-nbasis2*nocc) &
                & /(nbasis1*nocc+nbasis1*nbasis2*nthreads)))
           IF(DECinfo%PL.GT.0.OR.DECinfo%MemDebugPrint)THEN
              print*,'MaxNaux',MaxNaux
              InternalAOtoMO = (dim1*nocc*nvirt+MaxNaux*nbasis1*nbasis2*nthreads+MaxNaux*nbasis1*nocc &
                   & +nbasis1*nvirt+nbasis2*nocc)*8.0E-9_realk
              print*,'InternalAOtoMO',InternalAOtoMO
              print*,'InternalAOtoMO(Elements)',InternalAOtoMO/8.0E-9_realk,' of the allowed',maxsize,' GB'
           ENDIF
        ENDIF
     ENDIF
  ENDIF
END subroutine DetermineMaxNauxRI

!This should be call my master and slaves
subroutine Build_RobustERImatU(myLSitem,master,nbasisAux,&
     & LUPRI,FORCEPRINT,&
     & CollaborateWithSlaves,mynum,numnodes,&
     & AlphaBetaDecomp,intspec,Umat)
  implicit none
  type(lsitem), intent(inout) :: mylsitem
  integer,intent(in) :: nbasisAux,LUPRI,mynum,numnodes
  logical,intent(in) :: master,FORCEPRINT,CollaborateWithSlaves
  real(realk),intent(in) :: AlphaBetaDecomp(nbasisAux,nbasisAux)
  real(realk),intent(inout) :: Umat(nbasisAux,nbasisAux)
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
!  print*,'Build_RobustERImatU intspec=',intspec
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
     call mem_alloc(TMP,N,N,'RobustERImatU:TMP')
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
  call LSDPOTRF('U', N, A, N, INFO ) !U=Cholesky factor
  IF(INFO.ne. 0) THEN
     print *, 'DPOTRF NR 1 Failed in Get_InverseCholeskyFactor',INFO
     call lsquit('DPOTRF NR 1 Failed in Get_InverseCholeskyFactor',-1)
  ENDIF
  call LSDPOTRI('U', N, A, N, INFO ) !U=inverse of a original U
  IF(INFO.ne. 0) THEN
     print *, 'DPOTRI Failed in Get_InverseCholeskyFactor',INFO
     call lsquit('DPOTRI Failed in Get_InverseCholeskyFactor',-1)
  ENDIF
  call LSDPOTRF('U', N, A, N, INFO ) !U=Cholesky factor of inverse of a original U
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
offset2 = numnodes*ndimMax1 + mynum 
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
  call mem_alloc(NBAR,nred,N,'NAFSVDW:NBAR')
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
  call mem_alloc(TMP,N,N,'NAF_buildW:TMP')
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
  !$OMP PARALLEL DO DEFAULT(NONE) &
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
     & Calpha,TMPAlphaBeta_inv,nbasisAux,MynbasisAuxMPI)
  implicit none
  integer,intent(in) :: nocc,nvirt,myOriginalRank,numnodes,natoms,MynbasisAuxMPI
  integer,intent(in) :: OriginalRanknbasisAuxMPI,nbasisAux
  integer,intent(in) :: nAtomsMPI(numnodes)
  integer,intent(in) :: startAuxMPI(nAtoms,numnodes)
  integer,intent(in) :: nAuxMPI(nAtoms,numnodes)
  real(realk),intent(in) :: AlphaCD5(OriginalRanknbasisAuxMPI,nvirt*nocc)
  real(realk),intent(inout) :: Calpha(MynbasisAuxMPI,nvirt*nocc)
  real(realk),intent(in) :: TMPAlphaBeta_inv(MynbasisAuxMPI,nbasisAux)
  !
  integer :: IB,startB2,iAtomB,StartB,BETA,ALPHA
  real(realk) :: TMP
  !Step 2: Obtain part of Calpha from this contribution
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(IB,startB2,iAtomB,startB,TMP,ALPHA,BETA) &
!$OMP SHARED(nocc,nvirt,myOriginalRank,numnodes,natoms,MynbasisAuxMPI,&
!$OMP        OriginalRanknbasisAuxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,&
!$OMP        AlphaCD5,Calpha,TMPAlphaBeta_inv,nbasisAux)
  do IB = 1,nocc*nvirt
     startB2 = 0
     DO iAtomB=1,nAtomsMPI(myOriginalRank+1)
        StartB = startAuxMPI(iAtomB,myOriginalRank+1)
        do BETA = 1,nAuxMPI(iAtomB,myOriginalRank+1)
           TMP = AlphaCD5(startB2 + BETA,IB)
           do ALPHA = 1,MynbasisAuxMPI
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
     & TMPAlphaBeta_inv,nbasisAux)
  implicit none
  integer,intent(in) :: nocc,nvirt,mynum,numnodes,natoms
  integer,intent(in) :: MynbasisAuxMPI,nbasisAux
  integer,intent(in) :: nAtomsMPI(numnodes)
  integer,intent(in) :: startAuxMPI(nAtoms,numnodes)
  integer,intent(in) :: nAuxMPI(nAtoms,numnodes)
  real(realk),intent(in) :: AlphaCD3(MynbasisAuxMPI,nvirt*nocc)
  real(realk),intent(inout) :: Calpha(MynbasisAuxMPI,nvirt*nocc)
  real(realk),intent(in) :: TMPAlphaBeta_inv(MynbasisAuxMPI,nbasisAux)
  !
  integer :: IB,startB2,iAtomB,StartB,BETA,ALPHA
  real(realk) :: TMP

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(IB,startB2,iAtomB,StartB,BETA,ALPHA,TMP) &
!$OMP SHARED(nocc,nvirt,nAtomsMPI,startAuxMPI,mynum,&
!$OMP        nAuxMPI,AlphaCD3,Calpha,TMPAlphaBeta_inv,MynbasisAuxMPI) 
  do IB = 1,nocc*nvirt
     startB2 = 0
     DO iAtomB=1,nAtomsMPI(mynum+1)
        StartB = startAuxMPI(iAtomB,mynum+1)
        do BETA = 1,nAuxMPI(iAtomB,mynum+1)
           TMP = AlphaCD3(startB2 + BETA,IB)
           do ALPHA = 1,MynbasisAuxMPI
              Calpha(ALPHA,IB) = Calpha(ALPHA,IB) + &
                   & TMPAlphaBeta_inv(ALPHA,startB + BETA)*TMP
           enddo
        enddo
        startB2 = startB2 + nAuxMPI(iAtomB,mynum+1)
     ENDDO
  enddo
!$OMP END PARALLEL DO
end subroutine RIMP2_buildOwnCalphaFromAlphaCD

subroutine Build_CalphaMO(myLSitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
     & CollaborateWithSlaves,Cvirt,nvirt,Cocc,nocc,mynum,numnodes,nAtomsAux,Calpha,&
     & NBA,AlphaBetaDecomp,AlphaBetaDecompCreate,SymDecomp,Oper)
  implicit none
  type(lsitem), intent(inout) :: mylsitem
  integer,intent(inout) :: NBA
  integer,intent(in) :: nAtomsAux,nocc,nvirt
  integer,intent(in) :: nbasisAux,LUPRI,nbasis,mynum,numnodes
  logical,intent(in) :: master,FORCEPRINT,CollaborateWithSlaves,AlphaBetaDecompCreate
  logical,intent(in) :: SymDecomp
  ! If SymDecomp=True apply the full apply the decomposition of (alpha|beta)^{-1} to (alpha|ai) 
  ! in a symmetric fashion:  
  ! C(alpha,a,i) = (alpha|beta)^{-1/2} (beta|ai) (or the corresponding Cholesky)
  ! so that the Integral can be obtained from g(a,i,b,j) = C(alpha,a,i)*C(alpha,b,j)
  ! If SymDecomp=False apply the full decomposition to the right: 
  ! C(alpha,a,i) = (alpha|beta)^{-1}(beta|ai)
  real(realk) :: AlphaBetaDecomp(nbasisAux,nbasisAux)
  real(realk),intent(in) :: Cocc(nbasis,nocc),Cvirt(nbasis,nvirt)
  real(realk),pointer :: Calpha(:)
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
  integer :: ndimMax,rimp2_nodtot,nP
  integer :: nbuf1,nbuf2,nbuf3,inode,J
  logical :: useAlphaCD5,useAlphaCD6,ChangedDefault,first_order,MessageRecieved
  logical :: PerformReduction,RIMPSubGroupCreated,UseSubGroupCommunicator
  integer :: startB2,StartB,iatomB
  PerformReduction = .TRUE.
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
           WRITE(DECinfo%output,'(A,F8.1,A)')'RIUTIL: Full (alpha|cd) integral requires ',SizeCalpha,' GB'
           WRITE(DECinfo%output,'(A,F8.1,A)')'RIUTIL: Memory available                  ',MemInGBCollected,' GB'
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
     call mem_alloc(nbasisAuxMPI,numnodes,'RIUTIL1:nbasisAuxMPI')  !number of Aux basis func assigned to rank
     nbasisAuxMPI = 0 
     call mem_alloc(nAtomsMPI,numnodes,'RIUTIL1:nAtomsMPI')        !atoms assign to rank
     call mem_alloc(startAuxMPI,nAtomsAux,numnodes,'RIUTIL1:startAuxMPI')  !startindex in full (nbasisAux)
     call mem_alloc(AtomsMPI,nAtomsAux,numnodes,'RIUTIL1:AtomsMPI')!identity of atoms in full molecule
     call mem_alloc(nAuxMPI,nAtomsAux,numnodes,'RIUTIL1:nAuxMPI')  !nauxBasis functions for each of the nAtomsMPI
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
        call mem_alloc(AlphaBeta,nbasisAux,nbasisAux,'RIUTIL1:AlphaBeta')
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
        IF(SymDecomp)THEN
           ! Create the inverse square root AlphaBeta = (alpha|beta)^(-1/2)
           ! Warning the inverse is not unique so in order to make sure all slaves have the same
           ! inverse matrix we calculate it on the master a BCAST to slaves
           call lowdin_diag_S_minus_sqrt(nbasisAux, AlphaBeta,AlphaBetaDecomp, lupri)
        ELSE
           call lowdin_diag_S_minus1(nbasisAux, AlphaBeta,AlphaBetaDecomp, lupri)
        ENDIF
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

  IF(CollaborateWithSlaves)then 
     call mem_alloc(TMPAlphaBetaDecomp,MynbasisAuxMPI,nbasisAux,'RIUTIL1:TMPAlphaBetaDecomp')
     StartB2 = 0
     DO iAtomB=1,nAtomsMPI(mynum+1)
        StartB = startAuxMPI(iAtomB,mynum+1)
        nP = nAuxMPI(iAtomB,mynum+1)
        do I=1,nbasisAux           
           DO J=1,nP
              TMPAlphaBetaDecomp(StartB2+J,I) = AlphaBetaDecomp(StartB+J,I)
           enddo
        enddo
        StartB2 = StartB2 + nP
     ENDDO
     NBA = MynbasisAuxMPI
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
        call mem_alloc(alphaCDFull,nbasisAux,nvirt,nocc,'RIUTIL1:alphaCDFull')
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
        M =  MynbasisAuxMPI   !rows of Output Matrix
        N =  nvirt*nocc       !columns of Output Matrix
        K =  nbasisAux        !summation dimension
        call mem_alloc(Calpha,MynbasisAuxMPI*(i8*nvirt)*nocc,'RIUTIL1:Calpha3')
        call dgemm('N','N',M,N,K,1.0E0_realk,TMPAlphaBetaDecomp,&
             & M,alphaCDFull,K,0.0E0_realk,Calpha,M)
        call mem_dealloc(alphaCDFull)
        call mem_dealloc(TMPAlphaBetaDecomp)
     ELSE
        !Serial version
        M =  MynbasisAuxMPI   !rows of Output Matrix
        N =  nvirt*nocc       !columns of Output Matrix
        K =  nbasisAux        !summation dimension
        call mem_alloc(Calpha,MynbasisAuxMPI*(i8*nvirt)*nocc,'RIUTIL1:Calpha4')
        call dgemm('N','N',M,N,K,1.0E0_realk,AlphaBetaDecomp,&
             & M,AlphaCD3,K,0.0E0_realk,Calpha,M)
        call mem_dealloc(AlphaCD3)
     ENDIF
  ELSE

     !=====================================================================================
     ! MPI scheme:  Bcast Routine
     !=====================================================================================
     call mem_alloc(Calpha,MynbasisAuxMPI*(i8*nvirt)*nocc,'RIUTIL1:Calpha6')
     nsize = MynbasisAuxMPI*nvirt*nocc
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
                & AlphaCD3,Calpha,TMPAlphaBetaDecomp,nbasisAux)
           call mem_dealloc(AlphaCD3)
        ELSE
           nbuf1 = nbasisAuxMPI(inode)
           nbuf2 = nvirt
           nbuf3 = nocc
           node = inode-1
           !recieve
           call mem_alloc(AlphaCD5,nbasisAuxMPI(inode),nvirt,nocc,'RIUTIL1:AlphaCD5')
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
                & Calpha,TMPAlphaBetaDecomp,nbasisAux,MynbasisAuxMPI)
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
     call mem_dealloc(AtomsMPI) !not used in this subroutine 
  ENDIF
!  call sleep(mynum*5)
!  PRINT*,'MynbasisAuxMPI2',MynbasisAuxMPI2
!  WRITE(6,*)'Final Calph(NBA=',NBA,',nvirt=',nvirt,',nocc=',nocc,')'
!  WRITE(6,*)'Print Subset Final Calph(NBA=',NBA,',1:4)  MYNUM',MYNUM
!  call ls_output(Calpha,1,NBA,1,4,NBA,nvirt*nocc,1,6)

end subroutine Build_CalphaMO

!This subroutine calculates the DEC-RI-MP2 gradient contribution:
!the derivatie of the 3 term energy contribution
!sum_{a,i,b,j} C(alpha,a,i) (alpha|bj) Thetha(a,i,b,j) +
!sum_{a,i,b,j} (ai|alpha) C(alpha,b,j) Thetha(a,i,b,j) -
!sum_{a,i,b,j} C(alpha,a,i) (alpha|beta) C(beta,b,j) Thetha(a,i,b,j)
! =
! the 6 term gradient contribution 
!sum_{a,i,b,j} C(alpha,a,i) (alpha^x|bj) Thetha(a,i,b,j) +
!sum_{a,i,b,j} C(alpha,a,i) (alpha|(bj)^x) Thetha(a,i,b,j) + 
!sum_{a,i,b,j} ((ai)^x|alpha) C(alpha,b,j) Thetha(a,i,b,j) + 
!sum_{a,i,b,j} (ai|alpha^x) C(alpha,b,j) Thetha(a,i,b,j) - 
!sum_{a,i,b,j} C(alpha,a,i) (alpha^x|beta) C(beta,b,j) Thetha(a,i,b,j) -
!sum_{a,i,b,j} C(alpha,a,i) (alpha|beta^x) C(beta,b,j) Thetha(a,i,b,j)
!due to symmetry this can be written as A,B,C (the last 2 terms are collected into C)
! A. C(alpha,a,i) (alpha^x|bj) Thetha(a,i,b,j) +
! B. C(alpha,a,i) (alpha|(bj)^x) Thetha(a,i,b,j) + 
! C. - (P|Q)^x*Cpq(P,Q) 
! However we do the evaluation in AO basis so it becomes
! A. (P^x|bj)*CalphaTheta(P,b,j)
! B. (P|(beta nu)^x)*Cvirt(beta,b)*Cocc(nu,J)*CalphaTheta(P,b,j)
! C. - (P|Q)^x*Cpq(P,Q) 
! using the following intermediate
! CalphaTheta(P,b,j) = C(P,a,i) * Thetha(a,i,b,j)
! Cpq = CalphaTheta(P,b,j)*C(Q,b,j)
subroutine Build_RIMP2grad(myLSitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
     & CollaborateWithSlaves,Cvirt,nvirt,Cocc,nocc,mynum,numnodes,nAtomsAux,&
     & natoms,ThetaOcc,RIMP2grad,dopair_occ,use_bg_buf)
  implicit none
  type(lsitem), intent(inout) :: mylsitem
  integer,intent(in) :: nAtomsAux,nocc,nvirt,natoms
  integer,intent(in) :: nbasisAux,LUPRI,nbasis,mynum,numnodes
  logical,intent(in) :: master,FORCEPRINT,CollaborateWithSlaves
  logical,intent(in) :: dopair_occ(nocc,nocc),use_bg_buf
  real(realk),intent(in) :: Cocc(nbasis,nocc),Cvirt(nbasis,nvirt)
  real(realk),intent(inout) :: RIMP2grad(3*natoms)
  real(realk),intent(inout) :: ThetaOcc(nvirt*(nocc*i8*nocc)*nvirt)  
  !local variables
  real(realk) :: AlphaBetaDecomp(nbasisAux,nbasisAux),maxsize
  real(realk),pointer :: Calpha(:),CalphaTheta(:),AlphaBetaDeriv(:,:,:)
  real(realk),pointer :: Cpq(:,:),CalphaTmp(:)
  logical :: SymDecomp,AlphaBetaDecompCreate
  integer(kind=8) :: nsize,natoms8,nbasisAux8
  integer :: MynbasisAuxMPI,MYNATOMSMPI,I,Oper,M,N,K,NBA
  integer,pointer :: nbasisAuxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
  integer :: inode,myOriginalRank,MynbasisAuxMPI2,nbasisAuxMPI2(numnodes)
  integer(kind=ls_mpik) :: node 

  nullify(Cpq)
  nullify(Calpha)
  nullify(CalphaTheta)
  nullify(CalphaTmp)
  nullify(AlphaBetaDeriv)
!  real(realk) :: RIMP2gradC(3*natoms)
  natoms8 = natoms
  nbasisAux8 = nbasisAux  
  Oper = CoulombOperator
  SymDecomp = .FALSE.
  AlphaBetaDecompCreate = .TRUE.
  !===========================================================
  !   Determine Sizes1: used to calc 3 center integrals
  !===========================================================
  IF(CollaborateWithSlaves)then 
     !all nodes have info about all nodes 
     call mem_alloc(nbasisAuxMPI,numnodes,'RIGRAD:nbasisAuxMPI')   !number of Aux basis func assigned to rank
     nbasisAuxMPI = 0 
     call mem_alloc(nAtomsMPI,numnodes,'RIGRAD:nAtomsMPI')         !atoms assign to rank
     call mem_alloc(startAuxMPI,nAtomsAux,numnodes,'RIGRAD:startAuxMPI')  !startindex in full (nbasisAux)
     call mem_alloc(AtomsMPI,nAtomsAux,numnodes,'RIGRAD:AtomsMPI') !identity of atoms in full molecule
     call mem_alloc(nAuxMPI,nAtomsAux,numnodes,'RIGRAD:nAuxMPI')   !nauxBasis functions for each of the nAtomsMPI
     IF(DECinfo%AuxAtomicExtent)THEN   
        call getRIbasisMPI(mylsitem%INPUT%AUXMOLECULE,nAtomsAux,numnodes,&
             & nbasisAuxMPI,startAuxMPI,AtomsMPI,nAtomsMPI,nAuxMPI)
     ELSE
        call getRIbasisMPI(mylsitem%SETTING%MOLECULE(1)%p,nAtomsAux,numnodes,&
             & nbasisAuxMPI,startAuxMPI,AtomsMPI,nAtomsMPI,nAuxMPI)
     ENDIF
     MynAtomsMPI = nAtomsMPI(mynum+1)
     MynbasisAuxMPI = nbasisAuxMPI(mynum+1)
     call mem_dealloc(AtomsMPI) !not used in this subroutine 
  ELSE
     MynbasisAuxMPI = nbasisAux     
  ENDIF
  
  !1. build Calpha(P,a,i) = (P|Q)^-1 (Q|ai)
  ! Note this is (nbasisAuxMPI(mynum+1),nvirtAOS,noccEOS) so fairly small quantity 
  ! we assume this fits on node. 
  call Build_CalphaMO(myLSitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
       & CollaborateWithSlaves,Cvirt,nvirt,Cocc,nocc,mynum,numnodes,nAtomsAux,Calpha,&
       & NBA,AlphaBetaDecomp,AlphaBetaDecompCreate,SymDecomp,Oper)
  !2. Construct CalphaTheta(P,j,b) = Calpha(P,a,i) Theta(a,i,j,b)
  !  Size:  (nbasisAuxMPI(mynum+1),noccEOS,nvirtAOS)
  nsize = nvirt*nocc*MynbasisAuxMPI
  IF(use_bg_buf)THEN
     call mem_pseudo_alloc(Cpq,nbasisAux*i8,MynbasisAuxMPI*i8)
     call mem_pseudo_alloc(CalphaTheta,nsize)
  ELSE
     call mem_alloc(CalphaTheta,nsize,'RIGRAD:CalphaTheta')
  ENDIF
  !Build CalphaTheta(P,j,b) =  Calpha(P,a,i)*ThetaOcc(a,i,j,b)
  IF(COUNT(dopair_occ).EQ.nocc*nocc)THEN
     !Atomic fragment
     M =  MynbasisAuxMPI   !rows of Output Matrix
     N =  nocc*nvirt       !columns of Output Matrix
     K =  nvirt*nocc       !summation dimension
     call DGEMM('N','N',M,N,K,1.0E0_realk,Calpha,M,ThetaOcc,K,0.0E0_realk,CalphaTheta,M)     
  ELSE
     !pair fragment more complicated due to dopair_occ
     !only contributions for i on center P and j on center Q and vise versa
     CALL BuildCalphaTheta(nocc,nvirt,MynbasisAuxMPI,Calpha,ThetaOcc,CalphaTheta,dopair_occ)
  ENDIF

  !3. Calculate 2 of 3 gradient contribution: 
  ! A. (P^x|bj)*CalphaTheta(P,b,j)
  ! B. (P|(beta nu)^x)*Cvirt(beta,b)*Cocc(nu,J)*CalphaTheta(P,b,j)
  call II_get_RIMP2_grad(LUPRI,LUPRI,RIMP2grad,Mylsitem%SETTING,&
       & nbasisAux,nbasis,nvirt,nocc,Cvirt,Cocc,mynum,numnodes,natoms,&
       & MynbasisAuxMPI,CalphaTheta,use_bg_buf)

  !This only contains contributions from the Auxiliary functions assinged to this node. 
  !A reduction will performed at the end. 

#ifdef VAR_MPI
!  call lsmpi_reduction(RIMP2grad,3*natoms,infpar%master,infpar%lg_comm)
!  print*,'RIMP2 GRAD A+B'
!  call ls_output(RIMP2grad,1,3,1,natoms,3,natoms,1,6)
#endif

  IF(numnodes.EQ.1)THEN
     ! Build Cpq(P,Q) = CalphaTheta(P,j,b)*Calpha(Q,j,b)  
     IF(.NOT.use_bg_buf)call mem_alloc(Cpq,nbasisAux,nbasisAux,'RIGRAD:Cpq')
     call ConstructCpq(nbasisAux,nocc,nvirt,Calpha,CalphaTheta,Cpq)
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(CalphaTheta)
     ELSE
        call mem_dealloc(CalphaTheta)
     ENDIF
     call mem_dealloc(Calpha)

     !4. Calculate the remianing gradient contribution: CalphaTheta(P,b,j)*(P|Q)^x*Calpha(Q,b,j)
     ! Calculate (P|Q)^x Full Aux
     IF(use_bg_buf)THEN
        !FIXME: test memory if not enough memory have to make special fullcontraction routine. 
        call mem_pseudo_alloc(AlphaBetaDeriv,nbasisAux*i8,nbasisAux*i8,3*natoms*i8)
     ELSE
        call mem_alloc(AlphaBetaDeriv,nbasisAux,nbasisAux,3*natoms,'RIGRAD:AlphaBetaDeriv')
     ENDIF
     call II_get_RI_AlphaBeta_geoderiv2CenterInt(DECinfo%output,DECinfo%output,&
          & AlphaBetaDeriv,mylsitem%setting,nbasisAux,natoms)
     ! Contract grad(x) = grad(x) - (P|Q)^x*Cpq(P,Q) using MPI on outer loop , OpenMP inner loop
     call get_PQ_RIMP2_grad(Cpq,AlphaBetaDeriv,nbasisAux,3*natoms,RIMP2grad)
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(AlphaBetaDeriv)
        call mem_pseudo_dealloc(Cpq)
     ELSE
        call mem_dealloc(AlphaBetaDeriv)
        call mem_dealloc(Cpq)
     ENDIF
  ELSE
#ifdef VAR_MPI
     !We want to construct Cpq(P,Q) = Calpha(P,j,b)*CalphaTheta(Q,j,b)
     !However each node only have part of Calpha(Ppartial,j,b) and same part of 
     !CalphaTheta(Qpartial,j,b) so we broadcast Calpha(Ppartial,j,b)
     !building Cpq(Pfull,Qpartial) = sum_nodes Calpha(PpartialOnnode,j,b)*CalphaTheta(Qpartial,j,b)  
     !this is the part we need to calculate part of 
     !the remaining gradient contribution: Cpq(P,Q)*(P|Q)^x
     !meaning we do Cpq(Pfull,Qpartial)*(Pfull|Qpartial)^x
     !a reduction will give the full contribution. 
     
     IF(.NOT.use_bg_buf)call mem_alloc(Cpq,nbasisAux,MynbasisAuxMPI,'RIGRAD:Cpq2')
     call ls_dzero8(Cpq,nbasisAux*MynbasisAuxMPI*i8)
     DO inode = 1,numnodes
        nsize = nbasisAuxMPI(inode)*nvirt*nocc
        node = inode-1
        IF(mynum.EQ.inode-1)THEN
           !Bcast
           call ls_mpibcast(Calpha,nsize,node,infpar%lg_comm)
           call ConstructCpqMPI(nocc,nvirt,mynum,&
                & numnodes,natoms,MynbasisAuxMPI,&
                & nAtomsMPI,startAuxMPI,nAuxMPI,Calpha,&
                & CalphaTheta,nbasisAux,MynbasisAuxMPI,Cpq)
        ELSE
           IF(use_bg_buf)THEN
              call mem_pseudo_alloc(CalphaTmp,nsize)
           ELSE
              call mem_alloc(CalphaTmp,nsize,'RIGRAD:CalphaTmp')
           ENDIF
           call ls_mpibcast(CalphaTmp,nsize,node,infpar%lg_comm)
           myOriginalRank = inode-1
           call ConstructCpqMPI(nocc,nvirt,myOriginalRank,&
                & numnodes,natoms,nbasisAuxMPI(inode),&
                & nAtomsMPI,startAuxMPI,nAuxMPI,CalphaTmp,&
                & CalphaTheta,nbasisAux,MynbasisAuxMPI,Cpq)
           IF(use_bg_buf)THEN
              call mem_pseudo_dealloc(CalphaTmp)
           ELSE
              call mem_dealloc(CalphaTmp)
           ENDIF
        ENDIF
     ENDDO
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(CalphaTheta)
     ELSE
        call mem_dealloc(CalphaTheta)
     ENDIF
     call mem_dealloc(Calpha)
     !4. Calculate the remianing gradient contribution: CalphaTheta(P,b,j)*(P|Q)^x*Calpha(Q,b,j)
     ! Calculate (P|Q)^x Full Aux
     IF(use_bg_buf)THEN
        !FIXME: test memory if not enough memory have to make special fullcontraction routine. 
        call mem_pseudo_alloc(AlphaBetaDeriv,nbasisAux*i8,nbasisAux*i8,3*natoms*i8)
     ELSE
        call mem_alloc(AlphaBetaDeriv,nbasisAux,nbasisAux,3*natoms,'RIGRAD:AlphaBetaDeriv')
     ENDIF
     call II_get_RI_AlphaBeta_geoderiv2CenterInt(DECinfo%output,DECinfo%output,&
          & AlphaBetaDeriv,mylsitem%setting,nbasisAux,natoms)
     ! Contract grad(x) = grad(x) - (P|Q)^x*Cpq(P,Q) using MPI on outer loop , OpenMP inner loop
     call get_PQ_RIMP2_gradMPI(Cpq,AlphaBetaDeriv,nbasisAux,MynbasisAuxMPI,3*natoms,RIMP2grad,&
          & nAtoms,numnodes,nAtomsMPI,startAuxMPI,nAuxMPI,mynum)

     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(AlphaBetaDeriv)
        call mem_pseudo_dealloc(Cpq)
     ELSE
        call mem_dealloc(AlphaBetaDeriv)
        call mem_dealloc(Cpq)
     ENDIF
     !REDUCTION
     call lsmpi_reduction(RIMP2grad,3*natoms,infpar%master,infpar%lg_comm)
#else
     call lsquit('numnodes not equal to 1, but no MPI',-1)
#endif
  ENDIF

!  IF(mynum.EQ.0)THEN
!     print*,'RIMP2 GRAD A+B+C'
!     call ls_output(RIMP2grad,1,3,1,natoms,3,natoms,1,6)
!  ENDIF

  IF(CollaborateWithSlaves)then 
     call mem_dealloc(nbasisAuxMPI)
     call mem_dealloc(nAtomsMPI)
     call mem_dealloc(startAuxMPI)
     call mem_dealloc(nAuxMPI)
  ENDIF

end subroutine Build_RIMP2grad

subroutine ConstructCpq(nbasisAux,nocc,nvirt,Calpha,CalphaTheta,Cpq)
  implicit none
  integer,intent(in) :: nbasisAux,nocc,nvirt
  real(realk),intent(in) :: Calpha(nbasisAux,nvirt,nocc),CalphaTheta(nbasisAux,nocc,nvirt)
  real(realk),intent(inout) :: Cpq(nbasisAux,nbasisAux)
  !
  integer :: P,Q,I,A
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(P,Q,A,&
  !$OMP I) SHARED(nbasisAux,nocc,nvirt,Calpha,CalphaTheta,Cpq)
  !$OMP DO COLLAPSE(2)
  DO Q=1,nbasisAux
     DO P=1,nbasisAux
        Cpq(P,Q) = Calpha(P,1,1)*CalphaTheta(Q,1,1)
     ENDDO
  ENDDO
  !$OMP END DO 
  DO A=2,nvirt
     !$OMP DO COLLAPSE(2)
     DO Q=1,nbasisAux
        DO P=1,nbasisAux
           Cpq(P,Q) = Cpq(P,Q) + Calpha(P,A,1)*CalphaTheta(Q,1,A)
        ENDDO
     ENDDO
     !$OMP END DO 
  ENDDO
  DO I=2,nocc
     DO A=1,nvirt
        !$OMP DO COLLAPSE(2)
        DO Q=1,nbasisAux
           DO P=1,nbasisAux
              Cpq(P,Q) = Cpq(P,Q) + Calpha(P,A,I)*CalphaTheta(Q,I,A)
           ENDDO
        ENDDO
        !$OMP END DO 
     ENDDO
  ENDDO
  !$OMP END PARALLEL

end subroutine ConstructCpq

subroutine ConstructCpqMPI(nocc,nvirt,myOriginalRank,numnodes,natoms,&
     & nbasisAuxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,CalphaTmp,CalphaTheta,&
     & nbasisAux,MynbasisAuxMPI,Cpq)
  implicit none
  integer,intent(in) :: nocc,nvirt,myOriginalRank,numnodes,natoms,nbasisAuxMPI
  integer,intent(in) :: nAtomsMPI(numnodes),nbasisAux,MynbasisAuxMPI
  integer,intent(in) :: startAuxMPI(nAtoms,numnodes)
  integer,intent(in) :: nAuxMPI(nAtoms,numnodes)
  real(realk),intent(in) :: CalphaTmp(nbasisAuxMPI,nvirt,nocc)
  real(realk),intent(in) :: CalphaTheta(MynbasisAuxMPI,nocc,nvirt)
  real(realk),intent(inout) :: Cpq(nbasisAux,MynbasisAuxMPI)
  !
  integer :: I,A,Q,startB,startB2,iatomB,P,nP
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(I,A,Q,startB,startB2,iatomB,nP,&
  !$OMP P) SHARED(nocc,nvirt,myOriginalRank,numnodes,natoms,nbasisAuxMPI,&
  !$OMP nAtomsMPI,startAuxMPI,nAuxMPI,CalphaTmp,CalphaTheta,Cpq,MynbasisAuxMPI)
  DO I=1,nocc
     DO A=1,nvirt
        startB2 = 0
        DO iAtomB=1,nAtomsMPI(myOriginalRank+1)
           StartB = startAuxMPI(iAtomB,myOriginalRank+1)
           nP = nAuxMPI(iAtomB,myOriginalRank+1)
           !$OMP DO COLLAPSE(2)
           DO Q=1,MynbasisAuxMPI
              DO P = 1,nP
                 Cpq(startB+P,Q) = Cpq(startB+P,Q) + CalphaTmp(startB2+P,A,I)*CalphaTheta(Q,I,A)
              ENDDO
           ENDDO
           !$OMP END DO
           startB2 = startB2 + nP
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL

end subroutine ConstructCpqMPI

subroutine BuildCalphaTheta(nocc,nvirt,MynbasisAuxMPI,Calpha,ThetaOcc,&
     & CalphaTheta,dopair_occ)
  implicit none
  integer,intent(in) :: nocc,nvirt,MynbasisAuxMPI
  real(realk),intent(in) :: Calpha(MynbasisAuxMPI,nvirt,nocc)
  real(realk),intent(in) :: ThetaOcc(nvirt,nocc,nocc,nvirt)
  real(realk),intent(inout) :: CalphaTheta(MynbasisAuxMPI,nocc,nvirt)
  logical,intent(in) :: dopair_occ(nocc,nocc)
  !local variables
  integer :: I,J,A,B,ALPHA
  real(realk) :: TMP(MynbasisAuxMPI),TMP2
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(I,J,A,B,ALPHA,TMP,&
  !$OMP TMP2) SHARED(nocc,nvirt,MynbasisAuxMPI,Calpha,ThetaOcc,CalphaTheta,&
  !$OMP dopair_occ)
  DO B=1,nvirt
     DO J=1,nocc
        DO ALPHA=1,MynbasisAuxMPI
           TMP(ALPHA) = 0.0E0_realk
        ENDDO        
        DO I=1,nocc
           IF(dopair_occ(I,J))THEN
              DO A=1,nvirt
                 TMP2 = ThetaOcc(A,I,J,B)
                 DO ALPHA=1,MynbasisAuxMPI
                    TMP(ALPHA) = TMP(ALPHA) + Calpha(ALPHA,A,I)*TMP2
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
        DO ALPHA=1,MynbasisAuxMPI
           CalphaTheta(ALPHA,J,B) = TMP(ALPHA)
        ENDDO
     ENDDO
  ENDDO
  !OMP END PARALLEL DO
end subroutine BuildCalphaTheta

subroutine get_PQ_RIMP2_gradMPI(Cpq,AlphaBetaDeriv,nbasisAux,MynbasisAuxMPI,ngeoComp,RIMP2grad,&
     & nAtoms,numnodes,nAtomsMPI,startAuxMPI,nAuxMPI,mynum)
implicit none
integer :: nbasisAux,ngeoComp,nAtoms,numnodes,mynum
integer,intent(in) :: nAtomsMPI(numnodes),MynbasisAuxMPI
integer,intent(in) :: startAuxMPI(nAtoms,numnodes)
integer,intent(in) :: nAuxMPI(nAtoms,numnodes)
real(realk),intent(in) :: Cpq(nbasisAux,MynbasisAuxMPI)
real(realk),intent(in) :: AlphaBetaDeriv(nbasisAux,nbasisAux,ngeoComp)
real(realk),intent(inout) :: RIMP2grad(ngeoComp)
!
integer :: P,Q,X,startQ2,startQ,iAtomQ,nQ
real(realk) TMP

!$OMP PARALLEL DEFAULT(none) PRIVATE(P,Q,X,startQ2,startQ,iAtomQ,&
!$OMP nQ) SHARED(Cpq,AlphaBetaDeriv,nbasisAux,MynbasisAuxMPI,ngeoComp,&
!$OMP RIMP2grad,nAtoms,numnodes,nAtomsMPI,startAuxMPI,nAuxMPI,mynum,TMP)
DO X=1,ngeoComp
   !$OMP MASTER
   TMP = 0.0E0_realk
   !$OMP END MASTER
   !$OMP BARRIER 
   startQ2 = 0
   DO iAtomQ=1,nAtomsMPI(mynum+1)
      StartQ = startAuxMPI(iAtomQ,mynum+1) !global startindex
      nQ = nAuxMPI(iAtomQ,mynum+1)
      !$OMP DO COLLAPSE(2) REDUCTION(+:TMP)
      DO Q=1,nQ
         DO P = 1,nbasisAux
            TMP = TMP + Cpq(P,startQ2+Q)*AlphaBetaDeriv(P,startQ+Q,X)   
         ENDDO
      ENDDO
      !$OMP END DO
      startQ2 = startQ2 + nQ   
   ENDDO
   !$OMP MASTER
   RIMP2grad(X) = RIMP2grad(X) - TMP
   !$OMP END MASTER
ENDDO
!$OMP END PARALLEL

end subroutine get_PQ_RIMP2_gradMPI

subroutine get_PQ_RIMP2_grad(Cpq,AlphaBetaDeriv,nbasisAux,ngeoComp,RIMP2grad)
implicit none
integer :: nbasisAux,ngeoComp
real(realk),intent(in) :: Cpq(nbasisAux*nbasisAux),AlphaBetaDeriv(nbasisAux*nbasisAux,ngeoComp)
real(realk),intent(inout) :: RIMP2grad(ngeoComp)
!
integer :: P,X
real(realk) TMP

!$OMP PARALLEL DEFAULT(none) PRIVATE(P,X) SHARED(nbasisAux,ngeoComp,&
!$OMP Cpq,AlphaBetaDeriv,RIMP2grad,TMP)
DO X=1,ngeoComp
!$OMP MASTER
   TMP = 0.0E0_realk
!$OMP END MASTER
!$OMP BARRIER 
!$OMP DO REDUCTION(+:TMP)
   DO P=1,nbasisAux*nbasisAux
      TMP = TMP + Cpq(P)*AlphaBetaDeriv(P,X)   
   ENDDO
!$OMP END DO
!$OMP MASTER
   RIMP2grad(X) = RIMP2grad(X) - TMP
!$OMP END MASTER
ENDDO
!$OMP END PARALLEL

end subroutine get_PQ_RIMP2_grad

subroutine BuilDUmatTmpRIF12(Umat,nAux,UmatTmp,NBA,NBA2,AuxMPIstartMy,iAuxMPIextraMy,&
     & AuxMPIstartMPI,iAuxMPIextraMPI)
  implicit none
  integer,intent(in) :: nAux,NBA,NBA2,AuxMPIstartMy,iAuxMPIextraMy,AuxMPIstartMPI,iAuxMPIextraMPI
  real(realk),intent(in) :: Umat(nAux,nAux)
  real(realk),intent(inout) :: UmatTmp(NBA,NBA2)
  !local variables                              
  integer :: J,I                                
  IF(iAuxMPIextraMy.EQ.0)THEN                                                                                                        
     do J=1,NBA2
        do I=1,NBA
           UmatTmp(I,J) = Umat(AuxMPIstartMy+I,AuxMPIstartMPI+J)
        enddo
     enddo
  ELSE
     do J=1,NBA2
        do I=1,NBA
           UmatTmp(I,J) = Umat(AuxMPIstartMy+I,AuxMPIstartMPI+J)                                                                     
        enddo
        UmatTmp(NBA,J) = Umat(iAuxMPIextraMy,AuxMPIstartMPI+J)                                                                       
     enddo
  ENDIF
  IF(iAuxMPIextraMPI.NE.0)THEN
     IF(iAuxMPIextraMy.EQ.0)THEN                                                                                                     
        do I=1,NBA
           UmatTmp(I,NBA2) = Umat(AuxMPIstartMy+I,iAuxMPIextraMPI)
        enddo
     ELSE
        do I=1,NBA
           UmatTmp(I,NBA2) = Umat(AuxMPIstartMy+I,iAuxMPIextraMPI)                
        enddo
        UmatTmp(NBA,NBA2) = Umat(iAuxMPIextraMy,iAuxMPIextraMPI)                                                                     
     ENDIF
  ENDIF
end subroutine BuilDUmatTmpRIF12

end module ri_util_module

