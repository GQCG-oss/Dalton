!> @file
!> Full calculation
!> This file is mainly a playground for new developments, not intended to be included in a release.

module fullmp2

#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
  use decmpi_module, only: mpi_bcast_fullmolecule
  use lsmpi_op
#endif
  use fundamental
  use precision
  use typedeftype!,only:lsitem
  use typedef
  use dec_typedef_module
  use matrix_module
  use matrix_operations
  use memory_handling
  use MemoryLeakToolMod
  !  DEC DEPENDENCIES (within deccc directory)   
  !  *****************************************
  use dec_tools_module
  use dec_fragment_utils
  use array4_simple_operations
  use array3_simple_operations
  use array2_simple_operations
  use mp2_module
  !  use orbital_operations
  use full_molecule

  public :: full_canonical_mp2, &
       & full_canonical_mpmp2, canonical_mpmp2_memreq_test
  private

contains
  !> \brief Calculate canonical MP2 energy for full molecular system
  !> This is a MP2 version that distribute the integrals
  !> across the nodes. The code does less integral recalculation but 
  !> requires more nodes/more memory
  !> \author Thomas Kjaergaard
  !> \date October 2014
  subroutine full_canonical_mpmp2(MyMolecule,MyLsitem,mp2_energy)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Canonical MP2 correlation energy
    real(realk),intent(inout) :: mp2_energy    
    !
    integer :: nbasis,nocc,nvirt,naux,noccfull,numnodes
    logical :: master,wakeslaves
    real(realk),pointer :: EpsOcc(:),EpsVirt(:)
    integer :: J,I,node,natoms,A,lupri,restart_lun
    logical :: MessageRecieved,MessageRecievedW,FORCEPRINT,file_exists
    real(realk) :: TS,TE,TS2,TE2,TS3,TE3,CPU1,CPU2,WALL1,WALL2,CPU_MPICOMM,WALL_MPICOMM
    real(realk) :: CPU_MPIWAIT,WALL_MPIWAIT,MemInGBCollected,epsIJ,tmp_mp2_energy
    real(realk) :: CPU3,CPU4,WALL3,WALL4,CPU_AOINT,WALL_AOINT,tmp_mp2_energy2
    real(realk) :: CPU_AOTOMO,WALL_AOTOMO,CPU_ECONT,WALL_ECONT,numnodesR
    real(realk),pointer :: tmp1(:),tmp2(:),tmp3(:),tmp4(:)
!    real(realk),pointer :: tmp5(:,:),tmp6(:,:),tmp7(:,:),Amat(:,:),Bmat(:,:)
    real(realk),pointer :: CoBatchA(:,:),CoBatchB(:,:)
!    real(realk),pointer :: CoBI(:,:),CoBJ(:,:),CoI(:,:),CoJ(:,:)
!    real(realk),pointer :: CAV(:,:),CgammaMPI(:,:),tmp1b(:,:,:,:),CgammaMPI2(:,:)
!    real(realk),pointer :: VGVO(:,:),CvA(:,:),CoIG(:,:)
!    real(realk),pointer :: VOVO(:,:),VOVO2(:,:)
    !
    type(DecAObatchinfo),pointer :: Alphabatchinfo(:),Gammabatchinfo(:)
    type(batchtoorb), pointer :: batch2orbAlpha(:),batch2orbGamma(:)
    integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:),DimAOGamma(:),JobGamma(:)
    integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:),DimAOAlpha(:),JobAlpha(:)
    integer, pointer :: batchdimGamma(:),dimGammaArray(:),GammaIndexArray(:),offsetGamma(:,:)
    integer, pointer :: AOstartGamma(:,:),AOendGamma(:,:),AOdimGamma(:,:),AOdimGammaMPI(:)
    integer, pointer :: batchdimAlpha(:),dimAlphaArray(:),AlphaIndexArray(:),offsetAlpha(:,:)
    integer, pointer :: AOstartAlpha(:,:),AOendAlpha(:,:),AOdimAlpha(:,:),AOdimAlphaMPI(:)
    integer, pointer :: nBlocksGamma(:),nBlocksAlpha(:),nOccBatchDimJrank(:),OccIndexJrank(:,:)
    integer :: inodeLoop,jnodeLoop,nodeLoop,nrownodes,ncolnodes,nBlocks,MAXnBlocksGamma,MAXnBlocksAlpha
    integer :: MaxAllowedDimAlpha,MaxAllowedDimGamma,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI
    integer :: MaxActualDimAlpha,MaxActualDimGamma,dimAOoffsetG,dimAOoffsetA
    integer :: nOccBatchDimImax,K,iorb,nbatchesGamma,nbatchesAlpha,nb,inode,jnode,ibatch,offsetG,offsetA
    integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint,M,N,iB
    integer :: MinAObatch,gammaB,alphaB,nOccBatchesI,jB,nOccBatchDimI,nOccBatchDimJ
    integer :: nbatchesGammaRestart,nbatchesAlphaRestart,dimGamma,GammaStart,GammaEnd,dimAlpha
    integer :: AlphaStart,AlphaEnd,B,noccRestartI,noccRestartJ,nJobs,startjB,idxx
    integer :: sqrtnumnodes,gB,idx(1),dimAOoffset,kk,nBlocksG,nBlocksA
    integer :: dimGammaMPI,dimAlphaMPI,ibatchG,ibatchA,dimAlpha2,dimGamma2
    integer :: IMYNUMNBATCHES1,IMYNUMNBATCHES2,nOccBatchDimJmax,offset
    integer :: nOccbatchesIrestart,noccIstart,nbuf1,nbuf2,MAXnOccBatchDimJrank
    logical :: MoTrans, NoSymmetry,SameMol,JobDone,JobInfo1Free,FullRHS,doscreen,NotAllMessagesRecieved
    logical :: PermutationalSymmetryIJ,SetdimGamma
    logical,pointer :: JobsCompleted(:,:)
    integer(kind=8) :: maxsize,Ibuf(12),sizetmp1,sizetmp2,sizetmp3,sizetmp4
    TYPE(DECscreenITEM)   :: DecScreen
    Character            :: intSpec(5)
    integer(kind=ls_mpik)  :: nMPI,TAG,IERR,request,Receiver,sender,comm,TAG1,TAG2
    integer(kind=ls_mpik)  :: request1,request2,masterrank,senderID,mynum,lsmpinode
    !  integer(kind=ls_mpik)  :: Sender
#ifdef VAR_MPI
    integer(kind=ls_mpik)  :: mpistatus(MPI_STATUS_SIZE)
#endif
    integer(kind=4) :: JobInfo1(2)
    !  Character(80)        :: FilenameCS,FilenamePS

    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(full_canonical_mpmp2): does not work with distributed&
       & moleule structure",-1)
    endif

#ifdef VAR_TIME    
    FORCEPRINT = .TRUE.
#else
    IF(LSTIME_PRINT)THEN
       ForcePrint = .TRUE.
    ELSE
       ForcePrint = .FALSE.
    ENDIF
#endif

    CALL LSTIMER('START ',TS,TE,DECINFO%OUTPUT)
    mp2_energy = 0.0E0_realk

    CPU_AOINT = 0.0E0_realk
    WALL_AOINT = 0.0E0_realk
    CPU_AOTOMO = 0.0E0_realk
    WALL_AOTOMO = 0.0E0_realk
    CPU_ECONT = 0.0E0_realk
    WALL_ECONT = 0.0E0_realk

    CPU_MPICOMM = 0.0E0_realk
    WALL_MPICOMM = 0.0E0_realk
    CPU_MPIWAIT = 0.0E0_realk
    WALL_MPIWAIT = 0.0E0_realk
    TAG = 1411; TAG1 = 1412; TAG2 = 1413

    !sanity check
    if(.NOT.DECinfo%use_canonical) then
       call lsquit('Error: full_canonical_mpmp2 require canonical Orbitals',-1)
    endif
    ! Init stuff
    ! **********
    nbasis = MyMolecule%nbasis
    nb = nbasis
    nvirt  = MyMolecule%nvirt
    !MyMolecule%Co is allocated (nbasis,MyMolecule%nocc)
    !with MyMolecule%nocc = Valence + Core 
    !In case of Frozen core we only need Valence and will access
    !Co(nbasis,nval) = MyMolecule%Co(nbasis,MyMolecule%ncore+1:MyMolecule%nocc)
    noccfull = MyMolecule%nocc
    IF(DECinfo%Frozencore)THEN
       nocc   = MyMolecule%nval
       offset = MyMolecule%ncore
    ELSE
       nocc   = MyMolecule%nocc
       offset = 0
    ENDIF
    nAtoms = MyMolecule%nAtoms
    LUPRI = DECinfo%output

#ifdef VAR_MPI
    comm = MPI_COMM_LSDALTON
    master= (infpar%mynum == infpar%master)
    mynum = infpar%mynum
    numnodes = infpar%nodtot
    wakeslaves = infpar%nodtot.GT.1
    IF(.NOT.master)LUPRI = 6 !standard Output
#else
    ! If MPI is not used, consider the single node to be "master"
    master=.true.
    mynum = 0
    numnodes = 1
    wakeslaves = .false.
#endif

#ifdef VAR_MPI 
    ! Master starts up slave
    StartUpSlaves: if(wakeslaves .and. master) then
       ! Wake up slaves to do the job: slaves awoken up with (RIMP2FULL)
       ! and call full_canonical_rimp2_slave which communicate info 
       ! then calls full_canonical_rimp2.
       CALL LS_GETTIM(CPU1,WALL1)
       call ls_mpibcast(CANONMP2FULL,infpar%master,comm)
       ! Communicate fragment information to slaves
       call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
       call mpicopy_lsitem(MyLsitem,comm)
       call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
       call mpi_bcast_fullmolecule(MyMolecule)    
       CALL LS_GETTIM(CPU2,WALL2)
       CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
       WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
    endif StartUpSlaves
#endif

    ! Set integral info
    ! *****************
    !R = Regular Basis set on centers 1-4, C = Coulomb operator
    INTSPEC = ['R','R','R','R','C']

    doscreen = mylsitem%setting%scheme%cs_screen.OR.&
         & mylsitem%setting%scheme%ps_screen

    !determine MinAObatch: the minimum allowed AObatch size + number of AO batches
    IF(DECinfo%useIchor)THEN
       !Determine the minimum allowed AObatch size MinAObatch
       !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
       !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch=6 (the 2*(Px,Py,Pz))
       !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
       !'R'  !Specifies that it is the Regular AO basis that should be batched
       iAO = 4 !the center that the batching should occur on (they are all the same in this case)  
       call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
       iprint = 0           !print level for Ichor Integral code
       MoTrans = .FALSE.    !Do not transform to MO basis! 
       NoSymmetry = .FALSE. !Use Permutational Symmetry! 
       SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
       !Determine the full number of AO batches - not to be confused with the batches of AOs
       !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
       iAO = 1
       call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
       nullify(Gammabatchinfo)
       nullify(Alphabatchinfo)    
    ELSE
       ! The smallest possible AO batch depends on the basis set
       ! (More precisely, if all batches are made as small as possible, then the
       !  call below determines the largest of these small batches).
       call determine_maxBatchOrbitalsize(DECinfo%output,mylsitem%setting,MinAObatch,'R')
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
    ENDIF


    ! ***************************************************************************************
    ! Get optimal values of: MaxAllowedDimAlpha,MaxAllowedDimGamma,nOccBatchDimImax,nOccBatchDimJmax
    ! ***************************************************************************************
#ifdef VAR_MPI
    !use the numbers obtained by master  
    IF(master)THEN
       call get_optimal_batch_sizes_for_canonical_mpmp2(MinAObatch,nbasis,nocc,nvirt,&
            & numnodes,nrownodes,ncolnodes,MaxAllowedDimAlpha,MaxAllowedDimGamma,&
            & MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,nOccBatchDimImax,nOccBatchDimJmax,&
            & sizetmp1,sizetmp2,sizetmp3,sizetmp4)

       Ibuf(1) = nrownodes
       Ibuf(2) = ncolnodes
       Ibuf(3) = MaxAllowedDimAlpha
       Ibuf(4) = MaxAllowedDimGamma
       Ibuf(5) = MaxAllowedDimAlphaMPI
       Ibuf(6) = MaxAllowedDimGammaMPI
       Ibuf(7) = nOccBatchDimImax
       Ibuf(8) = nOccBatchDimJmax
       Ibuf(9) = sizetmp1
       Ibuf(10) = sizetmp2
       Ibuf(11) = sizetmp3
       Ibuf(12) = sizetmp4
    ENDIF
    nbuf1 = 12
    call ls_mpibcast(Ibuf,nbuf1,infpar%master,comm)
    IF(.NOT.master)THEN
       nrownodes = Ibuf(1) 
       ncolnodes = Ibuf(2)
       MaxAllowedDimAlpha = Ibuf(3)
       MaxAllowedDimGamma = Ibuf(4)
       MaxAllowedDimAlphaMPI = Ibuf(5)
       MaxAllowedDimGammaMPI = Ibuf(6)
       nOccBatchDimImax = Ibuf(7)
       nOccBatchDimJmax = Ibuf(8)
       sizetmp1 = Ibuf(9)
       sizetmp2 = Ibuf(10)
       sizetmp3 = Ibuf(11)
       sizetmp4 = Ibuf(12)
    ENDIF
#else
    call get_optimal_batch_sizes_for_canonical_mpmp2(MinAObatch,nbasis,nocc,nvirt,&
         & numnodes,nrownodes,ncolnodes,MaxAllowedDimAlpha,MaxAllowedDimGamma,&
         & MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,nOccBatchDimImax,nOccBatchDimJmax,&
         & sizetmp1,sizetmp2,sizetmp3,sizetmp4)
#endif

    IF(master)THEN
       write(DECinfo%output,*)'nbasis               ',nbasis
       write(DECinfo%output,*)'nocc                 ',nocc
       write(DECinfo%output,*)'nvirt                ',nvirt
       write(DECinfo%output,*)'MinAObatch           ',MinAObatch
       write(DECinfo%output,*)'numnodes             ',numnodes
       write(DECinfo%output,*)'nrownodes            ',nrownodes
       write(DECinfo%output,*)'ncolnodes            ',ncolnodes
       write(DECinfo%output,*)'MaxAllowedDimAlpha   ',MaxAllowedDimAlpha
       write(DECinfo%output,*)'MaxAllowedDimGamma   ',MaxAllowedDimGamma
       write(DECinfo%output,*)'MaxAllowedDimAlphaMPI',MaxAllowedDimAlphaMPI
       write(DECinfo%output,*)'MaxAllowedDimGammaMPI',MaxAllowedDimGammaMPI
       write(DECinfo%output,*)'nOccBatchDimImax     ',nOccBatchDimImax
       write(DECinfo%output,*)'nOccBatchDimJmax     ',nOccBatchDimJmax
       write(DECinfo%output,*)'sizetmp1             ',sizetmp1,' = ',sizetmp1*8.0E-9_realk,' GB'
       write(DECinfo%output,*)'sizetmp2             ',sizetmp2,' = ',sizetmp2*8.0E-9_realk,' GB'
       write(DECinfo%output,*)'sizetmp3             ',sizetmp3,' = ',sizetmp3*8.0E-9_realk,' GB'
       write(DECinfo%output,*)'sizetmp4             ',sizetmp4,' = ',sizetmp4*8.0E-9_realk,' GB'
    ELSE
       print*,'SLAVE nbasis               ',nbasis
       print*,'SLAVE nocc                 ',nocc
       print*,'SLAVE nvirt                ',nvirt
       print*,'SLAVE MinAObatch           ',MinAObatch
       print*,'SLAVE numnodes             ',numnodes
       print*,'SLAVE nrownodes            ',nrownodes
       print*,'SLAVE ncolnodes            ',ncolnodes
       print*,'SLAVE MaxAllowedDimAlpha   ',MaxAllowedDimAlpha
       print*,'SLAVE MaxAllowedDimGamma   ',MaxAllowedDimGamma
       print*,'SLAVE MaxAllowedDimAlphaMPI',MaxAllowedDimAlphaMPI
       print*,'SLAVE MaxAllowedDimGammaMPI',MaxAllowedDimGammaMPI
       print*,'SLAVE nOccBatchDimImax     ',nOccBatchDimImax
       print*,'SLAVE nOccBatchDimJmax     ',nOccBatchDimJmax
       print*,'SLAVE sizetmp1             ',sizetmp1,' = ',sizetmp1*8.0E-9_realk,' GB'
       print*,'SLAVE sizetmp2             ',sizetmp2,' = ',sizetmp2*8.0E-9_realk,' GB'
       print*,'SLAVE sizetmp3             ',sizetmp3,' = ',sizetmp3*8.0E-9_realk,' GB'
       print*,'SLAVE sizetmp4             ',sizetmp4,' = ',sizetmp4*8.0E-9_realk,' GB'
    ENDIF

    write(DECinfo%output,*)'Allocating tmp1'; call flush(DECinfo%output)
    call mem_alloc(tmp1,sizetmp1)
    write(DECinfo%output,*)'Allocating tmp2'; call flush(DECinfo%output)
    call mem_alloc(tmp2,sizetmp2)
    write(DECinfo%output,*)'Allocating tmp3'; call flush(DECinfo%output)
    call mem_alloc(tmp3,sizetmp3)
    write(DECinfo%output,*)'Allocating tmp4'; call flush(DECinfo%output)
    call mem_alloc(tmp4,sizetmp4)

    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! * And 
    ! * Determine batch information for Alpha batch  *
    ! ************************************************

    IF(DECinfo%useIchor)THEN
       iAO = 2 !Gamma is the 2. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the bat%MaxAllowedDimGamma, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
            & nbatchesGamma,DECinfo%output)
       call mem_alloc(Gammabatchinfo,nbatchesGamma)
       !Construct the batches of AOS based on the bat%MaxAllowedDimGamma, the requested
       !size of the AO batches - bat%MaxAllowedDimGamma must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimGamma must be less og equal to bat%MaxAllowedDimGamma
       call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
            & nbatchesGamma,Gammabatchinfo,MaxActualDimGamma,DECinfo%output)

       iAO = 1 !Alpha is the 1. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the bat%MaxAllowedDimAlpha, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
            & nbatchesAlpha,DECinfo%output)
       call mem_alloc(Alphabatchinfo,nbatchesAlpha)
       !Construct the batches of AOS based on the bat%MaxAllowedDimAlpha, the requested
       !size of the AO batches - bat%MaxAllowedDimAlpha must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimAlpha must be less og equal to bat%MaxAllowedDimAlpha
       call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
            & nbatchesAlpha,Alphabatchinfo,MaxActualDimAlpha,DECinfo%output)
    ELSE
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchGamma,nbasis)
       call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
            & nbasis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,&
            & orb2BatchGamma,'R')
       ! Translate batchindex to orbital index
       ! -------------------------------------
       call mem_alloc(batch2orbGamma,nbatchesGamma)
       do idxx=1,nbatchesGamma
          call mem_alloc(batch2orbGamma(idxx)%orbindex,batchdimGamma(idxx) )
          batch2orbGamma(idxx)%orbindex = 0
          batch2orbGamma(idxx)%norbindex = 0
       end do
       do iorb=1,nbasis
          idxx = orb2batchGamma(iorb)
          batch2orbGamma(idxx)%norbindex = batch2orbGamma(idxx)%norbindex+1
          K = batch2orbGamma(idxx)%norbindex
          batch2orbGamma(idxx)%orbindex(K) = iorb
       end do
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchAlpha,nbasis)
       call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
            & nbasis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,&
            & orb2BatchAlpha,'R')
       ! Translate batchindex to orbital index
       ! -------------------------------------
       call mem_alloc(batch2orbAlpha,nbatchesAlpha)
       do idxx=1,nbatchesAlpha
          call mem_alloc(batch2orbAlpha(idxx)%orbindex,batchdimAlpha(idxx) )
          batch2orbAlpha(idxx)%orbindex = 0
          batch2orbAlpha(idxx)%norbindex = 0
       end do
       do iorb=1,nbasis
          idxx = orb2batchAlpha(iorb)
          batch2orbAlpha(idxx)%norbindex = batch2orbAlpha(idxx)%norbindex+1
          K = batch2orbAlpha(idxx)%norbindex
          batch2orbAlpha(idxx)%orbindex(K) = iorb
       end do
    ENDIF

    IF(master)THEN
       write(DECinfo%output,*)'MaxActualDimAlpha    ',MaxActualDimAlpha
       write(DECinfo%output,*)'MaxActualDimGamma    ',MaxActualDimGamma
    ENDIF

    ! ************************************************
    ! * Screening                                    *
    ! ************************************************
    CALL LS_GETTIM(CPU3,WALL3)
    IF(DECinfo%useIchor)THEN
       !Calculate Screening integrals 
       SameMOL = .TRUE. !Specifies same molecule on all centers 
       call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
    ELSE
       ! This subroutine builds the full screening matrix.
       call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
            & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%IntegralThreshold)
       IF(doscreen)THEN
          call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
               & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
               & batchindexAlpha,batchindexGamma,&
               & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
       ENDIF
    ENDIF
    CALL LS_GETTIM(CPU4,WALL4)
    CPU_AOINT = CPU_AOINT + (CPU4-CPU3)
    WALL_AOINT = WALL_AOINT + (WALL4-WALL3)

    !GAMMA SET (ncolnodes)
    call mem_alloc(dimGammaArray,nbatchesGamma)
    call mem_alloc(GammaIndexArray,nbatchesGamma)
    !determine which node have which gammaB batch 
    do gammaB = 1,nbatchesGamma  ! AO batches
       IF(DECinfo%useIchor)THEN
          dimGamma = Gammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
       ELSE
          dimGamma = batchdimGamma(gammaB)           ! Dimension of gamma batch
       ENDIF
       dimGammaArray(gammaB) = dimGamma
       GammaIndexArray(gammaB) = gammaB
    enddo
    call integer_inv_sort_with_tracking(dimGammaArray,GammaIndexArray,nbatchesGamma)
    call mem_alloc(DimAOGamma,ncolnodes)
    call mem_alloc(JobGamma,nbatchesGamma)
    call mem_alloc(nBlocksGamma,ncolnodes)
    DimAOGamma = 0
    nBlocksGamma = 0
    do gB = 1,nbatchesGamma  ! AO batches
       gammaB = GammaIndexArray(gB)
       dimGamma = dimGammaArray(gB)
       idx = MINLOC(DimAOGamma)
       DimAOGamma(idx(1)) = DimAOGamma(idx(1)) + DimGamma
       JobGamma(gammaB) = idx(1)
       nBlocksGamma(idx(1)) = nBlocksGamma(idx(1)) + 1
    enddo
    call mem_dealloc(DimAOGamma)
    call mem_dealloc(dimGammaArray)
    call mem_dealloc(GammaIndexArray)

    !  print*,'nBlocksGamma',nBlocksGamma,'MYNUM',MYNUM
    MAXnBlocksGamma = MAXVAL(nBlocksGamma)

    call mem_alloc(offsetGamma,MAXnBlocksGamma,ncolnodes)
    offsetGamma=0
    call mem_alloc(AOstartGamma,MAXnBlocksGamma,ncolnodes)
    AOstartGamma=0
    call mem_alloc(AOendGamma,MAXnBlocksGamma,ncolnodes)
    AOendGamma=0
    call mem_alloc(AOdimGamma,MAXnBlocksGamma,ncolnodes)
    AOdimGamma=0
    call mem_alloc(AOdimGammaMPI,ncolnodes)
    AOdimGammaMPI=0

    do jnode = 1,ncolnodes
       dimAOoffsetG = 0 
       nBlocks = 0
       do gammaB = 1,nbatchesGamma
          IF(JobGamma(gammaB).EQ.jnode)THEN 
             nBlocks = nBlocks + 1
             offsetGamma(nBlocks,jnode) = dimAOoffsetG
             IF(DECinfo%useIchor)THEN
                dimGamma = Gammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
                GammaStart = Gammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
                GammaEnd = Gammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
             ELSE
                dimGamma = batchdimGamma(gammaB)                     ! Dimension of gamma batch
                GammaStart = batch2orbGamma(gammaB)%orbindex(1)      ! First index in gamma batch
                GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma) ! Last index in gamma batch
             ENDIF
             IF(dimAOoffsetG + dimGamma.GT.MaxAllowedDimGammaMPI)&
                  & call lsquit('dimAOoffsetG + dimGamma.GT.MaxAllowedDimGammaMPI',-1)
             dimAOoffsetG = dimAOoffsetG + dimGamma
             AOstartGamma(nBlocks,jnode) = GammaStart
             AOendGamma(nBlocks,jnode) = GammaEnd
             AOdimGamma(nBlocks,jnode) = dimGamma
             AOdimGammaMPI(jnode) =  dimAOoffsetG
          endif
       enddo
    enddo
    !  do jnode = 1,ncolnodes
    !     print*,'AOstartGamma(:,jnode=',jnode,')',AOstartGamma(1:nBlocks,jnode)
    !     print*,'AOendGamma(:,jnode=',jnode,')',AOendGamma(1:nBlocks,jnode)
    !     print*,'AOdimGamma(:,jnode=',jnode,')',AOdimGamma(1:nBlocks,jnode)
    !     print*,'AOdimGammaMPI(jnode=',jnode,')',AOdimGammaMPI(jnode)
    !  enddo

    call mem_alloc(dimAlphaArray,nbatchesAlpha)
    call mem_alloc(AlphaIndexArray,nbatchesAlpha)
    !ALPHA SET (nrownodes)
    !determine which node have which alphaB batch 
    do alphaB = 1,nbatchesAlpha  ! AO batches
       IF(DECinfo%useIchor)THEN
          dimAlpha = Alphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
       ELSE
          dimAlpha = batchdimAlpha(alphaB)           ! Dimension of alpha batch
       ENDIF
       dimAlphaArray(alphaB) = dimAlpha
       AlphaIndexArray(alphaB) = alphaB
    enddo
    call integer_inv_sort_with_tracking(dimAlphaArray,AlphaIndexArray,nbatchesAlpha)
    call mem_alloc(DimAOAlpha,nrownodes)
    call mem_alloc(JobAlpha,nbatchesAlpha)
    call mem_alloc(nBlocksAlpha,nrownodes)
    DimAOAlpha = 0
    nBlocksALpha = 0
    do gB = 1,nbatchesAlpha  ! AO batches
       alphaB = AlphaIndexArray(gB)
       dimAlpha = dimAlphaArray(gB)
       idx = MINLOC(DimAOAlpha)
       DimAOAlpha(idx(1)) = DimAOAlpha(idx(1)) + DimAlpha
       JobAlpha(alphaB) = idx(1)
       nBlocksAlpha(idx(1)) = nBlocksAlpha(idx(1)) + 1
    enddo
    DO gB=1,nrownodes
       IF(DimAOAlpha(gB).GT.MaxAllowedDimAlphaMPI)THEN
          print*,'Warning: DimAOAlpha.GT.MaxAllowedDimAlphaMPI'
          print*,'Warning: DimAOAlpha=',DimAOAlpha(gB)
          print*,'Warning: MaxAllowedDimAlphaMPI=',MaxAllowedDimAlphaMPI
       ENDIF
    ENDDO
    call mem_dealloc(DimAOAlpha)
    call mem_dealloc(dimAlphaArray)
    call mem_dealloc(AlphaIndexArray)

    !  print*,'nBlocksAlpha',nBlocksAlpha,'MYNUM',MYNUM
    MAXnBlocksAlpha = MAXVAL(nBlocksAlpha)

    call mem_alloc(offsetAlpha,MAXnBlocksAlpha,nrownodes)
    offsetAlpha=0
    call mem_alloc(AOstartAlpha,MAXnBlocksAlpha,nrownodes)
    AOstartAlpha=0
    call mem_alloc(AOendAlpha,MAXnBlocksAlpha,nrownodes)
    AOendAlpha=0
    call mem_alloc(AOdimAlpha,MAXnBlocksAlpha,nrownodes)
    AOdimAlpha=0
    call mem_alloc(AOdimAlphaMPI,nrownodes)
    AOdimAlphaMPI=0

    do inode = 1,nrownodes
       dimAOoffsetA = 0 
       nBlocks = 0
       do alphaB = 1,nbatchesAlpha
          IF(JobAlpha(alphaB).EQ.inode)THEN 
             nBlocks = nBlocks + 1
             offsetAlpha(nBlocks,inode) = dimAOoffsetA
             IF(DECinfo%useIchor)THEN
                dimAlpha = Alphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
                AlphaStart = Alphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
                AlphaEnd = Alphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
             ELSE
                dimAlpha = batchdimAlpha(alphaB)                ! Dimension of alpha batch
                AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)            ! First index in alpha batch
                AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)       ! Last index in alpha batch
             ENDIF
             IF(dimAOoffsetA + dimAlpha.GT.MaxAllowedDimAlphaMPI)THEN
!                call lsquit('dimAOoffsetA + dimAlpha.GT.MaxAllowedDimAlphaMPI',-1)
                print*,'Warning: dimAOoffsetA + dimAlpha.GT.MaxAllowedDimAlphaMPI'
             ENDIF
             dimAOoffsetA = dimAOoffsetA + dimAlpha
             AOstartAlpha(nBlocks,inode) = AlphaStart
             AOendAlpha(nBlocks,inode) = AlphaEnd
             AOdimAlpha(nBlocks,inode) = dimAlpha
             AOdimAlphaMPI(inode) =  dimAOoffsetA
          endif
       enddo
    enddo
    !  do inode = 1,nrownodes
    !     print*,'AOstartAlpha(:,inode=',inode,')',AOstartAlpha(1:nBlocks,inode)
    !     print*,'AOendAlpha(:,inode=',inode,')',AOendAlpha(1:nBlocks,inode)
    !     print*,'AOdimAlpha(:,inode=',inode,')',AOdimAlpha(1:nBlocks,inode)
    !     print*,'AOdimAlphaMPI(inode=',inode,')',AOdimAlphaMPI(inode)
    !  enddo

    inode = mod(mynum,nrownodes)+1   !0 => 1, 1 => 2, 2 => 1, 3 => 2
    jnode = (mynum)/(nrownodes) + 1  !0 => 1, 1 => 1, 2 => 2, 3 => 2
    IF(inode+(jnode-1)*nrownodes.NE.mynum+1)call lsquit('node mismatch in MP2',-1)

    dimAlphaMPI = AOdimAlphaMPI(inode)
    dimGammaMPI = AOdimGammaMPI(jnode)
    !building
    !nOccBatchDimJrank(mynum)
    !OccIndexJrank(J,mynum)              
    IF(numnodes.GT.nocc)THEN
       call lsquit('Error more nodes then occupied orbitals: FIXME',-1)
    ENDIF

    call mem_alloc(nOccBatchDimJrank,ncolnodes)
    nOccBatchDimJrank = 0
    do J=1,nocc
       idx = MINLOC(nOccBatchDimJrank)
       nOccBatchDimJrank(idx(1)) = nOccBatchDimJrank(idx(1)) + 1
    enddo
    IF(MAXVAL(nOccBatchDimJrank).GT.nOccBatchDimJmax.OR.MAXVAL(nOccBatchDimJrank).LT.0)THEN
       print*,'nOccBatchDimJrank',nOccBatchDimJrank
       print*,'nOccBatchDimJmax',nOccBatchDimJmax
       print*,'MAXVAL(nOccBatchDimJrank)',MAXVAL(nOccBatchDimJrank)
       call lsquit('miscalc nOccBatchDimJrank full canon mpmp2 ',-1)
    ENDIF
    MAXnOccBatchDimJrank = MAXVAL(nOccBatchDimJrank)
    call mem_alloc(OccIndexJrank,MAXnOccBatchDimJrank,ncolnodes)   
    nOccBatchDimJrank = 0
    do J=1,nocc
       idx = MINLOC(nOccBatchDimJrank)
       nOccBatchDimJrank(idx(1)) = nOccBatchDimJrank(idx(1)) + 1
       OccIndexJrank(nOccBatchDimJrank(idx(1)),idx(1)) = J
    enddo

    nOccBatchesI = nOcc/nOccBatchDimImax
    IF(MOD(nOcc,nOccBatchDimImax).NE.0)nOccBatchesI = nOccBatchesI + 1  
    PermutationalSymmetryIJ = .FALSE.
    write(DECinfo%output,*)'nOccBatchesI         ',nOccBatchesI
    !IF(nrownodes.NE.nOccBatchesI)call lsquit('error in nocbatches in mp2',-1)

    call mem_alloc(EpsOcc,nocc)
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
    !$OMP SHARED(noccfull,MyMolecule,EpsOcc,offset)
    do I=1+offset,noccfull
       EpsOcc(I-offset) = MyMolecule%oofock%elm2(I,I)
    enddo
    !$OMP END PARALLEL DO
    call mem_alloc(EpsVirt,nvirt)
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
    !$OMP SHARED(nvirt,MyMolecule,EpsVirt)
    do A=1,nvirt
       EpsVirt(A) = MyMolecule%vvfock%elm2(A,A)
    enddo
    !$OMP END PARALLEL DO

    ! ***************************************************************************************
    ! Restart
    ! ***************************************************************************************
    IF(DECinfo%DECrestart)THEN
       !CHECK IF THERE ARE ENERGY CONTRIBUTIONS AVAILABLE
       INQUIRE(FILE='FULLMPMP2.restart',EXIST=file_exists)
       IF(file_exists)THEN
          IF(master)THEN
             WRITE(DECinfo%output,*)'Restart of Full molecular MP2(MPMP2) calculation:'
          ENDIF
          restart_lun = -1  !initialization
          call lsopen(restart_lun,'FULLMPMP2.restart','OLD','FORMATTED')
          rewind restart_lun
          read(restart_lun,'(I9)') nOccbatchesIrestart
          read(restart_lun,'(I9)') noccIstart
          IF(nOccbatchesIrestart.NE.nOccbatchesI)THEN
             print*,'Restart Error: nOccbatchesIrestart=',nOccbatchesIrestart
             print*,'Restart Error: nOccbatchesI       =',nOccbatchesI
             call lsquit('MP2 restart error first integer is wrong')
          ELSE
             IF(noccIstart.EQ.nOccbatchesI)THEN
                IF(master)WRITE(DECinfo%output,*)'All energies is on file'
                noccIstart = nocc+1
                read(restart_lun,'(F28.16)') mp2_energy
             ELSEIF(noccIstart.GT.nOccbatchesI.OR.noccIstart.LT.1)THEN
                IF(master)THEN
                   WRITE(DECinfo%output,*)'MP2 restart error, second integer is wrong. Read:',noccIstart
                ENDIF
                call lsquit('MP2 restart error second integer is wrong')             
             ELSE
                noccIstart = noccIstart + 1
                read(restart_lun,'(F28.16)') mp2_energy
             ENDIF
          ENDIF
          call lsclose(restart_lun,'KEEP')
       ELSE
          noccIstart=1
          mp2_energy = 0.0E0_realk
       ENDIF
    ELSE
       noccIstart=1
       mp2_energy = 0.0E0_realk
    ENDIF

    ! ************************************************
    ! * Main Loop                                    *
    ! ************************************************
    JobInfo1Free = .FALSE.
    FullRHS = (nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)

    BatchOccI: do iB = noccIstart,nOccbatchesI
       nOccBatchDimI = nOccBatchDimImax
       IF(MOD(nOcc,nOccBatchDimI).NE.0.AND.iB.EQ.nOccBatchesI)THEN
          !the remainder
          nOccBatchDimI = MOD(nOcc,nOccBatchDimImax)
       ENDIF

       !construct tmp4(nb,nOccBatchDimI)
!       call mem_alloc(tmp4,nb*nOccBatchDimI)
!       IF(nb*nOccBatchDimI.GT.size(tmp4))CALL LSQUIT('TEST1 tmp4',-1)
       call buildCoIMPMP2(tmp4,nb,nOccBatchDimI,nocc,Mymolecule%Co%elm2,&
            & offset,iB,nOccBatchDimImax)
!       IF(dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI.GT.size(tmp2))CALL LSQUIT('TEST1 tmp2',-1)
!       call mem_alloc(tmp2,dimAlphaMPI,dimGammaMPI,nb,nOccBatchDimI)
       !MemoryBookkeeping: tmp2,tmp4 = dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI+nb*nOccBatchDimI
       !dgemm:   tmp2(dimAlpha,dimGamma,nbasis,noccB)        
       nBlocksG = 0 
       BatchGamma: do gammaB = 1,nbatchesGamma
          nBlocksA = 0 
          BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
             IF(JobAlpha(alphaB).EQ.inode.AND.JobGamma(gammaB).EQ.jnode)THEN 
                IF(nBlocksA.EQ.0) nBlocksG = nBlocksG + 1
                nBlocksA = nBlocksA + 1

                offsetA = offsetAlpha(nBlocksA,inode)
                offsetG = offsetGamma(nBlocksG,jnode)
                IF(DECinfo%useIchor)THEN
                   dimGamma = Gammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
                   GammaStart = Gammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
                   GammaEnd = Gammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
                   AOGammaStart = Gammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
                   AOGammaEnd = Gammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
                   dimAlpha = Alphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
                   AlphaStart = Alphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
                   AlphaEnd = Alphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
                   AOAlphaStart = Alphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
                   AOAlphaEnd = Alphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
                ELSE
                   dimGamma = batchdimGamma(gammaB)                      ! Dimension of gamma batch
                   GammaStart = batch2orbGamma(gammaB)%orbindex(1)       ! First index in gamma batch
                   GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)  ! Last index in gamma batch
                   dimAlpha = batchdimAlpha(alphaB)                      ! Dimension of alpha batch
                   AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)       ! First index in alpha batch
                   AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)  ! Last index in alpha batch
                ENDIF

                CALL LS_GETTIM(CPU4,WALL4)
!                IF(dimAlpha*dimGamma*nb*nb.GT.size(tmp1))CALL LSQUIT('TEST1 tmp1',-1)

!                call mem_alloc(tmp1,dimAlpha*dimGamma,nb*nb)
                !MemoryBookkeeping: tmp2,tmp4,tmp1 = dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI
                !+nb*nOccBatchDimI+MaxdimAlpha*MaxdimGamma*nb*nb
                IF(DECinfo%useIchor)THEN
                   !(dimAlpha,dimGamma,nb,nb)
                   call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,dimAlpha,dimGamma,nb,nb,&
                        & tmp1,INTSPEC,FULLRHS,AOAlphaStart,AOAlphaEnd,AOGammaStart,AOGammaEnd,&
                        & 1,nAObatches,1,nAObatches,MoTrans,dimAlpha,dimGamma,nb,nb,NoSymmetry,DECinfo%IntegralThreshold)
                   CALL LS_GETTIM(CPU3,WALL3)
                   CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
                   WALL_AOINT = WALL_AOINT + (WALL3-WALL4)
                ELSE
                   IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%masterGabRHS
                   IF(doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p
                   call II_GET_DECPACKED4CENTER_J_ERI2(DECinfo%output,DECinfo%output, &
                        & mylsitem%setting,tmp1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
                        & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),dimAlpha,&
                        & dimGamma,nbasis,nbasis,FullRHS,INTSPEC,DECinfo%IntegralThreshold)
                ENDIF
                CALL LS_GETTIM(CPU3,WALL3)
                CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
                WALL_AOINT = WALL_AOINT + (WALL3-WALL4)

                !tmp3(dimAlpha,dimGamma,nb,nOccBatchDimI)=tmp1(dimAlpha,dimGamma,nb,nb)*tmp4(nb,nOccBatchDimI)
!                IF(dimAlpha*dimGamma*nb*nOccBatchDimI.GT.size(tmp1))CALL LSQUIT('TEST1 tmp3',-1)
!                call mem_alloc(tmp3,dimAlpha,dimGamma,nb,nOccBatchDimI)
                !MemoryBookkeeping: tmp2,tmp4,tmp1,tmp3 = dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI
                !+nb*nOccBatchDimI+MaxdimAlpha*MaxdimGamma*nb*nb
                !+MaxdimAlpha*MaxdimGamma*nb*nOccBatchDimI 

                M = dimAlpha*dimGamma*nb     !rows of Output Matrix
                N = nOccBatchDimI            !columns of Output Matrix
                K = nb                       !summation dimension
                call dgemm('N','N',M,N,K,1.0E0_realk,tmp1,M,tmp4,K,0.0E0_realk,tmp3,M)
!                call mem_dealloc(tmp1)
                call AddIntegralToCollectMPMP2(tmp2,tmp3,nOccBatchDimI,nb,dimGamma,dimAlpha,&
                     & dimAlphaMPI,dimGammaMPI,offsetA,offsetG)
!                call mem_dealloc(tmp3)
                CALL LS_GETTIM(CPU4,WALL4)
                CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
                WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
             ENDIF
          enddo BatchAlpha
       enddo BatchGamma
       CALL LS_GETTIM(CPU3,WALL3)
!       call mem_dealloc(tmp4)

       !reorder: tmp3(nb,nOccBatchDimJ,dimAlphaMPI,dimGammaMPI) = tmp2(dimAlphaMPI,dimGammaMPI,nb,nOccBatchDimI)
!       IF(nb*nOccBatchDimI*i8*dimAlphaMPI*dimGammaMPI.GT.size(tmp1))CALL LSQUIT('TEST2 tmp3',-1)
!       call mem_alloc(tmp3,nb*nOccBatchDimI*i8*dimAlphaMPI*dimGammaMPI)
       !MemoryBookkeeping: tmp2 ,tmp3= dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI
       !+nb*nOccBatchDimI*dimAlphaMPI*dimGammaMPI
       M = dimAlphaMPI*dimGammaMPI   !row of Input Matrix
       N = nb*nOccBatchDimI    !columns of Input Matrix
       call mat_transpose(M,N,1.0E0_realk,tmp2,0.0E0_realk,tmp3)
!       call mem_dealloc(tmp2)

       !tmp2(nvirt,nOccBatchDimI,dimAlphaMPI,dimGammaMPI) = Cv(nb,nvirt)*tmp3(nb,nOccBatchDimI,dimAlphaMPI,dimGammaMPI)
!       IF(nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI.GT.size(tmp2))CALL LSQUIT('TEST2 tmp2',-1)
!       call mem_alloc(tmp2,nvirt*nOccBatchDimI,dimAlphaMPI*dimGammaMPI)
       !MemoryBookkeeping: tmp3,tmp2 = nb*nOccBatchDimI*dimAlphaMPI*dimGammaMPI
       !+ nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI
       M = nvirt                                 !rows of Output Matrix
       N = nOccBatchDimI*dimAlphaMPI*dimGammaMPI !columns of Output Matrix
       K = nb                                    !summation dimension
       call dgemm('T','N',M,N,K,1.0E0_realk,MyMolecule%Cv%elm2,K,tmp3,K,0.0E0_realk,tmp2,M)
!       call mem_dealloc(tmp3)

       nOccBatchDimJ = nOccBatchDimJrank(JNODE)

!       call mem_alloc(tmp3,nvirt*nOccBatchDimI,dimAlphaMPI*nOccBatchDimJ)
       !MemoryBookkeeping: tmp2,tmp3 = nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI
       !+ nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ
       tmp3 = 0.0E0_realk
       CALL LS_GETTIM(CPU4,WALL4)
       CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
       WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)

       do jnodeLoop = 1,ncolnodes
          do inodeLoop = 1,nrownodes !alpha 
             nodeLoop = inodeLoop+(jnodeLoop-1)*nrownodes
             IF(nodeLoop.EQ.mynum+1)THEN
#ifdef VAR_MPI
                CALL LS_GETTIM(CPU1,WALL1)
                nbuf1 = nvirt*nOccBatchDimI
                nbuf2 = dimAlphaMPI*dimGammaMPI
                call ls_mpibcast(tmp2,nbuf1*i8*nbuf2,mynum,comm)
                CALL LS_GETTIM(CPU2,WALL2)
                CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
                WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
#endif
                !construct CgammaMPI(dimGammaMPI,nOccBatchDimI)
                IF(jnodeLoop.NE.jnode)call lsquit('jnodeLoop error',-1)
                CALL LS_GETTIM(CPU3,WALL3)
!                call mem_alloc(tmp1,dimGammaMPI,nOccBatchDimJ)
                !MemoryBookkeeping: tmp2,tmp3,tmp1 = nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI
                !+ nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ + dimGammaMPI*nOccBatchDimJ 
                call BuildCgammaMPIMPMP2(nOccBatchDimJ,MAXnOccBatchDimJrank,ncolnodes,OccIndexJrank,&
                     & MAXnBlocksGamma,nBlocksGamma,offsetGamma,AOdimGamma,AOstartGamma,AOendGamma,&
                     & Mymolecule%Co%elm2,nb,nocc,dimGammaMPI,tmp1,jnode,jnodeloop,offset)
                !tmp3(nvirt,noccBI,dimAlpha,noccBJ)=tmp2(nvirt,noccI,dimAlpha,dimGammaMPI)*C(dimGammaMPI,noccBJ)
                M = nvirt*nOccBatchDimI*dimAlphaMPI !rows of Output Matrix
                N = nOccBatchDimJ  !columns of Output Matrix
                K = dimGammaMPI    !summation dimension
                call dgemm('N','N',M,N,K,1.0E0_realk,tmp2,M,tmp1,K,1.0E0_realk,tmp3,M)
!                call mem_dealloc(tmp1)
                CALL LS_GETTIM(CPU4,WALL4)
                CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
                WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
             ELSE
                dimAlpha2 = AOdimAlphaMPI(inodeLoop)
                dimGamma2 = AOdimGammaMPI(jnodeLoop)
 !               call mem_alloc(tmp1,nvirt*nOccBatchDimI,dimAlpha2*dimGamma2)
                !MemoryBookkeeping: tmp2,tmp3,tmp1 = nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI
                !+ nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ 
                !+ nvirt*nOccBatchDimI,Max(dimAlpha2)*Max(dimGamma2)

                !recv
#ifdef VAR_MPI
                CALL LS_GETTIM(CPU1,WALL1)
                nbuf1 = nvirt*nOccBatchDimI
                nbuf2 = dimAlpha2*dimGamma2
                lsmpinode = nodeLoop-1
                call ls_mpibcast(tmp1,nbuf1*nbuf2*i8,lsmpinode,comm)
                CALL LS_GETTIM(CPU2,WALL2)
                CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
                WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
#endif
                CALL LS_GETTIM(CPU3,WALL3)
                IF(inodeLoop.EQ.inode)THEN
                   !same Alpha Batch - receiving a Gamma Batch and contract to OccJ
                   IF(dimAlpha2.NE.dimAlphaMPI)call lsquit('dimAlpha2.NE.dimAlphaMPI',-1)
!                   call mem_alloc(tmp4,dimGamma2*i8*nOccBatchDimJ)         
!                   IF(dimGamma2*i8*nOccBatchDimJ.GT.size(tmp4))CALL LSQUIT('TEST2 tmp4',-1)

                   !MemoryBookkeeping: tmp2,tmp3,tmp1,tmp4 = 
                   !+ nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI
                   !+ nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ 
                   !+ nvirt*nOccBatchDimI,Max(dimAlpha2)*Max(dimGamma2) + dimGamma2*nOccBatchDimJ
                   call BuildCgammaMPIMPMP2(nOccBatchDimJ,MAXnOccBatchDimJrank,ncolnodes,&
                        & OccIndexJrank,MAXnBlocksGamma,nBlocksGamma,offsetGamma,AOdimGamma,&
                        & AOstartGamma,AOendGamma,Mymolecule%Co%elm2,nb,nocc,dimGamma2,tmp4,&
                        & jnode,jnodeloop,offset)
                   !share alpha batches
                   !tmp3(nvirt,noccBI,dimAlpha,noccBJ)=tmp1(nvirt,noccBI,dimAlpha,dimGammaMPI)*C(dimGammaMPI,noccBJ)
                   M = nvirt*nOccBatchDimI*dimAlphaMPI !rows of Output Matrix
                   N = nOccBatchDimJ                   !columns of Output Matrix
                   K = dimGamma2                       !summation dimension
                   call dgemm('N','N',M,N,K,1.0E0_realk,tmp1,M,tmp4,K,1.0E0_realk,tmp3,M)
!                   call mem_dealloc(tmp4)
                ELSE
                   !recv something in BCAST that was not needed - could make new comm for subset.
                ENDIF
!                call mem_dealloc(tmp1)
                CALL LS_GETTIM(CPU4,WALL4)
                CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
                WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
             ENDIF
          enddo
       enddo
       CALL LS_GETTIM(CPU3,WALL3)
!       call mem_dealloc(tmp2)

       !reorder: tmp2(dimAlpha,noccBJ,nvirt,noccBI) <= tmp3(nvirt,noccBI,dimAlpha,noccBJ)
!       call mem_alloc(tmp2,dimAlphaMPI*nOccBatchDimJ*i8*nvirt*nOccBatchDimI)
!       IF(dimAlphaMPI*nOccBatchDimJ*i8*nvirt*nOccBatchDimI.GT.size(tmp2))CALL LSQUIT('TEST3 tmp2',-1)
       !MemoryBookkeeping: tmp3,tmp2 = 
       !+ nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ
       !+ dimAlphaMPI*nOccBatchDimJ,nvirt*nOccBatchDimI

       M = nvirt*nOccBatchDimI          !row of Input Matrix
       N = dimAlphaMPI*nOccBatchDimJ    !columns of Input Matrix
       call mat_transpose(M,N,1.0E0_realk,tmp3,0.0E0_realk,tmp2)
!       call mem_dealloc(tmp3)

!       call mem_alloc(tmp3,nvirt*nOccBatchDimJ*i8*nvirt*nOccBatchDimI)
       !MemoryBookkeeping: tmp2,tmp3 =        
       !+ dimAlphaMPI*nOccBatchDimJ,nvirt*nOccBatchDimI
       !+ nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI

       !Construct CAI(dimAlphaMPI,nvirt) 
!       call mem_alloc(tmp1,dimAlphaMPI,nvirt)     
       !MemoryBookkeeping: tmp2,tmp3,tmp1 =        
       !+ dimAlphaMPI*nOccBatchDimJ,nvirt*nOccBatchDimI
       !+ nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI
       !+ dimAlphaMPI*nvirt
       call BuildCAVMPMP2(nvirt,nBlocksAlpha,nrownodes,tmp1,Mymolecule%Cv%elm2,&
            & nb,dimAlphaMPI,inode,MAXnBlocksAlpha,offsetAlpha,AOdimAlpha,AOstartAlpha,AOendAlpha)
       !tmp3(nvirt,noccBJ,nvirt,noccBI) = tmp1(dimAlpha,nvirt)*tmp2(dimAlpha*nOccBatchDimJ,nvirt*nOccBatchDimI)
       M = nvirt                              !rows of Output Matrix
       N = nOccBatchDimJ*nvirt*nOccBatchDimI  !columns of Output Matrix
       K = dimAlphaMPI                        !summation dimension
       call dgemm('T','N',M,N,K,1.0E0_realk,tmp1,K,tmp2,K,0.0E0_realk,tmp3,M)     
!       call mem_dealloc(tmp1)
!       call mem_dealloc(tmp2)
       CALL LS_GETTIM(CPU4,WALL4)
       CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
       WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)

       !MemoryBookkeeping: tmp3 = nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI

       CALL LS_GETTIM(CPU1,WALL1)
       !at this point all have a tmp3(nvirt,noccBJ,nvirt,noccBI)
       !inode,jnode have the full tmp3(nvirt,noccBJ,nvirt,noccBI) 
       !with noccBJ determined by jnode but it only contains 
       !contributions from Alpha AO batch(inode)
       !so we now collect the Alpha AO batch(inode) to have a full tmp3(nvirt,noccBJ,nvirt,noccBI) 
       !on node inode=1 
       nbuf1 = nvirt*nOccBatchDimJ  
       nbuf2 = nvirt*nOccBatchDimI  
       !I share jnode and noccBJ with the rest (1:nrownodes,jnode)
       receiver = (1+(jnode-1)*nrownodes)-1   !(1,jnode)
       IF(inode.EQ.1)THEN !I collect results
          IF(receiver.NE.mynum) call lsquit('MPMP2 error in sendredv VOVO')
          do inodeLoop = 2,nrownodes !all send their contribution to inode=1,jnode=jnode
             sender = inodeLoop+(jnode-1)*nrownodes-1 !(inode,jnode), inode=2,nrownode
#ifdef VAR_MPI
!             call mem_alloc(tmp1,nvirt*nOccBatchDimJ*i8*nvirt*nOccBatchDimI)
             !MemoryBookkeeping: tmp3,tmp1 = nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI
             !+nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI
             call ls_mpisendrecv(tmp1,nbuf1*i8*nbuf2,comm,sender,receiver)
         
             call AddToVOVOMPMP2(tmp3,tmp1,nvirt,nOccBatchDimI,nOccBatchDimJ)
!             call mem_dealloc(tmp1)
#endif
          enddo
       ELSE
#ifdef VAR_MPI
          !I (inode,jnode) send to (1,jnode)
          call ls_mpisendrecv(tmp3,nbuf1*i8*nbuf2,comm,mynum,receiver)
#endif
       ENDIF
       CALL LS_GETTIM(CPU2,WALL2)
       CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
       WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)

       !MemoryBookkeeping: tmp3 = nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI

       CALL LS_GETTIM(CPU3,WALL3)
       IF(inode.EQ.1)THEN
!          call mem_alloc(tmp2,nvirt,nvirt)
!          call mem_alloc(Tmp1,nvirt,nvirt)
          !MemoryBookkeeping: tmp3,tmp2,tmp1 = nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI
          !+ nvirt*nvirt + nvirt*nvirt
          tmp_mp2_energy = 0.0E0_realk

          !tmp3 = VOVO(nvirt,nOccBatchDimJ,nvirt,nOccBatchDimI) 
          do I=1,nOccBatchDimI
             do J=1,nOccBatchDimJ
                JB = OccIndexJrank(J,JNODE)
                epsIJ = EpsOcc(I+(iB-1)*nOccBatchDimImax) + EpsOcc(JB)
                CALL CalcAmat2(nOccBatchDimJ,nOccBatchDimI,nvirt,tmp3,tmp2,J,I)
                call CalcBmat(nvirt,EpsIJ,EpsVirt,tmp2,Tmp1)
                tmp_mp2_energy2 = 0.0E0_realk
                call MP2_EnergyContribution(nvirt,tmp2,Tmp1,tmp_mp2_energy2)
                tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
             enddo
          enddo
          !        call lsmpi_reduction(tmp_mp2_energy,infpar%master,comm)
          receiver = 0 !master ?        
          IF(master)THEN 
             do jnodeLoop = 2,ncolnodes !for all jnodes, all noccBJ contributions
                tmp_mp2_energy2 = 0.0E0_realk
                sender = 1+(jnodeLoop-1)*nrownodes-1
#ifdef VAR_MPI
                call ls_mpisendrecv(tmp_mp2_energy2,comm,sender,receiver)
#endif
                WRITE(DECinfo%output,*)'MP2 Energy(iB=',iB,',mynum=',sender,') = ',tmp_mp2_energy
                tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
             enddo
          ELSE
#ifdef VAR_MPI
             call ls_mpisendrecv(tmp_mp2_energy,comm,mynum,receiver)
#endif
          ENDIF
!          call mem_dealloc(tmp2)
!          call mem_dealloc(Tmp1)
       ENDIF
!       call mem_dealloc(tmp3)
       !MemoryBookkeeping: 0

       CALL LS_GETTIM(CPU4,WALL4)
       CPU_ECONT = CPU_ECONT + (CPU4-CPU3)
       WALL_ECONT = WALL_ECONT + (WALL4-WALL3)           
       IF(master)THEN 
          WRITE(DECinfo%output,*)'MP2 Energy(iB=',iB,') = ',tmp_mp2_energy
          mp2_energy = mp2_energy + tmp_mp2_energy        
          !Write Restart File
          restart_lun = -1  !initialization
          call lsopen(restart_lun,'FULLMPMP2.restart','UNKNOWN','FORMATTED')
          rewind restart_lun
          write(restart_lun,'(I9)') nOccbatchesI
          write(restart_lun,'(I9)') iB
          write(restart_lun,'(F28.16)') mp2_energy
          call lsclose(restart_lun,'KEEP')
       ENDIF
    enddo BatchOccI
    IF(master)THEN 
       print*,'FINAL MP2 ENERGY',mp2_energy,'MYNUM',MYNUM
    ENDIF

    call mem_dealloc(tmp1)
    call mem_dealloc(tmp2)
    call mem_dealloc(tmp3)
    call mem_dealloc(tmp4)

    call mem_dealloc(nOccBatchDimJrank)
    call mem_dealloc(OccIndexJrank)
    call mem_dealloc(JobAlpha)
    call mem_dealloc(JobGamma)
    call mem_dealloc(nBlocksAlpha)
    call mem_dealloc(nBlocksGamma)

    call mem_dealloc(offsetGamma)
    call mem_dealloc(AOstartGamma)
    call mem_dealloc(AOendGamma)
    call mem_dealloc(AOdimGamma)
    call mem_dealloc(AOdimGammaMPI)

    call mem_dealloc(offsetAlpha)
    call mem_dealloc(AOstartAlpha)
    call mem_dealloc(AOendAlpha)
    call mem_dealloc(AOdimAlpha)
    call mem_dealloc(AOdimAlphaMPI)

    IF(DECinfo%useIchor)THEN
       call FREE_SCREEN_ICHORERI()
       call mem_dealloc(Alphabatchinfo)
       call mem_dealloc(Gammabatchinfo)
    ELSE
       nullify(mylsitem%setting%LST_GAB_LHS)
       nullify(mylsitem%setting%LST_GAB_RHS)
       call free_decscreen(DECSCREEN)
       ! Free gamma batch stuff
       call mem_dealloc(orb2batchAlpha)
       call mem_dealloc(batchdimAlpha)
       call mem_dealloc(batchsizeAlpha)
       call mem_dealloc(batchindexAlpha)
       orb2batchAlpha => null()
       batchdimAlpha => null()
       batchsizeAlpha => null()
       batchindexAlpha => null()
       do idxx=1,nbatchesAlpha
          call mem_dealloc(batch2orbAlpha(idxx)%orbindex)
          batch2orbAlpha(idxx)%orbindex => null()
       end do
       call mem_dealloc(batch2orbAlpha)
       batch2orbAlpha => null()

       call mem_dealloc(orb2batchGamma)
       call mem_dealloc(batchdimGamma)
       call mem_dealloc(batchsizeGamma)
       call mem_dealloc(batchindexGamma)
       orb2batchGamma => null()
       batchdimGamma => null()
       batchsizeGamma => null()
       batchindexGamma => null()
       do idxx=1,nbatchesGamma
          call mem_dealloc(batch2orbGamma(idxx)%orbindex)
          batch2orbGamma(idxx)%orbindex => null()
       end do
       call mem_dealloc(batch2orbGamma)
       batch2orbGamma => null()
    ENDIF
    call mem_dealloc(EpsOcc)
    call mem_dealloc(EpsVirt)
    !  call mem_dealloc(JobsCompleted)

    IF(MASTER)THEN
       write(lupri,*)  ''
       write(lupri,*)  'MP2 CORRELATION ENERGY = ', mp2_energy
       write(*,'(1X,a,f20.10)') 'MP2 CORRELATION ENERGY = ', mp2_energy
       write(lupri,*)  ''
    ENDIF
    CALL LSTIMER('FULL CANONICAL MPMP2',TS,TE,DECINFO%OUTPUT,FORCEPRINT)
    CALL ls_TIMTXT('>>>  WALL Time used in MP2 AO integral evaluation ',WALL_AOINT,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MP2 AO integral evaluation  ',CPU_AOINT,lupri)
    CALL ls_TIMTXT('>>>  WALL Time used in MP2 AO to MO transformation',WALL_AOTOMO,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MP2 AO to MO transformation ',CPU_AOTOMO,lupri)
    CALL ls_TIMTXT('>>>  WALL Time used in MP2 Energy evaluation      ',WALL_ECONT,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MP2 Energy evaluation       ',CPU_ECONT,lupri)
#ifdef VAR_MPI
    write(lupri,*)'Overall Time spent in MPI Communication and MPI Wait for rank=',infpar%mynum
    CALL ls_TIMTXT('>>>  WALL Time used MPI Communication inc. some Wait',WALL_MPICOMM,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used MPI Communication inc. some Wait',CPU_MPICOMM,lupri)
    CALL ls_TIMTXT('>>>  WALL Time used in MPI Wait',WALL_MPIWAIT,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MPI Wait',CPU_MPIWAIT,lupri)
#endif    
    write(lupri,*) ' '
  end subroutine full_canonical_mpmp2

  subroutine buildCoIMPMP2(CoI,nb,nOccBatchDimI,nocc,Co,offset,iB,nOccBatchDimImax)
    implicit none 
    integer,intent(in) :: nb,nOccBatchDimI,nocc,offset,iB,nOccBatchDimImax
    real(realk),intent(in) :: Co(nb,nocc)
    real(realk),intent(inout) :: CoI(nb,nOccBatchDimI)
    !
    integer :: I
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
    !$OMP SHARED(CoI,nb,nOccBatchDimI,nocc,Co,offset,iB,nOccBatchDimImax)
    do I=1,nOccBatchDimI
       CoI(:,I) = Co(:,offset+I+(iB-1)*nOccBatchDimImax) 
    enddo
    !$OMP END PARALLEL DO
  end subroutine buildCoIMPMP2

  subroutine BuildCgammaMPIMPMP2(nOccBatchDimJ,MAXnOccBatchDimJrank,ncolnodes,&
       & OccIndexJrank,MAXnBlocksGamma,&
       & nBlocksGamma,offsetGamma,AOdimGamma,AOstartGamma,AOendGamma,&
       & Co,nb,nocc,dimGammaMPI,tmp1,jnode,jnodeloop,offset)
    implicit none
    integer,intent(in) :: nOccBatchDimJ,MAXnOccBatchDimJrank,ncolnodes
    integer,intent(in) :: nb,nocc,dimGammaMPI,jnode,jnodeloop,offset
    integer,intent(in) :: MAXnBlocksGamma
    integer,intent(in) :: OccIndexJrank(MAXnOccBatchDimJrank,ncolnodes)
    integer,intent(in) :: nBlocksGamma(ncolnodes) 
    integer,intent(in) :: offsetGamma(MAXnBlocksGamma,ncolnodes)
    integer,intent(in) :: AOdimGamma(MAXnBlocksGamma,ncolnodes)
    integer,intent(in) :: AOstartGamma(MAXnBlocksGamma,ncolnodes)
    integer,intent(in) :: AOendGamma(MAXnBlocksGamma,ncolnodes)
    real(realk),intent(in) :: Co(nb,nocc)
    real(realk),intent(inout) :: tmp1(dimGammaMPI,nOccBatchDimJ)
    !
    integer :: J,JB,kk
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(J,JB,kk) SHARED(nOccBatchDimJ,MAXnOccBatchDimJrank,&
    !$OMP ncolnodes,OccIndexJrank,MAXnBlocksGamma,nBlocksGamma,offsetGamma,AOdimGamma,AOstartGamma,&
    !$OMP AOendGamma,Co,nb,nocc,dimGammaMPI,tmp1,jnode,jnodeloop,offset)
    do J=1,nOccBatchDimJ
       JB = OccIndexJrank(J,JNODE)
       do kk = 1,nBlocksGamma(jnodeLoop)
          tmp1(offsetGamma(kk,jnodeLoop)+1:offsetGamma(kk,jnodeLoop)+AOdimGamma(kk,jnodeLoop),J)=&
               &Co(AOstartGamma(kk,jnodeLoop):AOendGamma(kk,jnodeLoop),offset+JB)
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine BuildCgammaMPIMPMP2

  subroutine BuildCAVMPMP2(nvirt,nBlocksAlpha,nrownodes,tmp1,Cv,nb,dimAlphaMPI,&
       & inode,MAXnBlocksAlpha,offsetAlpha,AOdimAlpha,AOstartAlpha,AOendAlpha)
    implicit none
    integer,intent(in) :: nvirt,nrownodes,nb,dimAlphaMPI,inode,MAXnBlocksAlpha
    integer,intent(in) :: nBlocksAlpha(nrownodes)
    integer,intent(in) :: offsetAlpha(MAXnBlocksAlpha,nrownodes)
    integer,intent(in) :: AOdimAlpha(MAXnBlocksAlpha,nrownodes)
    integer,intent(in) :: AOstartAlpha(MAXnBlocksAlpha,nrownodes)
    integer,intent(in) :: AOendAlpha(MAXnBlocksAlpha,nrownodes)
    real(realk),intent(in) :: Cv(nb,nvirt)
    real(realk),intent(inout) :: tmp1(dimAlphaMPI,nvirt)
    !
    integer :: a,kk
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(a,kk) SHARED(nvirt,nBlocksAlpha,nrownodes,&
    !$OMP tmp1,Cv,nb,dimAlphaMPI,inode,MAXnBlocksAlpha,offsetAlpha,AOdimAlpha,AOstartAlpha,&
    !$OMP AOendAlpha)
    do a=1,nvirt
       do kk = 1,nBlocksAlpha(inode)
          tmp1(offsetAlpha(kk,inode)+1:offsetAlpha(kk,inode)+AOdimAlpha(kk,inode),A) = &
               Cv(AOstartAlpha(kk,inode):AOendAlpha(kk,inode),A) 
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine BuildCAVMPMP2

  subroutine AddIntegralToCollectMPMP2(tmp2,tmp3,nOccBatchDimI,nb,dimGamma,dimAlpha,&
       & dimAlphaMPI,dimGammaMPI,offsetA,offsetG)
    implicit none
    integer,intent(in) :: nOccBatchDimI,nb,dimGamma,dimAlpha,dimAlphaMPI,dimGammaMPI
    integer,intent(in) :: offsetA,offsetG
    real(realk),intent(in) :: tmp3(dimAlpha,dimGamma,nb,nOccBatchDimI)
    real(realk),intent(inout) :: tmp2(dimAlphaMPI,dimGammaMPI,nb,nOccBatchDimI)
    !
    integer :: j,b,i,a
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(I,J,A,&
    !$OMP B) SHARED(tmp2,tmp3,nOccBatchDimI,nb,dimGamma,dimAlpha,&
    !$OMP dimAlphaMPI,dimGammaMPI,offsetA,offsetG)
    do j=1,nOccBatchDimI
       do b=1,nb
          do i=1,dimGamma
             do a=1,dimAlpha
                tmp2(a+offsetA,i+offsetG,b,j) = tmp3(a,i,b,j)
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine AddIntegralToCollectMPMP2

  subroutine AddToVOVOMPMP2(VOVO,VOVO2,nvirt,nOccBatchDimI,nOccBatchDimJ)
    implicit none
    integer,intent(in) :: nvirt,nOccBatchDimI,nOccBatchDimJ
    real(realk),intent(in) :: VOVO2(nOccBatchDimJ*nvirt,nOccBatchDimI*nvirt)
    real(realk),intent(inout) :: VOVO(nOccBatchDimJ*nvirt,nOccBatchDimI*nvirt)
    !
    integer :: I,J               
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(nOccBatchDimI,nvirt,nOccBatchDimJ,VOVO,VOVO2)
    do J=1,nOccBatchDimI*nvirt
       do I=1,nOccBatchDimJ*nvirt
          VOVO(I,J) = VOVO(I,J) + VOVO2(I,J)
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine AddToVOVOMPMP2

  subroutine get_optimal_batch_sizes_for_canonical_mpmp2(MinAObatch,nbasis,nocc,nvirt,&
       & numnodes,nrownodes,ncolnodes,MaxAllowedDimAlpha,MaxAllowedDimGamma,&
       & MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,nOccBatchDimImax,nOccBatchDimJmax,&
       & sizetmp1,sizetmp2,sizetmp3,sizetmp4)
    implicit none
    integer,intent(in) :: MinAObatch,nbasis,nocc,nvirt,numnodes
    integer,intent(inout) :: MaxAllowedDimAlpha,MaxAllowedDimGamma
    integer,intent(inout) :: MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI
    integer,intent(inout) :: nOccBatchDimImax,nOccBatchDimJmax
    integer,intent(inout) :: nrownodes,ncolnodes
    integer(kind=8),intent(inout) :: sizetmp1,sizetmp2,sizetmp3,sizetmp4
    !local variables
    integer(kind=8) :: sizetmp1T,sizetmp2T,sizetmp3T,sizetmp4T
    real(realk) :: MemoryAvailable,GB,AG,maxsize
    integer :: iB,jB,iiB,jjB,nb,dimGamma,dimAlpha,AB,iGB,nTasks,K,tmprow,tmpcol,I,J
    integer :: nOccBatchDimImax2
    real(realk) :: nbasisR,noccR,nvirtR,numnodesR
    logical :: Success
    integer(kind=ls_mpik) :: nodtot 
    nbasisR = nbasis
    noccR = nocc
    nvirtR = nvirt
    numnodesR = numnodes   

    !The full_canonical_mpmp2 requires that (nb,nb,nb) can be distributed across all nodes 
    nodtot = INT(numnodes)
    call canonical_mpmp2_memreq_test(nbasis,nodtot,Success)

    call get_currently_available_memory(MemoryAvailable)
    ! Note: We multiply by 85 % to be on the safe side!
    MemoryAvailable = 0.85E0_realk*MemoryAvailable
    GB = 8.000E-9_realk 

    !assume you have 
    !numnodes = 144
    !nbasis = 3772
    !nvirt = 3508
    !nocc = 264
    !1. canonical_mpmp2_memreq_test already done
    !   nbasis*nbasis*nbasis = 430 GB can be distributed among 40 nodes (10.7 GB on each)  

    !2.Choose  nOccBatchDimImax as big as possible (nb*dimAlphaMPI*dimGammaMPI*nOccBatchDimImax) need to fit in mem!
    !          Same as (nb*nb*nb*nOccBatchDimImax/numnodes) 
    nOccBatchDimImax = MIN(nocc,FLOOR((MemoryAvailable*numnodes)/(nbasisR*nbasisR*nbasisR*GB))) 
    nOccBatchDimImax2 = nOccBatchDimImax
    do I=nOccBatchDimImax2,1,-1
       call memestimateCANONMPMP2(MinAObatch,MinAObatch,nbasis,nbasis,&
            & nbasis,nvirt,I,nOcc,maxsize,sizetmp1T,sizetmp2T,sizetmp3T,sizetmp4T)
       IF(maxsize.LT.MemoryAvailable)THEN
          nOccBatchDimImax = I
          sizetmp1 = sizetmp1T
          sizetmp2 = sizetmp2T
          sizetmp3 = sizetmp3T
          sizetmp4 = sizetmp4T
          EXIT
       ENDIF
       IF(I.EQ.1)THEN
!          print*,'Not enough Memory in MP2(canonical_mpmp2) '
!          print*,'(nbasis,nocc,nvirt,MinAObatch,numnodes):',nbasis,nocc,nvirt,MinAObatch,numnodes
!          print*,'With nOccBatchDimImax=1 the code will require'
!          print*,'maxsize = ',maxsize,' Which is less than MemoryAvailable=',MemoryAvailable
!          call lsquit('get_optimal_batch_sizes_for_canonical_mpmp2 Error MinAoBatches not enough',-1)
          nOccBatchDimImax = I
       ENDIF
    enddo

    !nOccBatchDimImax = 1 for this example 
    !This means recalculation of integrals 264 times for this example 

    !We divide the numnodes into an array of nodes (inode,jnode)
    !for numnodes=4 
    !mynum=0    means inode=1, jnode=1
    !mynum=1    means inode=2, jnode=1
    !mynum=2    means inode=1, jnode=2
    !mynum=3    means inode=2, jnode=2
    !Chose nrownodes and ncolnodes so that: 
    !numnodes = nrownodes*ncolnodes
    !and as square as possible (but VOVO should still fit in mem)
    !ncolnodes .GE. nrownodes
    ncolnodes = numnodes
    nrownodes = 1
    K=1
    do 
       K=K+1
       IF(numnodes+1.LE.K*K)EXIT
       tmprow = K
       tmpcol = numnodes/K
       IF(tmprow*tmpcol.EQ.numnodes)THEN
          IF(tmprow+tmpcol.LE.ncolnodes+ncolnodes)THEN
             !This is not correct - the AO batches of sizes (3,1,9,5) distributed among 2 nodes gives (10,8) not 9
             MaxAllowedDimAlphaMPI = CEILING(1.0E0_realk*nbasis/nrownodes)+MinAObatch 
             MaxAllowedDimGammaMPI = CEILING(1.0E0_realk*nbasis/ncolnodes)+MinAObatch 
             nOccBatchDimJmax = CEILING(1.0E0_realk*nocc/ncolnodes) 
             call memestimateCANONMPMP2(MinAObatch,MinAObatch,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,&
                  & nbasis,nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize,sizetmp1T,sizetmp2T,sizetmp3T,sizetmp4T)
             !           print*,'TRY tmprow,tmpcol',tmprow,tmpcol,'maxsize',maxsize,'MemoryAvailable',MemoryAvailable
             IF(maxsize.LT.MemoryAvailable)THEN
                nrownodes = tmprow
                ncolnodes = tmpcol
                sizetmp1 = sizetmp1T
                sizetmp2 = sizetmp2T
                sizetmp3 = sizetmp3T
                sizetmp4 = sizetmp4T
             ENDIF
          ENDIF
       ENDIF
    enddo
    !nOccBatchesI must match number of Alpha batches(1dim) = nrownodes
    MaxAllowedDimAlphaMPI = CEILING(1.0E0_realk*nbasis/nrownodes)+MinAObatch  !19/2 = 10
    MaxAllowedDimGammaMPI = CEILING(1.0E0_realk*nbasis/ncolnodes)+MinAObatch  !19/2 = 10
    nOccBatchDimJmax = CEILING(1.0E0_realk*nocc/ncolnodes) 
    !  print*,'MinAObatch',MinAObatch
    call memestimateCANONMPMP2(MinAObatch,MinAObatch,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,&
         & nbasis,nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize,sizetmp1T,sizetmp2T,sizetmp3T,sizetmp4T)
    IF(maxsize.GT.MemoryAvailable)THEN
       print*,'get_optimal_batch_sizes_for_canonical_mpmp2 Error MinAoBatches'
       print*,'(nbasis,nocc,nvirt,MinAObatch,numnodes):',nbasis,nocc,nvirt,MinAObatch,numnodes
       print*,'It was decided to use'
       print*,'nOccBatchDimImax=',nOccBatchDimImax
       print*,'MaxAllowedDimAlphaMPI=',MaxAllowedDimAlphaMPI
       print*,'MaxAllowedDimGammaMPI=',MaxAllowedDimGammaMPI       
       print*,'but even with  MaxAllowedDimGamma = MinAObatch = ',MinAObatch  
       print*,'and with MaxAllowedDimAlpha = MinAObatch = ',MinAObatch  
       print*,'the maxsize = ',maxsize,' is less than MemoryAvailable=',MemoryAvailable
       call lsquit('get_optimal_batch_sizes_for_canonical_mpmp2 Error MinAoBatches',-1)        
    ENDIF

    !assume 
    MaxAllowedDimGamma = MinAObatch  
    !find MaxAllowedDimAlpha as big as possible
    DO I=MaxAllowedDimAlphaMPI,MinAObatch,-1  
       call memestimateCANONMPMP2(I,MinAObatch,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,nbasis,&
            & nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize,sizetmp1T,sizetmp2T,sizetmp3T,sizetmp4T)
       MaxAllowedDimAlpha = I
       sizetmp1 = sizetmp1T
       sizetmp2 = sizetmp2T
       sizetmp3 = sizetmp3T
       sizetmp4 = sizetmp4T
       IF(maxsize.LT.MemoryAvailable)EXIT
       IF(I.EQ.MinAObatch)THEN
          print*,'get_optimal_batch_sizes_for_canonical_mpmp2 Error MinAoBatchAlpha'
          print*,'(nbasis,nocc,nvirt,MinAObatch,numnodes):',nbasis,nocc,nvirt,MinAObatch,numnodes
          print*,'It was decided to use'
          print*,'nOccBatchDimImax=',nOccBatchDimImax
          print*,'MaxAllowedDimAlphaMPI=',MaxAllowedDimAlphaMPI
          print*,'MaxAllowedDimGammaMPI=',MaxAllowedDimGammaMPI       
          print*,'but even with  MaxAllowedDimGamma = MinAObatch = ',MinAObatch  
          print*,'and with MaxAllowedDimAlpha = MinAObatch = ',MinAObatch  
          print*,'the maxsize = ',maxsize,' is less than MemoryAvailable=',MemoryAvailable
          call lsquit('get_optimal_batch_sizes_for_canonical_mpmp2 Error MinAoBatchAlpha',-1)        
       ENDIF
    ENDDO

    DO I=MaxAllowedDimGammaMPI,MinAObatch,-1
       call memestimateCANONMPMP2(MaxAllowedDimAlpha,I,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,nbasis,&
            & nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize,sizetmp1T,sizetmp2T,sizetmp3T,sizetmp4T)
       MaxAllowedDimGamma = I
       sizetmp1 = sizetmp1T
       sizetmp2 = sizetmp2T
       sizetmp3 = sizetmp3T
       sizetmp4 = sizetmp4T
       IF(maxsize.LT.MemoryAvailable)EXIT
       IF(I.EQ.MinAObatch)THEN
          print*,'get_optimal_batch_sizes_for_canonical_mpmp2 Error MinAoBatchGamma'
          print*,'(nbasis,nocc,nvirt,MinAObatch,numnodes):',nbasis,nocc,nvirt,MinAObatch,numnodes
          print*,'It was decided to use'
          print*,'nOccBatchDimImax=',nOccBatchDimImax
          print*,'MaxAllowedDimAlphaMPI=',MaxAllowedDimAlphaMPI
          print*,'MaxAllowedDimGammaMPI=',MaxAllowedDimGammaMPI       
          print*,'MaxAllowedDimAlpha   =',MaxAllowedDimAlpha
          print*,'and with MaxAllowedDimGamma = MinAObatch = ',MinAObatch  
          print*,'the maxsize = ',maxsize,' is less than MemoryAvailable=',MemoryAvailable
          call lsquit('get_optimal_batch_sizes_for_canonical_mpmp2 Error MinAoBatchGamma',-1)        
       ENDIF
    ENDDO

  end subroutine get_optimal_batch_sizes_for_canonical_mpmp2

  subroutine memestimateCANONMPMP2(MaxAllowedDimAlpha,MaxAllowedDimGamma,dimAlphaMPI,dimGammaMPI,nb,&
       & nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize,sizetmp1,sizetmp2,sizetmp3,sizetmp4)
    implicit none
    integer,intent(in) :: MaxAllowedDimAlpha,MaxAllowedDimGamma,dimAlphaMPI,dimGammaMPI,nb
    integer,intent(in) :: nvirt,nOccBatchDimImax,nOccBatchDimJmax
    real(realk),intent(inout) :: maxsize
    integer(kind=8),intent(inout) :: sizetmp1,sizetmp2,sizetmp3,sizetmp4
    integer(kind=8) :: nOccBatchDimI,nOccBI,nOccBatchDimJ,dimAlpha,dimGamma
    real(realk) :: GB,R
    R=1.0E0_realk
    nOccBatchDimI = nOccBatchDimImax
    nOccBI = nOccBatchDimI
    nOccBatchDimJ = nOccBatchDimJmax
    dimAlpha=MaxAllowedDimAlpha
    dimGamma=MaxAllowedDimGamma
    !construct CoI(nb,nOccBatchDimI)
    sizetmp1 = 0_8
    sizetmp2 = 0_8
    sizetmp3 = 0_8

    maxsize = nb*nOccBI
    sizetmp4 = nb*nOccBI*i8
    !call mem_alloc(tmp2,dimAlphaMPI,dimGammaMPI,nb,nOccBatchDimI)
    maxsize = MAX(maxsize,nb*nOccBI*R+dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI*R)
    sizetmp2 = MAX(sizetmp2,dimAlphaMPI*i8*dimGammaMPI*nb*i8*nOccBatchDimI)

    !BatchGamma: do gammaB = 1,nbatchesGamma
    ! BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
    !  call mem_alloc(tmp1,dimAlpha*dimGamma,nb*nb)
    maxsize = MAX(maxsize,nb*nOccBI*R+dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI*R+dimAlpha*dimGamma*nb*nb*R)
    sizetmp1 = MAX(sizetmp1,dimAlpha*i8*dimGamma*nb*i8*nb)

    !  call mem_alloc(tmp3,dimAlpha,dimGamma,nb,nOccBatchDimI)
    maxsize = MAX(maxsize,nb*nOccBI*R+dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI*R+dimAlpha*dimGamma*nb*nb*R+&
         & dimAlpha*dimGamma*nb*nOccBatchDimI*R)
    sizetmp3 = MAX(sizetmp3,dimAlpha*i8*dimGamma*nb*i8*nOccBatchDimI)

    !  call mem_dealloc(tmp1)
    !  call mem_dealloc(tmp3)
    ! enddo BatchAlpha
    !enddo BatchGamma
    !call mem_dealloc(CoI)
    !call mem_alloc(tmp3,nb*nOccBatchDimI,dimAlphaMPI*dimGammaMPI)
    maxsize = MAX(maxsize,dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI*R+nb*nOccBatchDimI*dimAlphaMPI*dimGammaMPI*R) 
    sizetmp3 = MAX(sizetmp3,nb*nOccBatchDimI*i8*dimAlphaMPI*i8*dimGammaMPI)

    !call mem_dealloc(tmp2)
    !call mem_alloc(tmp2,nvirt*nOccBatchDimI,dimAlphaMPI*dimGammaMPI)
    maxsize = MAX(maxsize,nb*nOccBatchDimI*dimAlphaMPI*dimGammaMPI*R+nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI*R) 
    sizetmp2 = MAX(sizetmp2,nvirt*nOccBatchDimI*i8*dimAlphaMPI*i8*dimGammaMPI)

    !call mem_dealloc(tmp3)
    !call mem_alloc(tmp3,nvirt*nOccBatchDimI,dimAlphaMPI*nOccBatchDimJ)
    maxsize = MAX(maxsize,nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI*R+nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ*R) 
    sizetmp3 = MAX(sizetmp3,nvirt*nOccBatchDimI*i8*dimAlphaMPI*i8*nOccBatchDimJ)

    !OPTION1
    !call mem_alloc(tmp1,dimGammaMPI,nOccBatchDimJ)
    maxsize = MAX(maxsize,nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI*R+nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ*R+& 
         & dimGammaMPI*nOccBatchDimJ*R)
    sizetmp1 = MAX(sizetmp1,dimGammaMPI*i8*nOccBatchDimJ)

    !call mem_dealloc(tmp1)
    !OPTION1
    !call mem_alloc(tmp1,nvirt*nOccBatchDimI,dimAlpha2*dimGamma2)              
    sizetmp1 = MAX(sizetmp1,nvirt*nOccBatchDimI*i8*dimAlphaMPI*dimGammaMPI)
    !call mem_alloc(tmp4,dimGamma2,nOccBatchDimJ)
    maxsize = MAX(maxsize,nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI*R+nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ*R+& !tmp5+tmp4
         & nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI*R + dimGammaMPI*nOccBatchDimJ*R) !tmp6
    sizetmp4 = MAX(sizetmp4,dimGammaMPI*i8*nOccBatchDimJ)

    !call mem_dealloc(tmp4)
    !call mem_dealloc(tmp1)
    !call mem_dealloc(tmp2)
    !call mem_alloc(tmp2,dimAlphaMPI*nOccBatchDimJ,nvirt*nOccBatchDimI)
    maxsize = MAX(maxsize,nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ*R & 
         & + dimAlphaMPI*nOccBatchDimJ*nvirt*nOccBatchDimI*R)
    sizetmp2 = MAX(sizetmp2,dimAlphaMPI*i8*nOccBatchDimJ*nvirt*i8*nOccBatchDimI)


    !call mem_dealloc(tmp3)
    !call mem_alloc(tmp3,nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
    maxsize = MAX(maxsize,dimAlphaMPI*nOccBatchDimJ*nvirt*nOccBatchDimI*R+nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI*R) 
    sizetmp3 = MAX(sizetmp3,nvirt*nOccBatchDimJ*i8*nvirt*i8*nOccBatchDimI)


    !call mem_alloc(tmp1,dimAlphaMPI,nvirt)     
    maxsize = MAX(maxsize,dimAlphaMPI*nOccBatchDimJ*nvirt*nOccBatchDimI*R+nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI*R+& 
         & dimAlphaMPI*nvirt*R)
    sizetmp1 = MAX(sizetmp1,dimAlphaMPI*i8*nvirt)


    !call mem_alloc(tmp1,nvirt*nOccBatchDimJ*i8*nvirt*nOccBatchDimI)
    maxsize = MAX(maxsize,nvirt*nOccBatchDimJ*i8*nvirt*nOccBatchDimI*R+nvirt*nOccBatchDimJ*i8*nvirt*nOccBatchDimI*R)
    sizetmp1 = MAX(sizetmp1,nvirt*nOccBatchDimJ*i8*nvirt*i8*nOccBatchDimI)


    !call mem_alloc(tmp1,nvirt*nOccBatchDimJ*i8*nvirt*nOccBatchDimI)
    maxsize = MAX(maxsize,nvirt*nOccBatchDimJ*i8*nvirt*nOccBatchDimI*R+nvirt*nvirt*R+nvirt*nvirt*R)
    sizetmp1 = MAX(sizetmp1,nvirt*nvirt*i8)
    sizetmp2 = MAX(sizetmp2,nvirt*nvirt*i8)

    GB = 8.000E-9_realk 
    maxsize = maxsize*GB
    
    WRITE(DECinfo%output,*)'Using reallocation gives maxsize=',maxsize,' GB'

    maxsize = (sizetmp1+sizetmp2+sizetmp3+sizetmp4)*GB
    WRITE(DECinfo%output,*)'Using 1 allocation at begining gives maxsize=',maxsize,' GB'

    !VOVO(nvirt,noccBJ,nvirt,noccBI) = CAV(dimAlpha,nvirt)*tmp7(dimAlpha*nOccBatchDimJ,nvirt*nOccBatchDimI)
    !call mem_dealloc(CAV)
    !call mem_dealloc(tmp7)
  end subroutine memestimateCANONMPMP2

  !The full_canonical_mpmp2 requires that (nb,nb,nb) can be distributed across all nodes 
  subroutine canonical_mpmp2_memreq_test(nbasis,numnodes,Success)
    implicit none
    integer(kind=ls_mpik),intent(in) :: numnodes
    integer,intent(in) :: nbasis
    logical,intent(inout) :: Success
    !
    real(realk) :: MemoryAvailable,GB,nbasisR,numnodesR
    ! Memory currently available
    ! **************************
    call get_currently_available_memory(MemoryAvailable)
    ! Note: We multiply by 85 % to be on the safe side!
    MemoryAvailable = 0.85E0_realk*MemoryAvailable
    GB = 1.000E-9_realk 
    nbasisR = nbasis
    numnodesR = numnodes
    IF(nbasisR*nbasisR*nbasisR*GB/numnodesR.GT.MemoryAvailable)THEN
       print*,'canonical_mpmp2_memreq_test  size(nb*nb*nb/numnodes) =',nbasisR*nbasisR*(nbasisR*GB)/numnodesR,' GB'
       print*,'canonical_mpmp2_memreq_test  MemoryAvailable         =',MemoryAvailable,' GB'
       call lsquit('canonical_mpmp2_memreq_test failure',-1)
    ENDIF
    Success = .TRUE.
  end subroutine canonical_mpmp2_memreq_test

  !> \brief Memory check for full_canonical_mp2 subroutine
  !> \author Thomas Kjaergaard
  !> \date October 2014
  subroutine full_canonical_mp2_memory_check(nbasis,nvirt,MinAtomic)
    implicit none
    integer,intent(in) :: nbasis,nvirt,MinAtomic
    real(realk) :: MemRequired,GB,MemStep1,MemStep2
    GB = 1.0E-9_realk
    ! Check that arrays fit in memory (hard-coded)
    MemRequired = real(MinAtomic*MinAtomic*nbasis*(nbasis+nvirt))
    MemRequired = MemRequired*realk*GB  
    if(MemRequired > 0.80E0_realk*DECinfo%memory) then
       call lsquit('full_canonical_mp2: Memory exceeded! ',-1)
    end if
  end subroutine full_canonical_mp2_memory_check

  !> \brief Calculate canonical MP2 energy for full molecular system
  !> \author Thomas Kjaergaard
  !> \date October 2014
  subroutine full_canonical_mp2(MyMolecule,MyLsitem,mp2_energy)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Canonical MP2 correlation energy
    real(realk),intent(inout) :: mp2_energy    
    !
    integer :: nbasis,nocc,nvirt,naux,noccfull,mynum,numnodes
    logical :: master,wakeslaves
    real(realk),pointer :: EpsOcc(:),EpsVirt(:)
    integer :: J,I,node,natoms,A,lupri,restart_lun
    logical :: MessageRecieved,MessageRecievedW,FORCEPRINT,file_exists
    real(realk) :: TS,TE,TS2,TE2,TS3,TE3,CPU1,CPU2,WALL1,WALL2,CPU_MPICOMM,WALL_MPICOMM
    real(realk) :: CPU_MPIWAIT,WALL_MPIWAIT,MemInGBCollected,epsIJ,tmp_mp2_energy
    real(realk) :: CPU3,CPU4,WALL3,WALL4,CPU_AOINT,WALL_AOINT,tmp_mp2_energy2
    real(realk) :: CPU_AOTOMO,WALL_AOTOMO,CPU_ECONT,WALL_ECONT
    real(realk),pointer :: Amat(:,:),Bmat(:,:),tmp1(:,:),tmp2(:,:),tmp3(:,:),tmp4(:,:)
    real(realk),pointer :: tmp5(:,:),tmp6(:,:),tmp7(:,:),CoBatchA(:,:),CoBatchB(:,:)
    real(realk),pointer :: CoBI(:,:),CoBJ(:,:),tmp62(:,:),tmp72(:,:),CoI(:,:),CoJ(:,:)
    real(realk),pointer :: VOVO(:,:),VGVO(:,:),CvA(:,:),CoIG(:,:)
    !
    type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
    type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:)
    integer, pointer :: batchdimAlpha(:), batchdimGamma(:)
    integer :: MaxAllowedDimAlpha,MaxAllowedDimGamma,MaxActualDimAlpha,MaxActualDimGamma
    integer :: nOccBatchDimImax,nOccBatchDimJmax,K,iorb,idx,nbatchesAlpha,nbatchesGamma,nb
    integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint,M,N,iB
    integer :: MinAObatch,gammaB,alphaB,nOccBatchesI,nOccBatchesJ,jB,nOccBatchDimI,nOccBatchDimJ
    integer :: nbatchesGammaRestart,nbatchesAlphaRestart,dimGamma,GammaStart,GammaEnd,dimAlpha
    integer :: AlphaStart,AlphaEnd,B,noccRestartI,noccRestartJ,nJobs,startjB,offset
    integer :: nOccBatchesIrestart,noccIstart,nbuf1,Ibuf(4)
    logical :: MoTrans, NoSymmetry,SameMol,JobDone,JobInfo1Free,FullRHS,doscreen,NotAllMessagesRecieved
    logical :: PermutationalSymmetryIJ
    logical,pointer :: JobsCompleted(:,:),JobsCompletedOrig(:,:)
    integer(kind=8) :: maxsize
    TYPE(DECscreenITEM)   :: DecScreen
    Character            :: intSpec(5)
    integer(kind=ls_mpik)  :: nMPI,TAG,IERR,request,Receiver,sender,comm,TAG1,TAG2
    integer(kind=ls_mpik)  :: request1,request2,masterrank,senderID
#ifdef VAR_MPI
    integer(kind=ls_mpik)  :: mpistatus(MPI_STATUS_SIZE) 
#endif
    integer(kind=4) :: JobInfo1(2)
    !  Character(80)        :: FilenameCS,FilenamePS

    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(full_canonical_mp2): this does not work with PDM&
       & molecular structure",-1)
    endif

#ifdef VAR_TIME    
    FORCEPRINT = .TRUE.
#endif

    CALL LSTIMER('START ',TS,TE,DECINFO%OUTPUT)
    mp2_energy = 0.0E0_realk

    CPU_AOINT = 0.0E0_realk
    WALL_AOINT = 0.0E0_realk
    CPU_AOTOMO = 0.0E0_realk
    WALL_AOTOMO = 0.0E0_realk
    CPU_ECONT = 0.0E0_realk
    WALL_ECONT = 0.0E0_realk

    CPU_MPICOMM = 0.0E0_realk
    WALL_MPICOMM = 0.0E0_realk
    CPU_MPIWAIT = 0.0E0_realk
    WALL_MPIWAIT = 0.0E0_realk
    TAG = 1411; TAG1 = 1412; TAG2 = 1413

    !sanity check
    if(.NOT.DECinfo%use_canonical) then
       call lsquit('Error: full_canonical_mp2 require canonical Orbitals',-1)
    endif
    ! Init stuff
    ! **********
    nbasis = MyMolecule%nbasis
    nb = nbasis
    !MyMolecule%Co is allocated (nbasis,MyMolecule%nocc)
    !with MyMolecule%nocc = Valence + Core 
    !In case of Frozen core we only need Valence and will access
    !Co(nbasis,nval) = MyMolecule%Co(nbasis,MyMolecule%ncore+1:MyMolecule%nocc)
    noccfull = MyMolecule%nocc
    IF(DECinfo%Frozencore)THEN
       nocc   = MyMolecule%nval
       offset = MyMolecule%ncore
    ELSE
       nocc   = MyMolecule%nocc
       offset = 0
    ENDIF
    nvirt  = MyMolecule%nvirt
    nAtoms = MyMolecule%nAtoms
    LUPRI = DECinfo%output

#ifdef VAR_MPI
    comm = MPI_COMM_LSDALTON
    master= (infpar%mynum == infpar%master)
    mynum = infpar%mynum
    numnodes = infpar%nodtot
    wakeslaves = infpar%nodtot.GT.1
    IF(.NOT.master)LUPRI = 6 !standard Output
#else
    ! If MPI is not used, consider the single node to be "master"
    master=.true.
    mynum = 0
    numnodes = 1
    wakeslaves = .false.
#endif

    ! Set integral info
    ! *****************
    !R = Regular Basis set on centers 1-4, C = Coulomb operator
    INTSPEC = ['R','R','R','R','C']

    !determine MinAObatch: the minimum allowed AObatch size + number of AO batches
    IF(DECinfo%useIchor)THEN
       !Determine the minimum allowed AObatch size MinAObatch
       !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
       !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch=6 (the 2*(Px,Py,Pz))
       !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
       !'R'  !Specifies that it is the Regular AO basis that should be batched
       iAO = 4 !the center that the batching should occur on (they are all the same in this case)  
       call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
       iprint = 0           !print level for Ichor Integral code
       MoTrans = .FALSE.    !Do not transform to MO basis! 
       NoSymmetry = .FALSE. !Use Permutational Symmetry! 
       SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
       !Determine the full number of AO batches - not to be confused with the batches of AOs
       !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
       iAO = 1
       call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
       nullify(AOGammabatchinfo)
       nullify(AOalphabatchinfo)    
    ELSE
       ! The smallest possible AO batch depends on the basis set
       ! (More precisely, if all batches are made as small as possible, then the
       !  call below determines the largest of these small batches).
       call determine_maxBatchOrbitalsize(DECinfo%output,mylsitem%setting,MinAObatch,'R')
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
    ENDIF

    doscreen = mylsitem%setting%scheme%cs_screen.OR.&
         & mylsitem%setting%scheme%ps_screen

#ifdef VAR_MPI 
    ! Master starts up slave
    StartUpSlaves: if(wakeslaves .and. master) then
       ! Wake up slaves to do the job: slaves awoken up with (RIMP2FULL)
       ! and call full_canonical_rimp2_slave which communicate info 
       ! then calls full_canonical_rimp2.
       CALL LS_GETTIM(CPU1,WALL1)
       call ls_mpibcast(CANONMP2FULL,infpar%master,comm)
       ! Communicate fragment information to slaves
       call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
       call mpicopy_lsitem(MyLsitem,comm)
       call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
       call mpi_bcast_fullmolecule(MyMolecule)    
       CALL LS_GETTIM(CPU2,WALL2)
       CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
       WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
    endif StartUpSlaves
#endif

    ! ***************************************************************************************
    !Get optimal values of: MaxAllowedDimAlpha,MaxAllowedDimGamma,nOccBatchDimImax,nOccBatchDimJmax
    ! ***************************************************************************************

#ifdef VAR_MPI
    !use the numbers obtained by master  
    IF(master)THEN
       call get_optimal_batch_sizes_for_canonical_mp2(MinAObatch,nbasis,nocc,nvirt,&
            & MaxAllowedDimAlpha,MaxAllowedDimGamma,nOccBatchDimImax,nOccBatchDimJmax,&
            & numnodes,DECinfo%output)
       
       Ibuf(1) = nOccBatchDimImax
       Ibuf(2) = nOccBatchDimJmax
       Ibuf(3) = MaxAllowedDimAlpha
       Ibuf(4) = MaxAllowedDimGamma
    ELSE
       Ibuf = 0
    ENDIF
    nbuf1 = 4
    call ls_mpibcast(Ibuf,nbuf1,infpar%master,comm)
    IF(.NOT.master)THEN
       nOccBatchDimImax = Ibuf(1)
       nOccBatchDimJmax = Ibuf(2)
       MaxAllowedDimAlpha = Ibuf(3)
       MaxAllowedDimGamma = Ibuf(4)
    ENDIF
#else
    call get_optimal_batch_sizes_for_canonical_mp2(MinAObatch,nbasis,nocc,nvirt,&
         & MaxAllowedDimAlpha,MaxAllowedDimGamma,nOccBatchDimImax,nOccBatchDimJmax,&
         & numnodes,DECinfo%output)
#endif

    write(DECinfo%output,*)'nbasis',nbasis
    write(DECinfo%output,*)'nocc  ',nocc
    write(DECinfo%output,*)'nvirt ',nvirt
    write(DECinfo%output,*)'nOccBatchDimImax  ',nOccBatchDimImax
    write(DECinfo%output,*)'nOccBatchDimJmax  ',nOccBatchDimJmax
    write(DECinfo%output,*)'MinAObatch        ',MinAObatch
    write(DECinfo%output,*)'MaxAllowedDimAlpha',MaxAllowedDimAlpha
    write(DECinfo%output,*)'MaxAllowedDimGamma',MaxAllowedDimGamma

    nOccBatchesI = nOcc/nOccBatchDimImax
    IF(MOD(nOcc,nOccBatchDimImax).NE.0)nOccBatchesI = nOccBatchesI + 1  
    nOccBatchesJ = nOcc/nOccBatchDimJmax
    IF(MOD(nOcc,nOccBatchDimJmax).NE.0)nOccBatchesJ = nOccBatchesJ + 1

    PermutationalSymmetryIJ = .FALSE.
    IF(nOccbatchesI.EQ.nOccbatchesJ)THEN
       PermutationalSymmetryIJ = .TRUE.
       WRITE(DECinfo%output,*)'Permutational Symmetry exploited'
    ENDIF
    write(DECinfo%output,*)'nOccBatchesI',nOccBatchesI
    write(DECinfo%output,*)'nOccBatchesJ',nOccBatchesJ

    ! ************************************************
    ! * Restart Option                               *
    ! ************************************************
    call mem_alloc(JobsCompleted,nOccBatchesI,nOccBatchesJ)
    call mem_alloc(JobsCompletedOrig,nOccBatchesI,nOccBatchesJ)
    JobDone = .FALSE.           
    IF(master)THEN
       IF(DECinfo%DECrestart)THEN
          !CHECK IF THERE ARE ENERGY CONTRIBUTIONS AVAILABLE
          INQUIRE(FILE='FULLMP2.restart',EXIST=file_exists)
          IF(file_exists)THEN
             WRITE(DECinfo%output,*)'Restart of Full canonical molecular MP2 calculation:'
             restart_lun = -1  !initialization
             call lsopen(restart_lun,'FULLMP2.restart','OLD','UNFORMATTED')
             rewind restart_lun
             read(restart_lun) noccRestartI,noccRestartJ
             IF(noccRestartI.NE.nOccBatchesI)THEN
                print*,'noccRestartI,nOccBatchesI',noccRestartI,nOccBatchesI
                call lsquit('nOccBatchesI must be same in Restart for MP2',-1)
             ENDIF
             IF(noccRestartJ.NE.nOccBatchesJ)THEN
                print*,'noccRestartJ,nOccBatchesJ',noccRestartJ,nOccBatchesJ
                call lsquit('nOccBatchesJ must be same in Restart for MP2',-1)
             ENDIF
             read(restart_lun) JobsCompleted
             IF(COUNT(JobsCompleted).EQ.nOccBatchesI*nOccBatchesJ)THEN
                WRITE(DECinfo%output,*)'All MP2 energies is on file JobsCompleted=',JobsCompleted
                JobDone = .TRUE.
             ELSE
                WRITE(DECinfo%output,*)'Restarting from file '
                WRITE(DECinfo%output,*) COUNT(JobsCompleted),' jobs completed out of ',nOccBatchesI*nOccBatchesJ
                JobDone = .FALSE.           
             ENDIF
             read(restart_lun) mp2_energy
             WRITE(DECinfo%output,*)'MP2 Energy Read From File: ',mp2_energy
             call lsclose(restart_lun,'KEEP')
          ELSE
             JobsCompleted = .FALSE.
          ENDIF
       ELSE
          JobDone = .FALSE.           
          JobsCompleted = .FALSE.
       ENDIF
    ELSE
       JobsCompleted = .FALSE.
       JobDone = .FALSE.
    ENDIF
    JobsCompletedOrig = JobsCompleted

#ifdef VAR_MPI
    IF(DECinfo%DECrestart)THEN
       CALL LS_GETTIM(CPU1,WALL1)
       call ls_mpibcast(JobsCompleted,nOccBatchesI,nOccBatchesJ,infpar%master,comm)
       CALL LS_GETTIM(CPU2,WALL2)
       CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
       WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
       IF(COUNT(JobsCompleted).EQ.nOccBatchesI*nOccBatchesJ)JobDone = .TRUE.
    ENDIF
#endif

    IF(JobDone)THEN
       !do nothing
    ELSE
       call mem_alloc(EpsOcc,nocc)
       !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
       !$OMP SHARED(noccfull,MyMolecule,EpsOcc,offset)
       do I=1+offset,noccfull
          EpsOcc(I-offset) = MyMolecule%oofock%elm2(I,I)
       enddo
       !$OMP END PARALLEL DO
       call mem_alloc(EpsVirt,nvirt)
       !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
       !$OMP SHARED(nvirt,MyMolecule,EpsVirt)
       do A=1,nvirt
          EpsVirt(A) = MyMolecule%vvfock%elm2(A,A)
       enddo
       !$OMP END PARALLEL DO

       ! ************************************************
       ! * Determine batch information for Gamma batch  *
       ! * And 
       ! * Determine batch information for Alpha batch  *
       ! ************************************************

       IF(DECinfo%useIchor)THEN
          iAO = 2 !Gamma is the 2. Center of the 4 center two electron coulomb integral
          !Determine how many batches of AOS based on the bat%MaxAllowedDimGamma, the requested
          !size of the AO batches. iAO is the center that the batching should occur on. 
          !'R'  !Specifies that it is the Regular AO basis that should be batched 
          call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
               & nbatchesGamma,DECinfo%output)
          call mem_alloc(AOGammabatchinfo,nbatchesGamma)
          !Construct the batches of AOS based on the bat%MaxAllowedDimGamma, the requested
          !size of the AO batches - bat%MaxAllowedDimGamma must be unchanged since the call 
          !to determine_Ichor_nbatchesofAOS
          !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
          !So MaxActualDimGamma must be less og equal to bat%MaxAllowedDimGamma
          call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
               & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)

          iAO = 1 !Alpha is the 1. Center of the 4 center two electron coulomb integral
          !Determine how many batches of AOS based on the bat%MaxAllowedDimAlpha, the requested
          !size of the AO batches. iAO is the center that the batching should occur on. 
          !'R'  !Specifies that it is the Regular AO basis that should be batched 
          call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
               & nbatchesAlpha,DECinfo%output)
          call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
          !Construct the batches of AOS based on the bat%MaxAllowedDimAlpha, the requested
          !size of the AO batches - bat%MaxAllowedDimAlpha must be unchanged since the call 
          !to determine_Ichor_nbatchesofAOS
          !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
          !So MaxActualDimAlpha must be less og equal to bat%MaxAllowedDimAlpha
          call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
               & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
       ELSE
          ! Orbital to batch information
          ! ----------------------------
          call mem_alloc(orb2batchGamma,nbasis)
          call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
               & nbasis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,&
               & orb2BatchGamma,'R')
          ! Translate batchindex to orbital index
          ! -------------------------------------
          call mem_alloc(batch2orbGamma,nbatchesGamma)
          do idx=1,nbatchesGamma
             call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
             batch2orbGamma(idx)%orbindex = 0
             batch2orbGamma(idx)%norbindex = 0
          end do
          do iorb=1,nbasis
             idx = orb2batchGamma(iorb)
             batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
             K = batch2orbGamma(idx)%norbindex
             batch2orbGamma(idx)%orbindex(K) = iorb
          end do
          ! Orbital to batch information
          ! ----------------------------
          call mem_alloc(orb2batchAlpha,nbasis)
          call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
               & nbasis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,&
               & orb2BatchAlpha,'R')
          ! Translate batchindex to orbital index
          ! -------------------------------------
          call mem_alloc(batch2orbAlpha,nbatchesAlpha)
          do idx=1,nbatchesAlpha
             call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
             batch2orbAlpha(idx)%orbindex = 0
             batch2orbAlpha(idx)%norbindex = 0
          end do
          do iorb=1,nbasis
             idx = orb2batchAlpha(iorb)
             batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
             K = batch2orbAlpha(idx)%norbindex
             batch2orbAlpha(idx)%orbindex(K) = iorb
          end do
       ENDIF

       CALL LS_GETTIM(CPU4,WALL4)

       ! ************************************************
       ! * Screening                                    *
       ! ************************************************
       IF(DECinfo%useIchor)THEN
          !Calculate Screening integrals 
          SameMOL = .TRUE. !Specifies same molecule on all centers 
          call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
       ELSE
          ! This subroutine builds the full screening matrix.
          call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
               & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%IntegralThreshold)
          IF(doscreen)THEN
             call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
                  & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
                  & batchindexAlpha,batchindexGamma,&
                  & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
          ENDIF
       ENDIF

       CALL LS_GETTIM(CPU3,WALL3)
       CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
       WALL_AOINT = WALL_AOINT + (WALL3-WALL4)

       ! ************************************************
       ! * Main Loop                                    *
       ! ************************************************
       JobInfo1Free = .FALSE.
       FullRHS = (nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)

       nJobs = 0 
       BatchOccI: do iB = 1,nOccbatchesI
          nOccBatchDimI = nOccBatchDimImax
          IF(MOD(nOcc,nOccBatchDimI).NE.0.AND.iB.EQ.nOccBatchesI)THEN
             !the remainder
             nOccBatchDimI = MOD(nOcc,nOccBatchDimImax)
          ENDIF

          !construct CoI(nb,nOccBatchDimI)
          call mem_alloc(CoI,nb,nOccBatchDimI)       
          !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
          !$OMP SHARED(nOccBatchDimI,MyMolecule,CoI,iB,nOccBatchDimImax,offset)
          do I=1,nOccBatchDimI
             CoI(:,I) = Mymolecule%Co%elm2(:,offset+I+(iB-1)*nOccBatchDimImax) 
          enddo
          !$OMP END PARALLEL DO

          startjB = 1
          IF(PermutationalSymmetryIJ) startjB = iB
          BatchOccJ: do jB = startjB,nOccbatchesJ
             nOccBatchDimJ = nOccBatchDimJmax
             IF(MOD(nOcc,nOccBatchDimJ).NE.0.AND.jB.EQ.nOccBatchesJ)THEN
                !the remainder
                nOccBatchDimJ = MOD(nOcc,nOccBatchDimJmax)
             ENDIF
             IF(JobsCompletedOrig(iB,jB))CYCLE BatchOccJ
             nJobs = nJobs + 1 
#ifdef VAR_MPI
             ! MPI: Only do this Job if this is a task for this particular rank
             if(MOD(nJobs,numnodes) .NE. mynum) cycle
#endif
             !construct CoJ(nb,nOccBatchDimJ)
             call mem_alloc(CoJ,nb,nOccBatchDimJ)       
             !TODO: OMP Workshare/Loop
             !$OMP PARALLEL DO DEFAULT(none) PRIVATE(J) &
             !$OMP SHARED(nOccBatchDimJ,MyMolecule,CoJ,jB,nOccBatchDimJmax,offset)
             do J=1,nOccBatchDimJ
                CoJ(:,J) = Mymolecule%Co%elm2(:,offset+J+(jB-1)*nOccBatchDimJmax) 
             enddo
             !$OMP END PARALLEL DO
             IF(.NOT.FullRHS)THEN
                call mem_alloc(VOVO,nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
             ENDIF
             BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
                IF(DECinfo%useIchor)THEN
                   dimGamma = AOGammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
                   GammaStart = AOGammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
                   GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
                   AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
                   AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
                ELSE
                   dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
                   GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
                   GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
                ENDIF
                IF(nbatchesAlpha.GT.1)THEN
                   call mem_alloc(VGVO,nvirt*dimGamma,nvirt*nOccBatchDimJ)
                ENDIF
                BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
                   IF(DECinfo%useIchor)THEN
                      dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
                      AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
                      AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
                      AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
                      AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
                   ELSE
                      dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
                      AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
                      AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch
                   ENDIF


                   CALL LS_GETTIM(CPU4,WALL4)
                   IF(DECinfo%useIchor)THEN
                      call mem_alloc(tmp1,dimAlpha*dimGamma,nb*nb)                    
                      !(dimAlpha,dimGamma,nb,nb)
                      call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,dimAlpha,dimGamma,nb,nb,&
                           & tmp1,INTSPEC,FULLRHS,AOAlphaStart,AOAlphaEnd,AOGammaStart,AOGammaEnd,&
                           & 1,nAObatches,1,nAObatches,MoTrans,dimAlpha,dimGamma,nb,nb,NoSymmetry,DECinfo%IntegralThreshold)
                      CALL LS_GETTIM(CPU3,WALL3)
                      CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
                      WALL_AOINT = WALL_AOINT + (WALL3-WALL4)
                   ELSE
                      call mem_alloc(tmp1,dimAlpha*dimGamma,nb*nb)
                      IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%masterGabRHS
                      IF(doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p
                      call II_GET_DECPACKED4CENTER_J_ERI2(DECinfo%output,DECinfo%output, &
                           & mylsitem%setting,tmp1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
                           & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),dimAlpha,&
                           & dimGamma,nbasis,nbasis,FullRHS,INTSPEC,DECinfo%IntegralThreshold)
                      CALL LS_GETTIM(CPU3,WALL3)
                      CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
                      WALL_AOINT = WALL_AOINT + (WALL3-WALL4)
                   ENDIF

                   !tmp2(dimAlpha,dimGamma,nb,nOccBatchDimJ) = tmp1(dimAlpha,dimGamma,nb,nb)*CoJ(nb,nOccBatchDimJ)
                   call mem_alloc(tmp2,dimAlpha*dimGamma,nb*nOccBatchDimJ)
                   M = dimAlpha*dimGamma*nb     !rows of Output Matrix
                   N = nOccBatchDimJ            !columns of Output Matrix
                   K = nb                       !summation dimension
                   call dgemm('N','N',M,N,K,1.0E0_realk,tmp1,M,CoJ,K,0.0E0_realk,tmp2,M)
                   call mem_dealloc(tmp1)

                   !reorder: tmp3(nb,nOccBatchDimJ,dimAlpha,dimGamma) = tmp2(dimAlpha,dimGamma,nb,nOccBatchDimJ)
                   call mem_alloc(tmp3,nb*nOccBatchDimJ,dimAlpha*dimGamma)
                   M = dimAlpha*dimGamma   !row of Input Matrix
                   N = nb*nOccBatchDimJ    !columns of Input Matrix
                   call mat_transpose(M,N,1.0E0_realk,tmp2,0.0E0_realk,tmp3)
                   call mem_dealloc(tmp2)

                   !tmp4(nvirt,nOccBatchDimJ,dimAlpha,dimGamma) = Cv(nb,nvirt)*tmp3(nb,nOccBatchDimJ,dimAlpha,dimGamma)
                   call mem_alloc(tmp4,nvirt*nOccBatchDimJ,dimAlpha*dimGamma)
                   M = nvirt                            !rows of Output Matrix
                   N = nOccBatchDimJ*dimAlpha*dimGamma  !columns of Output Matrix
                   K = nb                               !summation dimension
                   call dgemm('T','N',M,N,K,1.0E0_realk,MyMolecule%Cv%elm2,K,tmp3,K,0.0E0_realk,tmp4,M)
                   call mem_dealloc(tmp3)

                   !reorder: tmp5(dimAlpha,dimGamma,nvirt,nOccBatchDimJ) <= tmp4(nvirt,nOccBatchDimJ,dimAlpha,dimGamma)
                   call mem_alloc(tmp5,dimAlpha*dimGamma,nvirt*nOccBatchDimJ)
                   M = nvirt*nOccBatchDimJ   !row of Input Matrix
                   N = dimAlpha*dimGamma     !columns of Input Matrix
                   call mat_transpose(M,N,1.0E0_realk,tmp4,0.0E0_realk,tmp5)
                   call mem_dealloc(tmp4)

                   call mem_alloc(CvA,dimAlpha,nvirt)
                   do B=1,nvirt
                      CvA(1:dimAlpha,B) = Mymolecule%Cv%elm2(AlphaStart:AlphaEnd,B)
                   enddo
                   !VGVO(nvirt,dimGamma,nvirt,nOccBatchDimJ) = CvA(dimAlpha,nvirt)*tmp5(dimAlpha,dimGamma,nvirt,nOccBatchDimJ)
                   IF(nbatchesAlpha.EQ.1)THEN
                      call mem_alloc(VGVO,nvirt*dimGamma,nvirt*nOccBatchDimJ)
                   ENDIF
                   M = nvirt                         !rows of Output Matrix
                   N = dimGamma*nvirt*nOccBatchDimJ  !columns of Output Matrix
                   K = dimAlpha                      !summation dimension
                   IF(alphaB.EQ.1)THEN
                      call dgemm('T','N',M,N,K,1.0E0_realk,CvA,K,tmp5,K,0.0E0_realk,VGVO,M) 
                   ELSE
                      call dgemm('T','N',M,N,K,1.0E0_realk,CvA,K,tmp5,K,1.0E0_realk,VGVO,M) !ADD
                   ENDIF
                   call mem_dealloc(tmp5)
                   call mem_dealloc(CvA)

                   CALL LS_GETTIM(CPU4,WALL4)
                   CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
                   WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
                ENDDO BatchAlpha
                CALL LS_GETTIM(CPU3,WALL3)

                !reorder: tmp7(nvirt,nOccBatchDimJ,nvirt,dimGamma) <= VGVO(nvirt,dimGamma,nvirt,nOccBatchDimJ)
                call mem_alloc(tmp7,nvirt*nOccBatchDimJ,nvirt*dimGamma)
                M = nvirt*dimGamma          !row of Input Matrix
                N = nvirt*nOccBatchDimJ     !columns of Input Matrix
                call mat_transpose(M,N,1.0E0_realk,VGVO,0.0E0_realk,tmp7)
                call mem_dealloc(VGVO)

                call mem_alloc(CoIG,dimGamma,nOccBatchDimI)
                !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) SHARED(nOccBatchDimI,CoI,CoIG,dimGamma,GammaStart,GammaEnd)              
                do I=1,nOccBatchDimI
                   CoIG(1:dimGamma,I) = CoI(GammaStart:GammaEnd,I)
                enddo
                !$OMP END PARALLEL DO
                !VOVO(nvirt,nOccBatchDimJ,nvirt,nOccBatchDimI) = tmp7(nvirt,nOccBatchDimJ,nvirt,dimGamma)*CoIG(dimGamma,nOccBatchDimI)
                IF(FullRHS)THEN
                   call mem_alloc(VOVO,nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
                ENDIF
                M = nvirt*nOccBatchDimJ*nvirt     !rows of Output Matrix
                N = nOccBatchDimI                 !columns of Output Matrix
                K = dimGamma                      !summation dimension
                IF(gammaB.EQ.1)THEN
                   call dgemm('N','N',M,N,K,1.0E0_realk,tmp7,M,CoIG,K,0.0E0_realk,VOVO,M) 
                ELSE
                   call dgemm('N','N',M,N,K,1.0E0_realk,tmp7,M,CoIG,K,1.0E0_realk,VOVO,M) !ADD TO VOVO
                ENDIF
                call mem_dealloc(tmp7)     

                call mem_dealloc(CoIG)

                CALL LS_GETTIM(CPU4,WALL4)
                CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
                WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
             ENDDO BatchGamma
             CALL LS_GETTIM(CPU3,WALL3)
             call mem_dealloc(CoJ)
             call mem_alloc(Amat,nvirt,nvirt)
             call mem_alloc(Bmat,nvirt,nvirt)
             IF(PermutationalSymmetryIJ)THEN
                WRITE(DECinfo%output,*)'Occ Contribution iB=',iB,', jB=',jB
                IF(iB.NE.jB)THEN 
                   tmp_mp2_energy = 0.0E0_realk
                   do I=1,nOccBatchDimI
                      do J=1,nOccBatchDimJ
                         epsIJ = EpsOcc(I+(iB-1)*nOccBatchDimImax) + EpsOcc(J+(jB-1)*nOccBatchDimJmax)
                         CALL CalcAmat2(nOccBatchDimJ,nOccBatchDimI,nvirt,VOVO,Amat,J,I)
                         call CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
                         tmp_mp2_energy2 = 0.0E0_realk
                         call MP2_EnergyContribution(nvirt,Amat,Bmat,tmp_mp2_energy2)
                         WRITE(DECinfo%output,*)'E1(',iB,',',jB,')=',tmp_mp2_energy2
                         tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
                      enddo
                   enddo
                   !all these contributions appear twice due to permutational symmetry ( I <-> J )
                   tmp_mp2_energy = 2.0E0_realk*tmp_mp2_energy
                   WRITE(DECinfo%output,*)'PermutationalSymmetryIJ E(Triangular)',tmp_mp2_energy
                ELSE !iB = jB same block 
                   tmp_mp2_energy = 0.0E0_realk
                   do I=1,nOccBatchDimI
                      do J=I+1,nOccBatchDimJ
                         epsIJ = EpsOcc(I+(iB-1)*nOccBatchDimImax) + EpsOcc(J+(jB-1)*nOccBatchDimJmax)
                         CALL CalcAmat2(nOccBatchDimJ,nOccBatchDimI,nvirt,VOVO,Amat,J,I)
                         call CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
                         tmp_mp2_energy2 = 0.0E0_realk
                         call MP2_EnergyContribution(nvirt,Amat,Bmat,tmp_mp2_energy2)
                         WRITE(DECinfo%output,*)'E2(',iB,',',jB,')=',tmp_mp2_energy2
                         tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
                      enddo
                   enddo
                   !all these contributions appear twice due to permutational symmetry ( I <-> J )
                   tmp_mp2_energy = 2.0E0_realk*tmp_mp2_energy
                   WRITE(DECinfo%output,*)'PermutationalSymmetryIJ E(Triangular)',tmp_mp2_energy
                   !all these contributions only appear once since I=J in diagonal (iB,iB) block
                   do I=1,nOccBatchDimI
                      epsIJ = 2.0E0_realk*EpsOcc(I+(iB-1)*nOccBatchDimImax)
                      CALL CalcAmat2(nOccBatchDimI,nOccBatchDimI,nvirt,VOVO,Amat,I,I)
                      call CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
                      tmp_mp2_energy2 = 0.0E0_realk
                      call MP2_EnergyContribution(nvirt,Amat,Bmat,tmp_mp2_energy2)
                      WRITE(DECinfo%output,*)'E3(',iB,',',jB,')=',tmp_mp2_energy2
                      tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
                   enddo
                   WRITE(DECinfo%output,*)'PermutationalSymmetryIJ E(Triangular+diagonal)',tmp_mp2_energy
                ENDIF
                WRITE(DECinfo%output,*)'canon MP2 energy contribution =',tmp_mp2_energy
             ELSE
                tmp_mp2_energy = 0.0E0_realk
                do I=1,nOccBatchDimI
                   do J=1,nOccBatchDimJ
                      epsIJ = EpsOcc(I+(iB-1)*nOccBatchDimImax) + EpsOcc(J+(jB-1)*nOccBatchDimJmax)
                      CALL CalcAmat2(nOccBatchDimJ,nOccBatchDimI,nvirt,VOVO,Amat,J,I)
                      call CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
                      tmp_mp2_energy2 = 0.0E0_realk
                      call MP2_EnergyContribution(nvirt,Amat,Bmat,tmp_mp2_energy2)
!                      WRITE(DECinfo%output,*)'E4(I=',I+(iB-1)*nOccBatchDimImax,',J=',J+(jB-1)*nOccBatchDimJmax,')=',tmp_mp2_energy2
                      tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
                   enddo
                enddo
!                WRITE(DECinfo%output,*)'E4(iB=',iB,',jB=',jB,')=',tmp_mp2_energy
!                WRITE(DECinfo%output,*)'canon MP2 energy contribution =',tmp_mp2_energy,'Mynum=',mynum
             ENDIF

             call mem_dealloc(Amat)
             call mem_dealloc(Bmat)
             call mem_dealloc(VOVO)
             CALL LS_GETTIM(CPU4,WALL4)
             CPU_ECONT = CPU_ECONT + (CPU4-CPU3)
             WALL_ECONT = WALL_ECONT + (WALL4-WALL3)           
#ifdef VAR_MPI
             IF(master)THEN
                mp2_energy = mp2_energy + tmp_mp2_energy
                JobsCompleted(iB,jB) = .TRUE.       
                IF(PermutationalSymmetryIJ) JobsCompleted(jB,iB) = .TRUE.       
                !test to see if job info have been recieved from slaves
                MessageRecieved = .TRUE.
                MessageRecievedW = .TRUE.
                CALL LS_GETTIM(CPU1,WALL1)
                DO WHILE(MessageRecievedW)
                   call lsmpi_iprobe(comm,MessageRecieved,mpistatus)
                   MessageRecievedW = MessageRecieved
                   IF(MessageRecievedW)THEN
                      !get the sender ID
                      senderID = mpistatus(MPI_SOURCE)
                      nMPI = 2
                      call lsmpi_recv(JobInfo1,nMPI,senderID,TAG1,comm)
                      call lsmpi_recv(tmp_mp2_energy,senderID,TAG2,comm)
                      
!                      WRITE(DECinfo%output,*)'canon MP2 energy contribution =',tmp_mp2_energy,'From rank=',senderID
                      mp2_energy = mp2_energy + tmp_mp2_energy
                      JobsCompleted(JobInfo1(1),JobInfo1(2)) = .TRUE.
                      IF(PermutationalSymmetryIJ) JobsCompleted(JobInfo1(2),JobInfo1(1)) = .TRUE.
                   ENDIF
                ENDDO
                CALL LS_GETTIM(CPU2,WALL2)
                CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
                WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
             ELSE
                masterrank = 0
                !send info to master
                IF(JobInfo1Free)THEN
                   !wait for master to recieve the first mp2_energy
                   CALL LS_GETTIM(CPU1,WALL1)
                   call lsmpi_wait(request1)
                   call lsmpi_wait(request2)
                   CALL LS_GETTIM(CPU2,WALL2)
                   CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
                   WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
                ENDIF
                nMPI = 2
                JobInfo1(1) = iB
                JobInfo1(2) = jB
                call lsmpi_Isend(JobInfo1,nMPI,masterrank,TAG1,comm,request1)
                mp2_energy = tmp_mp2_energy
                call lsmpi_Isend(mp2_energy,masterrank,TAG2,comm,request2)
                JobInfo1Free = .TRUE.
             ENDIF
#else
             mp2_energy = mp2_energy + tmp_mp2_energy
             JobsCompleted(iB,jB) = .TRUE.
             IF(PermutationalSymmetryIJ) JobsCompleted(jB,iB) = .TRUE.
#endif
             IF(master)THEN
                !Restart File 
                restart_lun = -1  !initialization
                call lsopen(restart_lun,'FULLMP2.restart','UNKNOWN','UNFORMATTED')
                rewind restart_lun
                write(restart_lun) nOccbatchesI,nOccbatchesJ
                write(restart_lun) JobsCompleted
                write(restart_lun) mp2_energy
                call lsclose(restart_lun,'KEEP')
             ENDIF
          enddo BatchOccJ !Batched Occupied J
          call mem_dealloc(CoI)       
       enddo BatchOccI !Batched Occupied I

#ifdef VAR_MPI
       !Wait for all slaves to be finished
       IF(master)THEN
          NotAllMessagesRecieved = COUNT(JobsCompleted).NE.nOccbatchesI*nOccbatchesJ
          DO WHILE(NotAllMessagesRecieved)
             call lsmpi_probe(comm,mpistatus)
             senderID = mpistatus(MPI_SOURCE)
             nMPI = 2
             CALL LS_GETTIM(CPU1,WALL1)
             call lsmpi_recv(JobInfo1,nMPI,senderID,TAG1,comm)
             CALL LS_GETTIM(CPU2,WALL2)
             CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
             WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
             call lsmpi_recv(tmp_mp2_energy,senderID,TAG2,comm)
             CALL LS_GETTIM(CPU1,WALL1)
             CPU_MPICOMM = CPU_MPICOMM + (CPU1-CPU2)
             WALL_MPICOMM = WALL_MPICOMM + (WALL1-WALL2)
             mp2_energy = mp2_energy + tmp_mp2_energy
             JobsCompleted(JobInfo1(1),JobInfo1(2)) = .TRUE.             
             IF(PermutationalSymmetryIJ) JobsCompleted(JobInfo1(2),JobInfo1(1)) = .TRUE.
             NotAllMessagesRecieved = COUNT(JobsCompleted).NE.nOccbatchesI*nOccbatchesJ
          ENDDO
          !Restart File 
          restart_lun = -1  !initialization
          call lsopen(restart_lun,'FULLMP2.restart','UNKNOWN','UNFORMATTED')
          rewind restart_lun
          write(restart_lun) nOccbatchesI,nOccbatchesJ
          write(restart_lun) JobsCompleted
          write(restart_lun) mp2_energy
          call lsclose(restart_lun,'KEEP')
       ENDIF
#endif

       IF(DECinfo%useIchor)THEN
          call FREE_SCREEN_ICHORERI()
          call mem_dealloc(AOGammabatchinfo)
          call mem_dealloc(AOAlphabatchinfo)
       ELSE
          nullify(mylsitem%setting%LST_GAB_LHS)
          nullify(mylsitem%setting%LST_GAB_RHS)
          call free_decscreen(DECSCREEN)
          ! Free gamma batch stuff
          call mem_dealloc(orb2batchGamma)
          call mem_dealloc(batchdimGamma)
          call mem_dealloc(batchsizeGamma)
          call mem_dealloc(batchindexGamma)
          orb2batchGamma => null()
          batchdimGamma => null()
          batchsizeGamma => null()
          batchindexGamma => null()
          do idx=1,nbatchesGamma
             call mem_dealloc(batch2orbGamma(idx)%orbindex)
             batch2orbGamma(idx)%orbindex => null()
          end do

          call mem_dealloc(batch2orbGamma)
          batch2orbGamma => null()

          ! Free alpha batch stuff
          call mem_dealloc(orb2batchAlpha)
          call mem_dealloc(batchdimAlpha)
          call mem_dealloc(batchsizeAlpha)
          call mem_dealloc(batchindexAlpha)
          orb2batchAlpha => null()
          batchdimAlpha => null()
          batchsizeAlpha => null()
          batchindexAlpha => null()
          do idx=1,nbatchesAlpha
             call mem_dealloc(batch2orbAlpha(idx)%orbindex)
             batch2orbAlpha(idx)%orbindex => null()
          end do
          call mem_dealloc(batch2orbAlpha)
          batch2orbAlpha => null()
       ENDIF
       call mem_dealloc(EpsOcc)
       call mem_dealloc(EpsVirt)
    ENDIF
    call mem_dealloc(JobsCompleted)
    call mem_dealloc(JobsCompletedOrig)
    
    IF(MASTER)THEN
       write(lupri,*)  ''
       write(lupri,*)  'MP2 CORRELATION ENERGY = ', mp2_energy
       write(*,'(1X,a,f20.10)') 'MP2 CORRELATION ENERGY = ', mp2_energy
       write(lupri,*)  ''
    ENDIF
    CALL LSTIMER('FULL CANONICAL MP2 ',TS,TE,DECINFO%OUTPUT,FORCEPRINT)
    CALL ls_TIMTXT('>>>  WALL Time used in MP2 AO integral evaluation ',WALL_AOINT,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MP2 AO integral evaluation  ',CPU_AOINT,lupri)
    CALL ls_TIMTXT('>>>  WALL Time used in MP2 AO to MO transformation',WALL_AOTOMO,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MP2 AO to MO transformation ',CPU_AOTOMO,lupri)
    CALL ls_TIMTXT('>>>  WALL Time used in MP2 Energy evaluation      ',WALL_ECONT,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MP2 Energy evaluation       ',CPU_ECONT,lupri)
#ifdef VAR_MPI
    write(lupri,*)'Overall Time spent in MPI Communication and MPI Wait for rank=',infpar%mynum
    CALL ls_TIMTXT('>>>  WALL Time used MPI Communication inc. some Wait',WALL_MPICOMM,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used MPI Communication inc. some Wait',CPU_MPICOMM,lupri)
    CALL ls_TIMTXT('>>>  WALL Time used in MPI Wait',WALL_MPIWAIT,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MPI Wait',CPU_MPIWAIT,lupri)
#endif    
    write(lupri,*) ' '
  end subroutine full_canonical_mp2

  subroutine get_optimal_batch_sizes_for_canonical_mp2(MinAObatch,nbasis,nocc,nvirt,&
       & MaxAllowedDimAlpha,MaxAllowedDimGamma,nOccBatchDimI,nOccBatchDimJ,numnodes,lupri)
    implicit none
    integer,intent(in) :: MinAObatch,nbasis,nocc,nvirt,lupri,numnodes
    integer,intent(inout) :: MaxAllowedDimAlpha,MaxAllowedDimGamma
    integer,intent(inout) :: nOccBatchDimI,nOccBatchDimJ
    !local variables
    real(realk) :: MemoryAvailable,GB,AG,maxsize
    integer :: iB,jB,iiB,jjB,nb,dimGamma,dimAlpha,AB,iGB,nTasks
    real(realk) :: nOccTMP,nbasisTMP
    logical :: BOTH
    ! Memory currently available
    ! **************************
    call get_currently_available_memory(MemoryAvailable)
    ! Note: We multiply by 85 % to be on the safe side!
    MemoryAvailable = 0.85E0_realk*MemoryAvailable
    GB = 1.000E-9_realk 
    nb = nbasis
    nOccTMP = nocc
    nbasisTMP = nbasis
    !test if full is possible
    MaxAllowedDimAlpha = -1
    MaxAllowedDimGamma = -1
    if(DECinfo%manual_occbatchsizes)then
       nOccBatchDimI = DECinfo%batchOccI
       nOccBatchDimJ = DECinfo%batchOccJ
       WRITE(DECinfo%output,*)'CANONMP2MEM: Occupied batches chosen in input:',nOccBatchDimI,nOccBatchDimJ
       IF(DECinfo%manual_batchsizes)THEN
          MaxAllowedDimAlpha = DECinfo%ccsdAbatch
          MaxAllowedDimGamma = DECinfo%ccsdGbatch
          WRITE(DECinfo%output,*)'CANONMP2MEM: AO batches chosen in input:',MaxAllowedDimAlpha,MaxAllowedDimGamma
       ELSE
          iB = nOccBatchDimI
          jB = nOccBatchDimJ
          BatchAlpha3: do AB = 1,nbasis
             dimAlpha = MAX(MinAObatch,CEILING(nbasisTMP/AB)) !(nbasis/1,nbasis/2,..)
             BatchGamma3: do iGB = 1,nbasis
                dimGamma = MAX(MinAObatch,CEILING(nbasisTMP/iGB)) !(nbasis/1,nbasis/2,..)
                call memestimateCANONMP2(iB,jB,dimAlpha,dimGamma,nb,nvirt,maxsize)
                MaxAllowedDimAlpha = dimAlpha
                MaxAllowedDimGamma = dimGamma
                IF(maxsize.LT.MemoryAvailable)THEN
                   WRITE(DECinfo%output,*)'CANONMP2MEM: mem estimate means a batching of Occupied Index I'
                   WRITE(DECinfo%output,*)'Estimated Memory usage is ',maxsize,' GB'
                   exit BatchAlpha3
                ENDIF
                IF(dimGamma.EQ.MinAObatch)exit BatchGamma3
             enddo BatchGamma3
             IF(dimAlpha.EQ.MinAObatch)exit BatchAlpha3
          enddo BatchAlpha3
          WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
       ENDIF
    else
       WRITE(DECinfo%output,*)'Test Full N**4 possible'
       call memestimateCANONMP2(nOcc,nOcc,nb,nb,nb,nvirt,maxsize)
       IF(maxsize.LT.MemoryAvailable)THEN
          MaxAllowedDimAlpha = nb
          MaxAllowedDimGamma = nb
          nOccBatchDimI = nocc
          nOccBatchDimJ = nocc
          WRITE(DECinfo%output,*)'CANONMP2MEM: Final Occupied batches chosen:',nOccBatchDimI,nOccBatchDimJ
          WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
          WRITE(DECinfo%output,*)'CANONMP2MEM: mem estimate means full N**4 dim are used'
          WRITE(DECinfo%output,*)'Estimated Memory is ',maxsize,' GB'
       ELSE
          WRITE(DECinfo%output,*)'Full scheme would require',maxsize,' GB'
          IF(numnodes.GT.1)THEN
             WRITE(DECinfo%output,*)'Test scheme reducing the size of the nOccBatchDimI since numnodes=',numnodes
             !reduce the size of the nOccBatchDimI
             iB = CEILING(nOccTMP/numnodes)
             !assume minimum AO batches dimGamma = MinAObatch, dimAlpha = MinAObatch
             call memestimateCANONMP2(iB,nOcc,MinAObatch,MinAObatch,nb,nvirt,maxsize)
             IF(maxsize.LT.0.5E0_realk*MemoryAvailable)THEN
                nOccBatchDimI = iB
                nOccBatchDimJ = nocc
                WRITE(DECinfo%output,*)'CANONMP2MEM: Final Occupied batches chosen:',nOccBatchDimI,nOccBatchDimJ
                IF(DECinfo%manual_batchsizes)THEN
                   MaxAllowedDimAlpha = DECinfo%ccsdAbatch
                   MaxAllowedDimGamma = DECinfo%ccsdGbatch
                ELSE
                   BatchAlpha: do AB = 1,nbasis
                      dimAlpha = MAX(MinAObatch,CEILING(nbasisTMP/AB)) !(nbasis/1,nbasis/2,..)
                      BatchGamma: do iGB = 1,nbasis
                         dimGamma = MAX(MinAObatch,CEILING(nbasisTMP/iGB)) !(nbasis/1,nbasis/2,..)
                         call memestimateCANONMP2(iB,nOcc,dimAlpha,dimGamma,nb,nvirt,maxsize)
                         MaxAllowedDimAlpha = dimAlpha
                         MaxAllowedDimGamma = dimGamma
                         IF(maxsize.LT.MemoryAvailable)THEN
                            WRITE(DECinfo%output,*)'CANONMP2MEM: mem estimate means a batching of Occupied Index I'
                            WRITE(DECinfo%output,*)'Estimated Memory usage is ',maxsize,' GB'
                            exit BatchAlpha
                         ENDIF
                         IF(dimGamma.EQ.MinAObatch)exit BatchGamma
                      enddo BatchGamma
                      IF(dimAlpha.EQ.MinAObatch)exit BatchAlpha
                   enddo BatchAlpha
                   WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
                ENDIF
                BOTH = .FALSE.
             ELSE
                BOTH = .TRUE.
             ENDIF
          ELSE
             BOTH = .TRUE.
          ENDIF
          IF(BOTH)THEN
             WRITE(DECinfo%output,*)'Test scheme reducing both the size of the nOccBatchDimI and J'
             !reduce the size of the nOccBatchDimI and nOccBatchDimJ
             OccLoopI: do iiB = 1,nOcc
                iB = nOcc/iiB !(nOcc/2,nOcc/3)
                jB = iB
                nTasks = iiB*iiB 
                IF(nTasks.LT.numnodes)THEN
                   print*,'nTasks=',ntasks,'.LT.numnodes=',numnodes,' CYCLE '
                   CYCLE OccLoopI
                ENDIF
                WRITE(DECinfo%output,*)'Test Reduction in OccI,OccJ : nOccBatchDimI',IB,'nOccBatchDimJ',JB
                call memestimateCANONMP2(iB,jB,MinAObatch,MinAObatch,nb,nvirt,maxsize)
                WRITE(DECinfo%output,*)'mem for MinAObatch,MinAObatch:',maxsize,'GB'
                IF(maxsize.LT.0.5E0_realk*MemoryAvailable)THEN
                   nOccBatchDimI = iB
                   nOccBatchDimJ = jB
                   WRITE(DECinfo%output,*)'CANONMP2MEM: Final Occupied batches chosen:',nOccBatchDimI,nOccBatchDimJ
                   IF(DECinfo%manual_batchsizes)THEN
                      MaxAllowedDimAlpha = DECinfo%ccsdAbatch
                      MaxAllowedDimGamma = DECinfo%ccsdGbatch
                   ELSE
                      BatchAlpha2: do AB = 1,nbasis
                         dimAlpha = MAX(MinAObatch,CEILING(nbasisTMP/AB)) !(nbasis/1,nbasis/2,..)
                         BatchGamma2: do iGB = 1,nbasis
                            dimGamma = MAX(MinAObatch,CEILING(nbasisTMP/iGB)) !(nbasis/1,nbasis/2,..)
                            WRITE(DECinfo%output,*)'Obtain Alpha Gamma',dimAlpha,dimGamma
                            call memestimateCANONMP2(iB,jB,dimAlpha,dimGamma,nb,nvirt,maxsize)
                            MaxAllowedDimAlpha = dimAlpha
                            MaxAllowedDimGamma = dimGamma
                            IF(maxsize.LT.MemoryAvailable)THEN
                               WRITE(DECinfo%output,*)'CANONMP2MEM: mem estimate means a batching of Occupied Index I and J'
                               WRITE(DECinfo%output,*)'Estimated Memory usage is ',maxsize,' GB'
                               exit BatchAlpha2
                            ENDIF
                            IF(dimGamma.EQ.MinAObatch)exit BatchGamma2
                         enddo BatchGamma2
                         IF(dimAlpha.EQ.MinAObatch)exit BatchAlpha2
                      enddo BatchAlpha2
                      WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
                   ENDIF
                   EXIT OccLoopI
                ELSE
                   IF(iiB.EQ.nOcc)THEN
                      nOccBatchDimI = 1
                      nOccBatchDimJ = 1
                      MaxAllowedDimAlpha = MinAObatch
                      MaxAllowedDimGamma = MinAObatch
                      WRITE(DECinfo%output,*)'CANONMP2MEM: Final Occupied batches chosen:',nOccBatchDimI,nOccBatchDimJ
                      WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
                   ENDIF
                ENDIF
             enddo OccLoopI
          ENDIF
       ENDIF
    ENDIF
    IF(MaxAllowedDimAlpha.LT.MinAObatch)THEN
       call lsquit('Not enough memory in MP2',-1)
    ENDIF
    IF(MaxAllowedDimGamma.LT.MinAObatch)THEN
       call lsquit('Not enough memory in MP2',-1)
    ENDIF

  end subroutine get_optimal_batch_sizes_for_canonical_mp2

  subroutine memestimateCANONMP2(iB,jB,dimAlpha,dimGamma,nbasis,nvirt,maxsize)
    implicit none
    integer,intent(in) :: iB,jB,dimAlpha,dimGamma,nvirt,nbasis
    real(realk),intent(inout) :: maxsize
    !
    real(realk) :: step1,step2,step3,step4,step5,step6,step7,step8,step9
    real(realk) :: step10,step11,step12,GB,nb
    nb=nbasis
    GB = 1.000E-9_realk 
    !construct Co(nb,nOccBatchDimI)
    step1 = nbasis*iB
    !CoJ(nb,nOccBatchDimJ)       
    step2 = nbasis*iB + nbasis*jB

    IF(dimGamma.EQ.nbasis.AND.dimAlpha.EQ.nbasis)THEN
       step3 = nbasis*iB + nbasis*jB !noallocation of VOVO
    ELSE
       !VOVO(nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
       step3 = nbasis*iB + nbasis*jB + nvirt*nvirt*iB*jB
    ENDIF
    IF(dimAlpha.NE.nbasis)THEN
       !VGVO(nvirt*dimGamma,nvirt*nOccBatchDimJ)
       step4 = step3 + jB*nvirt*nvirt*dimGamma
    ELSE
       step4 = step3 !noallocation of VGVO
    ENDIF
    !tmp1(dimAlpha*dimGamma,nb*nb)
    step5 = step4 + dimAlpha*dimGamma*nb*nb
    !tmp2(dimAlpha,dimGamma,nb,nOccBatchDimJ)
    step6 = step4 + dimAlpha*dimGamma*nb*nb + dimAlpha*dimGamma*nb*jB
    !call mem_alloc(tmp3,nb*nOccBatchDimJ,dimAlpha*dimGamma)
    step7 = step4 + dimAlpha*dimGamma*nb*jB + nb*jB*dimAlpha*dimGamma
    !call mem_alloc(tmp4,nvirt*nOccBatchDimJ,dimAlpha*dimGamma)
    step8 = step4 + nb*jB*dimAlpha*dimGamma + nvirt*jB*dimAlpha*dimGamma
    !call mem_alloc(tmp5,dimAlpha*dimGamma,nvirt*nOccBatchDimJ)
    step9 = step4 + nvirt*jB*dimAlpha*dimGamma + dimAlpha*dimGamma*nvirt*jB
    IF(dimAlpha.EQ.nbasis)THEN
       !call mem_alloc(CvA,dimAlpha,nvirt)+VGVO(nvirt*dimGamma,nvirt*nOccBatchDimJ)
       step10 = step4 + dimAlpha*dimGamma*nvirt*jB + dimAlpha*nvirt + jB*nvirt*nvirt*dimGamma
    ELSE
       !call mem_alloc(CvA,dimAlpha,nvirt)
       step10 = step4 + dimAlpha*dimGamma*nvirt*jB + dimAlpha*nvirt
    ENDIF
    !call mem_alloc(tmp7,nvirt*nOccBatchDimJ,nvirt*dimGamma)
    step11 = step3 + nvirt*jB*nvirt*dimGamma 
    IF(dimGamma.EQ.nbasis.AND.dimAlpha.NE.nbasis)THEN
       !call mem_alloc(CoIG,dimGamma,nOccBatchDimI)
       step12 = step3 + nvirt*jB*nvirt*dimGamma + dimGamma*iB
    ELSE
       !alloc VOVO(nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI) + CoIG
       step12 = step3 + nvirt*jB*nvirt*dimGamma + dimGamma*iB + nvirt*nvirt*iB*jB
    ENDIF
    maxsize = MAX(step1,step2,step3,step4,step5,step6,step7,step8,step9,step10,step11,step12)*realk*GB
    !  WRITE(DECinfo%output,*)'MemRequired Step  1:',Step1*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step  2:',Step2*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step  3:',Step3*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step  4:',Step4*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step  5:',Step5*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step  6:',Step6*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step  7:',Step7*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step  8:',Step8*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step  9:',Step9*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step 10:',Step10*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step 11:',Step11*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'MemRequired Step 12:',Step12*realk*GB,' GB'
    !  WRITE(DECinfo%output,*)'DECinfo%memory',DECinfo%memory
  end subroutine memestimateCANONMP2

  subroutine CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
    implicit none
    integer,intent(in) :: nvirt
    real(realk),intent(in) :: EpsIJ,Amat(nvirt,nvirt),EpsVirt(nvirt)
    real(realk),intent(inout) :: Bmat(nvirt,nvirt)
    !
    integer :: A,B
    !permutational symmetry ?
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
    !$OMP PRIVATE(B,A) SHARED(nvirt,Amat,Bmat,epsIJ,EpsVirt)
    do B=1,nvirt        
       do A=1,nvirt
          Bmat(A,B) =(2.0E0_realk*Amat(A,B)-Amat(B,A))/(epsIJ-EpsVirt(A)-EpsVirt(B))
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine CalcBmat

  !FIXME USE A DOT - is that faster?
  subroutine MP2_EnergyContribution(nvirt,Amat,Bmat,rimp2_energy)
    implicit none
    integer,intent(in) :: nvirt
    real(realk),intent(in) :: Amat(nvirt*nvirt),Bmat(nvirt*nvirt)
    real(realk),intent(inout) :: rimp2_energy
    !
    integer :: A
    real(realk) :: TMP

    TMP = 0.0E0_realk
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP PRIVATE(A) REDUCTION(+:TMP) SHARED(nvirt,Amat,Bmat)
    DO A=1,nvirt*nvirt
       TMP = TMP + Amat(A)*Bmat(A)      
    ENDDO
    !$OMP END PARALLEL DO
    rimp2_energy = rimp2_energy + TMP
  end subroutine MP2_EnergyContribution

  subroutine CalcAmat2(nOccBatchDimI,nOccBatchDimJ,nvirt,tmp7,Amat,I,J)
    implicit none
    integer,intent(in) :: nOccBatchDimI,nOccBatchDimJ,nvirt,I,J
    real(realk),intent(in) :: tmp7(nvirt,nOccBatchDimI,nvirt,nOccBatchDimJ)
    real(realk),intent(inout) :: Amat(nvirt,nvirt)
    !local variables
    integer :: A,B
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
    !$OMP PRIVATE(B,A) SHARED(nvirt,tmp7,Amat,I,J)
    do B=1,nvirt
       do A=1,nvirt
          Amat(A,B) = tmp7(A,I,B,J) 
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine CalcAmat2

  !> \brief Full canonical MP2 calculation, not particularly efficient, mainly to be used for
  !> testing.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine full_canonical_mp2_correlation_energy(MyMolecule,mylsitem,Ecorr)

    implicit none
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS Dalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Canonical MP2 correlation energy
    real(realk),intent(inout) :: Ecorr
    real(realk),pointer :: Cocc(:,:), Cvirt(:,:)
    type(array4) :: g
    integer :: nbasis,i,j,a,b,ncore,offset,nocc,nvirt
    real(realk) :: eps
    real(realk), pointer :: ppfock(:,:)

    ! Sanity check
    if(.not. DECinfo%use_canonical) then
       call lsquit('full_canonical_mp2_correlation_energy requires canonical orbitals! &
            & Insert .CANONICAL keyword OR insert .PRINTFRAGS keyword to run test calculation,&
            & where the individual fragment energies are calculated',-1)
    end if

    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(full_canonical_mp2_correlation_energy): does not work&
       & with PDM distributed molecule structure",-1)
    endif

    ! Initialize stuff
    ! ****************

    if(DECinfo%frozencore) then
       ! Frozen core: Only valence orbitals
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%nocc
    end if

    nvirt = MyMolecule%nvirt
    ncore = MyMolecule%ncore
    nbasis=MyMolecule%nbasis
    call mem_alloc(ppfock,nocc,nocc)
    if(DECinfo%frozencore) then
       ! Only copy valence orbitals into array2 structure
       call mem_alloc(Cocc,nbasis,nocc)
       do i=1,nocc
          Cocc(:,i) = MyMolecule%Co%elm2(:,i+Ncore)
       end do

       ! Fock valence
       do j=1,nocc
          do i=1,nocc
             ppfock(i,j) = MyMolecule%oofock%elm2(i+Ncore,j+Ncore)
          end do
       end do
       offset = ncore
    else
       ! No frozen core, simply copy elements for all occupied orbitals
       call mem_alloc(Cocc,nbasis,nocc)
       Cocc=MyMolecule%Co%elm2
       ppfock = MyMolecule%oofock%elm2
       offset=0
    end if
    call mem_alloc(Cvirt,nbasis,nvirt)
    Cvirt = MyMolecule%Cv%elm2

    ! Get (AI|BJ) integrals stored in the order (A,I,B,J)
    ! ***************************************************
    call get_VOVO_integrals(mylsitem,nbasis,nocc,nvirt,Cvirt,Cocc,g)
    call mem_dealloc(Cocc)
    call mem_dealloc(Cvirt)

    ! Calculate canonical MP2 energy
    ! ******************************
    Ecorr = 0.0_realk
    do J=1,nocc
       do B=1,nvirt
          do I=1,nocc
             do A=1,nvirt
                ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
                eps = MyMolecule%oofock%elm2(I+offset,I+offset) + MyMolecule%oofock%elm2(J+offset,J+offset) &
                     & - MyMolecule%vvfock%elm2(A,A) - MyMolecule%vvfock%elm2(B,B)

                ! Ecorr = sum_{IJAB} (AI|BJ) * [ 2*(AI|BJ) - (BI|AJ) ] / [eps(I)+eps(J)-eps(A)-eps(B)]
                Ecorr = Ecorr + g%val(A,I,B,J)*(2E0_realk*g%val(A,I,B,J)-g%val(B,I,A,J))/eps
             enddo
          enddo
       enddo
    enddo

    call mem_dealloc(ppfock)
    call array4_free(g)

  end subroutine Full_canonical_mp2_correlation_energy

end module fullmp2

#ifdef VAR_MPI
subroutine full_canonical_mp2_slave
  use fullmp2,only: full_canonical_mp2,full_canonical_mpmp2
  use infpar_module !infpar
  use lsmpi_type,only:ls_mpiInitBuffer,ls_mpiFinalizeBuffer,&
       & LSMPIBROADCAST,MPI_COMM_LSDALTON 
  use lsmpi_op,only: mpicopy_lsitem
  use precision
  use typedeftype,only:lsitem
  use lsparameters
  use decmpi_module, only: mpi_bcast_fullmolecule
  use DALTONINFO, only: ls_free
!  use typedef
  use dec_typedef_module!,only DECinfo
  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
!  use dec_fragment_utils
  use full_molecule, only:fullmolecule, molecule_finalize
  implicit none
  !> Full molecule info
  type(fullmolecule) :: MyMolecule
  !> Lsitem structure
  type(lsitem) :: mylsitem
  !> Canonical MP2 correlation energy
  real(realk) :: mp2_energy    
  logical :: Success
  
  ! Init MPI buffer
  ! ***************
  ! Main master:  Prepare for writing to buffer
  ! Local master: Receive buffer
  call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,MPI_COMM_LSDALTON)
  ! Integral lsitem
  ! ---------------
  call mpicopy_lsitem(MyLsitem,MPI_COMM_LSDALTON)
  
  call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,MPI_COMM_LSDALTON)
  ! Full molecule bcasting
  ! **********************
  call mpi_bcast_fullmolecule(MyMolecule)
  
  ! Finalize MPI buffer
  ! *******************
  ! Main master:  Send stuff to local masters and deallocate temp. buffers
  ! Local master: Deallocate buffer etc.
  IF(DECinfo%MPMP2)THEN
     call full_canonical_mpmp2(MyMolecule,MyLsitem,mp2_energy)
  ELSE
     call full_canonical_mp2(MyMolecule,MyLsitem,mp2_energy)
  ENDIF
  call ls_free(MyLsitem)
  call molecule_finalize(MyMolecule,.false.)
  
end subroutine full_canonical_mp2_slave
#endif

