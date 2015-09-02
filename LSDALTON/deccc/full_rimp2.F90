!> @file
!> Full calculation
!> This file is mainly a playground for new developments, not intended to be included in a release.

module fullrimp2

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
  use BUILDAOBATCH

  !  DEC DEPENDENCIES (within deccc directory)   
  !  *****************************************
  use dec_fragment_utils
  use array4_simple_operations
  use array3_simple_operations
  use array2_simple_operations
  use ri_util_module
  !  use orbital_operations
  use full_molecule

  public :: full_canonical_rimp2

  private

contains
  !> \brief Memory check for full_canonical_rimp2 subroutine
  !> \author Thomas Kjaergaard
  !> \date October 2014
  subroutine full_canonical_rimp2_memory_check(nbasis,nocc,nvirt,nAux,&
       & numnodes,MinAtomicnAux,MemoryReduced)
    implicit none
    integer,intent(in) :: nbasis,nocc,nvirt,nAux,numnodes,MinAtomicnAux
    real(realk) :: MemRequired,GB,MemStep1,MemStep2
    logical :: MemoryReduced
    GB = 1.0E9_realk

    ! Check that arrays fit in memory (hard-coded)
    MemStep1=2*real(MinAtomicnAux*nbasis*nbasis)+real(nAux*nocc*nvirt)/numnodes
    MemStep2=2*real(nAux*nocc*nvirt)/numnodes + real(nAux*nAux)
    MemRequired = MAX(MemStep1,MemStep2)
!    MemRequired = 2*real(nAux*nocc*nvirt/numnodes) + real(nAux*nAux)
    MemRequired = MemRequired*realk/GB

    if(MemRequired > 0.80E0_realk*DECinfo%memory) then
       print *, 'RIMP2 Have 2 memory intensive steps:'
       print *, 'Step1:'
       print *, '2 arrays of size MinAtomicnAux*nbasis*nbasis=',&
            & MinAtomicnAux*nbasis*nbasis
       print *, '1 array of size (nAux/numnodes)*nocc*nvirt  =',&
            & (nAux/numnodes)*nocc*nvirt
       print *, 'Resulting in =',MemStep1/GB,' GB'
       print *, 'Step1:'
       print *, '2 arrays of size (nAux/numnodes)*nocc*nvirt  =',&
            & (nAux/numnodes)*nocc*nvirt
       print *, '1 array of size nAux*nAux  =',(nAux/numnodes)*nocc*nvirt
       print *, 'Resulting in =',MemStep2/GB,' GB'
       print *, 'With'
       print *, 'MinAtomicnAux=',MinAtomicnAux
       print *, 'nbasis       =',nbasis
       print *, 'nAux         =',nAux
       print *, 'RIMP2 Mem required (GB) = ', MemRequired
       print *, 'We restrict ourselves to 80 % of Max memory (GB)'
       print *, 'RIMP2 Max memory (GB)   = ', DECinfo%memory
       call lsquit('full_canonical_rimp2: Memory exceeded! Try to increase memory &
            & using the .MaxMemory (in GB) keyword in the *DEC section.',-1)
    end if

    IF(MemoryReduced)THEN
       WRITE(DECinfo%output,*)'MemoryReduced RIMP2 scheme chosen as default'
    ELSE
       !determine if MemoryReduced scheme must be used
       MemStep1=2*real(MinAtomicnAux*nbasis*nbasis)
       MemStep1=MemStep1+real(nAux*nocc*nvirt)/numnodes
       MemStep2=4*real(nAux*nocc*nvirt)/numnodes + real(nAux*nAux)
       MemRequired = MAX(MemStep1,MemStep2)
       MemRequired = MemRequired*realk/GB       
       if(MemRequired .GE. 0.80E0_realk*DECinfo%memory) then 
          MemoryReduced = .TRUE.
          WRITE(DECinfo%output,*)'MemoryReduced RIMP2 scheme chosen'
          WRITE(DECinfo%output,*)'MemRequired Step 1:',MemStep1*realk/GB,' GB'
          WRITE(DECinfo%output,*)'MemRequired Step 2:',MemStep2*realk/GB,' GB'
          WRITE(DECinfo%output,*)'DECinfo%memory',DECinfo%memory
       endif
    ENDIF

  end subroutine full_canonical_rimp2_memory_check

  !> \brief Calculate canonical RIMP2 energy for full molecular system
  !> \author Thomas Kjaergaard
  !> \date October 2014
  subroutine full_canonical_rimp2(MyMolecule,MyLsitem,rimp2_energy)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Canonical MP2 correlation energy
    real(realk),intent(inout) :: rimp2_energy    

    integer :: nbasis,nocc,nvirt,naux,noccfull,mynum,numnodes
    logical :: master,wakeslaves
    integer,pointer :: IPVT(:)
    real(realk), pointer   :: work1(:)
    real(realk)            :: RCOND,dummy(2)
    integer(kind=long) :: maxsize
    real(realk) :: tcpuTOT,twallTOT,tcpu_start,twall_start, tcpu_end,twall_end
    real(realk),pointer :: AlphaCD(:,:,:),AlphaCD5(:,:,:),AlphaCD6(:,:,:)
    real(realk),pointer :: Calpha(:),Calpha2(:,:,:),Calpha3(:,:,:)
    real(realk),pointer :: AlphaBeta(:,:),AlphaBeta_minus_sqrt(:,:),TMPAlphaBeta_minus_sqrt(:,:)
    real(realk),pointer :: EpsOcc(:),EpsVirt(:),ABdecomp(:,:)
    logical :: ABdecompCreate
    integer(kind=ls_mpik)  :: COUNT,TAG,IERR,request,Receiver,sender,COUNT2,comm,TAG1,TAG2
    integer :: J,CurrentWait(2),nAwaitDealloc,iAwaitDealloc,I,NBA,OriginalRanknauxMPI
    integer :: myOriginalRank,node,natoms,A,lupri,MinAtomicnAux,restart_lun
    integer :: noccJstart,offset
    logical :: useAlphaCD5,useAlphaCD6,MessageRecieved,RoundRobin,RoundRobin5,RoundRobin6
    logical :: RoundRobin2,RoundRobin3,useCalpha2,useCalpha3,FORCEPRINT
    integer(kind=ls_mpik)  :: request1,request2,request5,request6,request7,request8
    integer,pointer :: nbasisauxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
    integer(KIND=long) :: MaxMemAllocated,MemAllocated
    real(realk) :: TS,TE,TS2,TE2,TS3,TE3,CPU1,CPU2,WALL1,WALL2,CPU_MPICOMM,WALL_MPICOMM
    real(realk) :: CPU_MPIWAIT,WALL_MPIWAIT,MemInGBCollected,epsIJ
    logical :: MemoryReduced,AlphaCDAlloced,AlphaCD_Deallocate
    logical(kind=ls_mpik) :: TransferCompleted
    logical :: NotMatSet,file_exists,PerformCalc
    real(realk),pointer :: Amat(:,:),Bmat(:,:)
    character :: intspec(5)
    PerformCalc = .TRUE.
    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(full_canonical_rimp2): does not work with distributed&
       & molecular structure",-1)
    endif

    MemoryReduced = MyLsitem%setting%scheme%ForceRIMP2memReduced
    FORCEPRINT = .TRUE.
    CALL LSTIMER('START ',TS,TE,DECINFO%OUTPUT,FORCEPRINT)
#ifdef VAR_TIME
    FORCEPRINT = .TRUE.
#else
    FORCEPRINT = .FALSE.
#endif    
    CPU_MPICOMM = 0.0E0_realk
    WALL_MPICOMM = 0.0E0_realk
    CPU_MPIWAIT = 0.0E0_realk
    WALL_MPIWAIT = 0.0E0_realk
    !use Memory leak tool
    MaxMemAllocated = 0
    MemAllocated = 0
    call set_LeakTool_memvar(MaxMemAllocated,MemAllocated)
    TAG = 1411
    TAG1 = 1412
    TAG2 = 1413
    !sanity check
    if(.NOT.DECinfo%use_canonical) then
       call lsquit('Error: full_canonical_rimp2 require canonical Orbitals',-1)
    endif
    ! Init stuff
    ! **********
    nbasis = MyMolecule%nbasis
    nvirt  = MyMolecule%nvirt
    naux   = MyMolecule%nauxbasis
    !MyMolecule%Co is allocated (nbasis,MyMolecule%nocc)
    !with MyMolecule%nocc = Core + Valence 
    !In case of Frozen core we only need Valence and will access
    !Co(nbasis,ncore+nval) = MyMolecule%Co(nbasis,MyMolecule%ncore+1:MyMolecule%nocc)
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
    IF(master.AND.DECinfo%DECrestart)THEN
       !CHECK IF THERE ARE ENERGY CONTRIBUTIONS AVAILABLE
       INQUIRE(FILE='FULLRIMP2.restart',EXIST=file_exists)
       IF(file_exists)THEN
          IF(master)THEN
             WRITE(DECinfo%output,*)'Restart of Full molecular RIMP2 calculation:'
          ENDIF
          restart_lun = -1  !initialization
          call lsopen(restart_lun,'FULLRIMP2.restart','OLD','FORMATTED')
          rewind restart_lun
          read(restart_lun,'(I9)') noccJstart
          IF(noccJstart.EQ.nocc)THEN
             IF(master)WRITE(DECinfo%output,*)'All energies is on file'
             PerformCalc = .FALSE.
             noccJstart = nocc+1
             read(restart_lun,'(F28.16)') rimp2_energy
          ELSEIF(noccJstart.GT.nocc.OR.noccJstart.LT.1)THEN
             IF(master)THEN
                WRITE(DECinfo%output,*)'RIMP2 restart error first integer is wrong. Read:',noccJstart
             ENDIF
             call lsquit('RIMP2 restart error first integer is wrong')             
          ELSE
             noccJstart = noccJstart + 1
             read(restart_lun,'(F28.16)') rimp2_energy
          ENDIF
          call lsclose(restart_lun,'KEEP')
       ELSE
          noccJstart=1
          rimp2_energy = 0.0E0_realk
       ENDIF
    ELSE
       noccJstart=1
       rimp2_energy = 0.0E0_realk
#ifdef VAR_MPI
       IF(DECinfo%DECrestart)THEN
          call ls_mpibcast(noccJstart,infpar%master,infpar%lg_comm)
       ENDIF
#endif
    ENDIF

    IF(PerformCalc)THEN
       ! Memory check!
       ! ********************
       CALL LSTIMER('START ',TS2,TE2,LUPRI,FORCEPRINT)
       !    call getMaxAtomicnAux(Mylsitem%SETTING%MOLECULE(1)%p,MaxAtomicnAux,nAtoms)
       call determine_maxBatchOrbitalsize(DECinfo%output,&
            & Mylsitem%setting,MinAtomicnAux,'D')
       call full_canonical_rimp2_memory_check(nbasis,nocc,nvirt,nAux,&
            & numnodes,MinAtomicnAux,MemoryReduced)
       call Test_if_64bit_integer_required(nAux,nAux)
       CALL LSTIMER('RIMP2: MemCheck ',TS2,TE2,LUPRI,FORCEPRINT)
#ifdef VAR_MPI 
       ! Master starts up slave
       StartUpSlaves: if(wakeslaves .and. master) then
          ! Wake up slaves to do the job: slaves awoken up with (RIMP2FULL)
          ! and call full_canonical_rimp2_slave which communicate info 
          ! then calls full_canonical_rimp2.
          CALL LS_GETTIM(CPU1,WALL1)
          call ls_mpibcast(RIMP2FULL,infpar%master,comm)
          ! Communicate fragment information to slaves
          call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
          call mpicopy_lsitem(MyLsitem,comm)
          call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
          call mpi_bcast_fullmolecule(MyMolecule)    
          CALL LS_GETTIM(CPU2,WALL2)
          CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
          WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
          call ls_mpibcast(noccJstart,infpar%master,infpar%lg_comm)
       endif StartUpSlaves
#endif
       CALL LSTIMER('RIMP2: WakeSlaves ',TS2,TE2,LUPRI,FORCEPRINT)    
       
       call mem_alloc(ABdecomp,nAux,nAux)
       ABdecompCreate = .TRUE.
       intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
       intspec(2) = 'R' !Regular AO basis function on center 3
       intspec(3) = 'R' !Regular AO basis function on center 4
       intspec(4) = 'C' !Coulomb Operator
       intspec(5) = 'C' !Coulomb Operator
       call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
            & FORCEPRINT,wakeslaves,MyMolecule%Cv%elm2,nvirt,&
            & MyMolecule%Co%elm2(:,offset+1:offset+nocc),nocc,mynum,numnodes,&
            & Calpha,NBA,ABdecomp,ABdecompCreate,intspec,.FALSE.)
       !    PRINT*,'Build_CalphaMO2  nbasis,nAux,nvirt,nocc,NBA',NBA
       !    WRITE(6,*)'Final Calph2(NBA=',NBA,',nvirt=',nvirt,',nocc=',nocc,')'
       !    WRITE(6,*)'Print Subset Final Calph2(NBA=',NBA,',1:4)  MYNUM',MYNUM
       !    call ls_output(Calpha,1,NBA,1,4,NBA,nvirt*nocc,1,6)
       call mem_dealloc(ABdecomp)
       
       call mem_alloc(EpsOcc,nocc)
       call mem_leaktool_alloc(EpsOcc,LT_Eps)
       !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
       !$OMP SHARED(nocc,MyMolecule,EpsOcc,offset)
       do I=1,nocc
          EpsOcc(I) = MyMolecule%oofock%elm2(I+offset,I+offset)
       enddo
       !$OMP END PARALLEL DO
       call mem_alloc(EpsVirt,nvirt)
       call mem_leaktool_alloc(EpsVirt,LT_Eps)
       !$OMP PARALLEL DO DEFAULT(none) PRIVATE(A) &
       !$OMP SHARED(nvirt,MyMolecule,EpsVirt)
       do A=1,nvirt
          EpsVirt(A) = MyMolecule%vvfock%elm2(A,A)
       enddo
       !$OMP END PARALLEL DO
       
       IF(Wakeslaves)THEN
#ifdef VAR_MPI
          NotMatSet = .TRUE.
          call mem_alloc(Amat,nvirt,nvirt)
          call mem_alloc(Bmat,nvirt,nvirt)
          do J=noccJstart,nocc
             do I=1,nocc
                epsIJ = EpsOcc(I) + EpsOcc(J)
                IF(NBA.GT.0)THEN 
                   !A(nvirt,nvirt)
                   CALL CalcAmat(nocc,nvirt,NBA,Calpha,Amat,I,J)
                   CALL CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
                ELSE
                   IF(NotMatSet)THEN
                      NotMatSet = .FALSE.
                      call ls_dzero(Amat,nvirt*nvirt) 
                      call ls_dzero(Bmat,nvirt*nvirt) 
                   ENDIF
                ENDIF
                !The barrier is mostly here to detect the time spent 
                !waiting vs the time spent in communication. 
                CALL LS_GETTIM(CPU1,WALL1)
                call lsmpi_barrier(comm)
                CALL LS_GETTIM(CPU2,WALL2)
                CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
                WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
                !Reduce A,B
                call lsmpi_reduction(Amat,nvirt,nvirt,infpar%master,comm)
                call lsmpi_reduction(Bmat,nvirt,nvirt,infpar%master,comm)
                CALL LS_GETTIM(CPU1,WALL1)
                CPU_MPICOMM = CPU_MPICOMM + (CPU1-CPU2)
                WALL_MPICOMM = WALL_MPICOMM + (WALL1-WALL2)             
                IF(master)THEN
                   call MP2_EnergyContribution(nvirt,Amat,Bmat,rimp2_energy)
                ENDIF
             enddo
             !Write Restart File 
             restart_lun = -1  !initialization
             call lsopen(restart_lun,'FULLRIMP2.restart','UNKNOWN','FORMATTED')
             rewind restart_lun
             write(restart_lun,'(I9)') J
             write(restart_lun,'(F28.16)') rimp2_energy
             call lsclose(restart_lun,'KEEP')
          enddo
          call mem_dealloc(Amat)
          call mem_dealloc(Bmat)
          !       call mem_dealloc(nbasisauxMPI)
          !       call mem_dealloc(startAuxMPI)
          !       call mem_dealloc(nAtomsMPI)
          !       call mem_dealloc(nAuxMPI)
          IF(NBA.GT.0)THEN 
             call mem_dealloc(Calpha)              
          ENDIF
#endif
       ELSE
          !Energy = sum_{AIBJ} (AI|BJ)_N*[ 2(AI|BJ)_N - (BI|AJ)_N ]/(epsI+epsJ-epsA-epsB)
          call RIMP2_CalcOwnEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,&
               & NBA,Calpha,rimp2_energy)
          call mem_dealloc(Calpha)              
       ENDIF
       CALL LSTIMER('RIMP2: EnergyCont ',TS2,TE2,LUPRI,FORCEPRINT)
       call mem_leaktool_dealloc(EpsOcc,LT_Eps)
       call mem_dealloc(EpsOcc)
       call mem_leaktool_dealloc(EpsVirt,LT_Eps)
       call mem_dealloc(EpsVirt)
    ENDIF
    IF(MASTER)THEN
       write(lupri,*)  'RIMP2 CORRELATION ENERGY = ', rimp2_energy
       print*,'RIMP2 CORRELATION ENERGY = ', rimp2_energy
    ENDIF
    write(lupri,*)  'LEAK TOOL STATISTICS IN full_canonical_rimp2'
    call LeakTools_stat_mem(lupri)
    CALL LSTIMER('RIMP2: Finalize ',TS2,TE2,LUPRI,FORCEPRINT)
    FORCEPRINT = .TRUE.
    CALL LSTIMER('FULL RIMP2 ',TS,TE,DECINFO%OUTPUT,FORCEPRINT)
#ifdef VAR_MPI
    write(lupri,*)'Overall Time spent in MPI Communication and MPI Wait for rank=',infpar%mynum
    CALL ls_TIMTXT('>>>  WALL Time used MPI Communication inc. some Wait',WALL_MPICOMM,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used MPI Communication inc. some Wait',CPU_MPICOMM,lupri)
    CALL ls_TIMTXT('>>>  WALL Time used in MPI Wait',WALL_MPIWAIT,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MPI Wait',CPU_MPIWAIT,lupri)
#endif    
    write(lupri,*) ' '
  end subroutine full_canonical_rimp2

  subroutine CalcAmat(nocc,nvirt,NBA,Calpha,Amat,I,J)
    implicit none
    integer,intent(in) :: nocc,nvirt,NBA,I,J
    real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
    real(realk),intent(inout) :: Amat(nvirt,nvirt)
    !
    integer :: A,B,ALPHA
    real(realk) :: TMP
    !permutational symmetry
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
    !$OMP PRIVATE(B,A,ALPHA,TMP) SHARED(nocc,nvirt,NBA,Calpha,Amat,I,J)
    do B=1,nvirt        
       do A=1,nvirt      
          TMP = 0.0E0_realk
          DO ALPHA = 1,NBA
             TMP = TMP + Calpha(ALPHA,A,I)*Calpha(ALPHA,B,J)
          ENDDO
          Amat(A,B) = TMP
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine CalcAmat

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

! Calculate canonical MP2 energy
subroutine RIMP2_CalcOwnEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,NBA,Calpha,rimp2_energy)
  implicit none
  integer,intent(in) :: nocc,nvirt,NBA
  real(realk),intent(in) :: EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(inout) :: rimp2_energy
  !
  integer :: J,B,A,I,ALPHA
  real(realk) :: eps,gmoAIBJ,gmoBIAJ,TMP,epsIJB
  Tmp = 0.0E0_realk
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
  !$OMP PRIVATE(J,B,A,I,eps,gmoAIBJ,gmoBIAJ,epsIJB,ALPHA) REDUCTION(+:TMP) &
  !$OMP SHARED(nocc,nvirt,EpsOcc,EpsVirt,NBA,Calpha)
  do J=1,nocc
     do B=1,nvirt        
        epsIJB = EpsOcc(J) + EpsOcc(J) - EpsVirt(B)
        !==================================================================
        ! I=J AND A=B
        !==================================================================
        ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
        eps = epsIJB - EpsVirt(B)
        gmoAIBJ = 0.0E0_realk
        DO ALPHA = 1,NBA
           gmoAIBJ = gmoAIBJ + Calpha(ALPHA,B,J)*Calpha(ALPHA,B,J)
        ENDDO
        !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
        Tmp = Tmp + gmoAIBJ*gmoAIBJ/eps
        !==================================================================
        ! I=J AND A>B
        !==================================================================
        do A=B+1,nvirt
           ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
           eps = epsIJB - EpsVirt(A)
           gmoAIBJ = 0.0E0_realk
           DO ALPHA = 1,NBA
              gmoAIBJ = gmoAIBJ + Calpha(ALPHA,A,J)*Calpha(ALPHA,B,J)
           ENDDO
           !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
           Tmp = Tmp + 2.0E0_realk*gmoAIBJ*gmoAIBJ/eps
        end do
        do I=J+1,nocc
           epsIJB = EpsOcc(I) + EpsOcc(J) - EpsVirt(B)
           !==================================================================
           ! I>J AND A=B
           !==================================================================
           ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
           eps = epsIJB - EpsVirt(B)
           gmoAIBJ = 0.0E0_realk
           DO ALPHA = 1,NBA
              gmoAIBJ = gmoAIBJ + Calpha(ALPHA,B,I)*Calpha(ALPHA,B,J)
           ENDDO
           !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
           Tmp = Tmp + 2.0E0_realk*gmoAIBJ*gmoAIBJ/eps
           !==================================================================
           ! I>J AND A>B
           !==================================================================
           do A=B+1,nvirt
              ! Difference in orbital energies: eps(I) + eps(J) - eps(A) -eps(B)
              eps = epsIJB - EpsVirt(A)
              gmoAIBJ = 0.0E0_realk
              DO ALPHA = 1,NBA
                 gmoAIBJ = gmoAIBJ + Calpha(ALPHA,A,I)*Calpha(ALPHA,B,J)
              ENDDO
              gmoBIAJ = 0.0E0_realk                   
              DO ALPHA = 1,NBA
                 gmoBIAJ = gmoBIAJ + Calpha(ALPHA,B,I)*Calpha(ALPHA,A,J)
              ENDDO
              !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
!              Tmp = Tmp + 2.0E0_realk*gmoAIBJ*(2E0_realk*gmoAIBJ-gmoBIAJ)/eps
!              Tmp = Tmp + 2.0E0_realk*gmoBIAJ*(2E0_realk*gmoBIAJ-gmoAIBJ)/eps
              Tmp = Tmp + 4.0E0_realk*(gmoAIBJ*gmoAIBJ-gmoAIBJ*gmoBIAJ+gmoBIAJ*gmoBIAJ)/eps
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  rimp2_energy = rimp2_energy + Tmp
end subroutine RIMP2_CalcOwnEnergyContribution

! Calculate canonical MP2 energy
!Energy = sum_{AIBJ} (AI|BJ)_K*[ 2(AI|BJ)_N - (BI|AJ)_N ]/(epsI+epsJ-epsA-epsB)
subroutine RIMP2_CalcEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,NBA,&
     & Calpha,Calpha2,NBA2,rimp2_energy)
  implicit none
  integer,intent(in) :: nocc,nvirt,NBA,NBA2
  real(realk),intent(in) :: EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(in) :: Calpha2(NBA2,nvirt,nocc)
  real(realk),intent(inout) :: rimp2_energy
  !
  integer :: J,B,A,I,ALPHA
  real(realk) :: eps,gmoAIBJ,gmoBIAJ,TMP,epsIJB,gmoAIBJ2,gmoBIAJ2
  Tmp = 0.0E0_realk
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
  !$OMP PRIVATE(J,B,A,I,eps,gmoAIBJ,gmoBIAJ,epsIJB,gmoAIBJ2,ALPHA,gmoBIAJ2) &
  !$OMP REDUCTION(+:TMP) &
  !$OMP SHARED(nocc,nvirt,EpsOcc,EpsVirt,NBA,Calpha,NBA2,Calpha2)
  do J=1,nocc
     do B=1,nvirt
        epsIJB = EpsOcc(J) + EpsOcc(J) - EpsVirt(B)
        !================================================================
        ! I = J, A = B       
        !================================================================
        ! Difference in orbital energies: eps(I)+eps(J)-eps(A)-eps(B)
        eps = epsIJB - EpsVirt(B)
        gmoAIBJ2 = 0.0E0_realk
        DO ALPHA = 1,NBA2
           gmoAIBJ2 = gmoAIBJ2 + Calpha2(ALPHA,B,J)*Calpha2(ALPHA,B,J)
        ENDDO
        gmoAIBJ = 0.0E0_realk
        DO ALPHA = 1,NBA
           gmoAIBJ = gmoAIBJ + Calpha(ALPHA,B,J)*Calpha(ALPHA,B,J)
        ENDDO
        !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
        Tmp = Tmp + gmoAIBJ2*gmoAIBJ/eps
        !================================================================
        ! I = J, A > B       
        !================================================================
        do A=B+1,nvirt
           ! Difference in orbital energies: eps(I)+eps(J)-eps(A)-eps(B)
           eps = epsIJB - EpsVirt(A)
           gmoAIBJ2 = 0.0E0_realk
           DO ALPHA = 1,NBA2
              gmoAIBJ2 = gmoAIBJ2 + Calpha2(ALPHA,A,J)*Calpha2(ALPHA,B,J)
           ENDDO
           gmoAIBJ = 0.0E0_realk
           DO ALPHA = 1,NBA
              gmoAIBJ = gmoAIBJ + Calpha(ALPHA,A,J)*Calpha(ALPHA,B,J)
           ENDDO
           !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
           Tmp = Tmp + 2.0E0_realk*gmoAIBJ2*gmoAIBJ/eps
        end do
        !================================================================
        ! I > J, A = B       
        !================================================================
        do I=J+1,nocc
           epsIJB = EpsOcc(I) + EpsOcc(J) - EpsVirt(B)
           ! Difference in orbital energies: eps(I)+eps(J)-eps(A)-eps(B)
           eps = epsIJB - EpsVirt(B)
           gmoAIBJ2 = 0.0E0_realk
           DO ALPHA = 1,NBA2
              gmoAIBJ2 = gmoAIBJ2 + Calpha2(ALPHA,B,I)*Calpha2(ALPHA,B,J)
           ENDDO
           gmoAIBJ = 0.0E0_realk
           DO ALPHA = 1,NBA
              gmoAIBJ = gmoAIBJ + Calpha(ALPHA,B,I)*Calpha(ALPHA,B,J)
           ENDDO
           !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
           Tmp = Tmp + 2.0E0_realk*gmoAIBJ2*gmoAIBJ/eps
           !================================================================
           ! I > J, A > B       
           !================================================================
           do A=B+1,nvirt
              ! Difference in orbital energies: eps(I)+eps(J)-eps(A)-eps(B)
              eps = epsIJB - EpsVirt(A)
              gmoAIBJ2 = 0.0E0_realk
              gmoBIAJ2 = 0.0E0_realk
              DO ALPHA = 1,NBA2
                 gmoAIBJ2 = gmoAIBJ2 + Calpha2(ALPHA,A,I)*Calpha2(ALPHA,B,J)
                 gmoBIAJ2 = gmoBIAJ2 + Calpha2(ALPHA,B,I)*Calpha2(ALPHA,A,J)
              ENDDO
              gmoAIBJ = 0.0E0_realk
              gmoBIAJ = 0.0E0_realk                   
              DO ALPHA = 1,NBA
                 gmoAIBJ = gmoAIBJ + Calpha(ALPHA,A,I)*Calpha(ALPHA,B,J)
                 gmoBIAJ = gmoBIAJ + Calpha(ALPHA,B,I)*Calpha(ALPHA,A,J)
              ENDDO
              !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
              Tmp = Tmp + 2.0E0_realk*gmoAIBJ2*(2E0_realk*gmoAIBJ-gmoBIAJ)/eps
              Tmp = Tmp + 2.0E0_realk*gmoBIAJ2*(2E0_realk*gmoBIAJ-gmoAIBJ)/eps
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  rimp2_energy = rimp2_energy + Tmp
end subroutine RIMP2_CalcEnergyContribution

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

end module fullrimp2

#ifdef VAR_MPI
subroutine full_canonical_rimp2_slave
  use fullrimp2,only: full_canonical_rimp2
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
!  use dec_typedef_module
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
  real(realk) :: rimp2_energy    
  
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
  call full_canonical_rimp2(MyMolecule,MyLsitem,rimp2_energy)

  call ls_free(MyLsitem)
  call molecule_finalize(MyMolecule,.false.)
  
end subroutine full_canonical_rimp2_slave
#endif

