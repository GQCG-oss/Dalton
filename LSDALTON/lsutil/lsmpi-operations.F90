module lsmpi_op
  use precision
  use lstiming
  use molecule_typetype, only: molecularOrbitalInfo
  use io, only: io_init
  use typedeftype, only: DALTONINPUT,LSITEM,lssetting,reducedScreeningInfo,&
       & lsintscheme, integralconfig
  use typedef, only: integral_set_default_config, typedef_init_setting,&
       & init_reduced_screen_info, typedef_free_setting
  use basis_typetype, only: BASISSETINFO,nBasisBasParam,nullifyMainBasis,&
       & nullifyBasisset
  use basis_type, only: lsmpi_alloc_basissetinfo
  use memory_handling, only: mem_alloc,mem_dealloc, mem_shortintsize,&
       & mem_realsize, mem_intsize, mem_allocated_mem_lstensor
  use LSparameters
  use Matrix_Operations, only: mat_mpicopy, mtype_scalapack, matrix_type
  use matrix_operations_scalapack, only: pdm_matrixsync
  use molecule_typetype, only: MOLECULEINFO
  use LSTENSOR_OPERATIONSMOD, only: lstensor_NULLIFY, lstensor,&
       & slsaotensor_nullify,lsaotensor_nullify,&
       & lsaotensor,slsaotensor
  use LSTENSORmem, only: mem_LSTpointer_alloc, init_lstensorMem, &
       & retrieve_lstMemVal, free_lstensorMem, set_lstmemrealkbufferpointer
  use f12_module, only: GaussianGeminal
  use integraloutput_typetype, only: INTEGRALOUTPUT
  use dft_type,only: mpicopy_DFTparam
#ifdef VAR_MPI
  use lsmpi_type, only: ls_mpibcast, lsmpi_reduction, ls_mpisendrecv, &
       & ls_mpi_buffer, lsmpi_local_reduction, get_rank_for_comm, &
       & get_size_for_comm, ls_mpiInitBuffer, ls_MpiFinalizeBuffer, &
       & lsmpi_print, lsmpi_default_mpi_group, lsmpi_finalize, lsmpi_barrier,&
       & LSMPIREDUCTION, LSMPIREDUCTIONmaster, LSMPIBROADCAST,&
       & ls_mpiinitbufferaddtobuffer,printmpibuffersizes,ls_mpiModbuffersizes
  use infpar_module
  use lsmpi_module
#else
  use lsmpi_type, only: ls_mpi_buffer, get_rank_for_comm, get_size_for_comm, &
       & ls_mpiInitBuffer, ls_MpiFinalizeBuffer,LSMPIBROADCAST,&
       & ls_mpiinitbufferaddtobuffer,printmpibuffersizes,ls_mpiModbuffersizes
#endif
  use screen_mod
  !*****************************************
  !*
  !* OBJECTS CONTAINING INFORMATION ABOUT 
  !* THE (MPI) TASK PARTITONING 
  !*   - used for ls_IntegralInterface.f90
  !*
  !*****************************************
  TYPE LSTASK
  LOGICAL :: full_mol_row
  LOGICAL :: full_mol_col
  INTEGER :: nrow
  INTEGER :: ncol
  INTEGER :: nrow_full
  INTEGER :: ncol_full
  REAL(REALK) :: task_part
  INTEGER,POINTER :: row_atoms(:)
  INTEGER,POINTER :: col_atoms(:)
  TYPE(LSTASK),POINTER :: next
  TYPE(LSTASK),POINTER :: previous
  END TYPE LSTASK
  
  TYPE LSTASK_p
  TYPE(LSTASK),POINTER :: p
  END TYPE LSTASK_p
  
  TYPE LS_TASKS
  INTEGER    :: target_ntasks
  INTEGER    :: ntasks
  INTEGER    :: natoms_row
  INTEGER    :: natoms_col
  INTEGER    :: natoms_row_full
  INTEGER    :: natoms_col_full
  LOGICAL    :: no_row_part
  LOGICAL    :: no_col_part
  TYPE(LSTASK),POINTER :: first
  TYPE(LSTASK),POINTER :: last
  END TYPE LS_TASKS

  TYPE LSMPI_TASK_LIST
  INTEGER :: mynode        ! The index of my MPI node
  INTEGER :: numMPItasks   ! The number of tasks to be calculated by the current node
  LOGICAL :: ownLHS        ! Specifies wether all LHS tasks belong to the current node or not
  LOGICAL :: ownRHS        ! Specifies wether all RHS tasks belong to the current node or not
  INTEGER :: nuniqueLHS    ! The number of unique LHS tasks (should be one for (ab|cd)Dcd-type construction)
  INTEGER :: nuniqueRHS    ! The number of unique RHS tasks (should be equal to #nodes for (ab|cd)Dcd-type construction)
  INTEGER,pointer :: taskPair(:,:) !Give index to the LHS- and RHS-task lists
  INTEGER,pointer :: taskNode(:,:) !Give the node the LHS- and RHS-task belongs to
  TYPE(LSTASK_p),pointer :: LHS(:) !LHS-task list
  TYPE(LSTASK_p),pointer :: RHS(:) !RHS-task list
! REAL(realk),pointer  :: taskWall(:)
! REAL(realk),pointer  :: taskCPU(:)
! REAL(realk),pointer  :: taskEstimate(:)
  END TYPE LSMPI_TASK_LIST

  TYPE LS_TASK_MANAGER
  INTEGER :: target_nlhs
  INTEGER :: target_nrhs
  LOGICAL :: lhs_aux
  LOGICAL :: rhs_aux
  LOGICAL :: sameAOsLHS
  LOGICAL :: sameAOsRHS
  LOGICAL :: sameODs
  INTEGER :: itask
  INTEGER :: ntasks
  TYPE(LS_TASKS) :: lhs
  TYPE(LS_TASKS) :: rhs
  type(lstask),pointer :: lhs_current
  type(lstask),pointer :: rhs_current
  TYPE(molecularOrbitalInfo) :: orbInfo(4)
  END TYPE LS_TASK_MANAGER
contains
!> \brief mpi alloc the daltoninput structure
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param dalton the dalton input structure
SUBROUTINE LSMPI_ALLOC_DALTONINPUT(DALTON)
IMPLICIT NONE
TYPE(DALTONINPUT) :: DALTON
integer :: I
! THE MOLECULE
NULLIFY(DALTON%MOLECULE%ATOM)
call mem_alloc(DALTON%MOLECULE%ATOM,DALTON%MOLECULE%nAtoms)
! THE Secondary MOLECULE  (used in DEC for RI/CABS)
NULLIFY(DALTON%AUXMOLECULE%ATOM)
call mem_alloc(DALTON%AUXMOLECULE%ATOM,DALTON%AUXMOLECULE%nAtoms)
!THE BASISSET
do I=1,nBasisBasParam
   CALL LSMPI_ALLOC_BASISSETINFO(DALTON%BASIS%BINFO(I))
enddo
END SUBROUTINE LSMPI_ALLOC_DALTONINPUT

!> \brief mpi copy the lsitem structure - highly discouraged this should not be necessary
!> \author T. Kjaergaard
!> \date 2010
!> \param 
SUBROUTINE mpicopy_lsitem(ls,comm)
implicit none
type(lsitem),intent(inout) :: ls 
! communicator
integer(kind=ls_mpik),intent(in) :: comm
integer(kind=ls_mpik) :: mynum
#ifdef VAR_MPI
call get_rank_for_comm(comm,mynum)
call mpicopy_daltoninput(ls%input,comm)
call mpicopy_setting(ls%setting,comm,.FALSE.)
IF(mynum.NE.infpar%master)THEN !SLAVE
   ls%lupri = -1
   ls%luerr = -1
   ls%fopen = .false.
ENDIF
#endif
END SUBROUTINE mpicopy_lsitem

!> \brief mpi copy the lsitem structure - highly discouraged this should not be necessary
!> \author T. Kjaergaard
!> \date 2010
!> \param 
SUBROUTINE mpicopy_daltoninput(input,comm)
implicit none
type(daltoninput),intent(inout) :: input
integer(kind=ls_mpik),intent(in) :: comm  ! communicator 
#ifdef VAR_MPI
integer(kind=ls_mpik) :: MASTER
LOGICAL :: SLAVE
integer :: I
integer(kind=ls_mpik) :: mynum,nodtot,ierr
call get_rank_for_comm(comm, mynum)
CALL get_size_for_comm(comm, nodtot)
input%numNodes = nodtot 
input%node = mynum

Master = infpar%master
Slave = mynum.NE.infpar%master
CALL LS_MPI_BUFFER(input%DO_DFT,Master) 
CALL LS_MPI_BUFFER(input%POTNUC,Master) 
CALL LS_MPI_BUFFER(input%nfock,Master) 
CALL LS_MPI_BUFFER(input%numfragments,Master) 


IF(SLAVE)THEN
   NULLIFY(input%MOLECULE)
   NULLIFY(input%AUXMOLECULE)
   NULLIFY(input%BASIS)
   ALLOCATE(input%AuxMolecule)
   ALLOCATE(input%Molecule)
   ALLOCATE(input%Basis)
   call nullifyMainBasis(input%Basis)
   CALL io_init(input%IO)
   call integral_set_default_config(input%dalton)
ENDIF
call mpicopy_integralconfig(input%dalton,Slave,Master)
call mpicopy_molecule(input%Molecule,Slave,Master)
call mpicopy_molecule(input%AuxMolecule,Slave,Master)
call LS_MPI_BUFFER(input%basis%WBASIS,nBasisBasParam,Master)
do I=1,nBasisBasParam
   IF(input%basis%WBASIS(I))THEN
      call mpicopy_basissetinfo(input%basis%BINFO(I),Slave,Master)
   ENDIF
enddo
#endif
END SUBROUTINE mpicopy_daltoninput

!> \brief mpi copy the lssetting structure
!> \author T. Kjaergaard
!> \date 2010
!> \param setting the setting structure to broadcast
SUBROUTINE mpicopy_setting(setting,comm,rankslave)
  implicit none
  type(lssetting),intent(inout) :: setting
  integer(kind=ls_mpik),intent(in) :: comm  ! communicator
  logical,intent(in) :: rankslave
  !
  integer :: lupri
  integer :: I,nAO,ndmat,dim1,dim2,dim3,iAO,J
  logical :: SLAVE,nonemptyMolecule,DsymRHSassociated,DsymLHSassociated
  logical :: LHS_GAB,RHS_GAB

  real(realk) :: ts,te
  integer(kind=ls_mpik) :: mynum,nodtot,ierr,master
#ifdef VAR_MPI
  Master  = infpar%master
#else
  Master  = 0  
#endif
  call get_rank_for_comm(comm,mynum)
  CALL get_size_for_comm(comm, nodtot)
  setting%numNodes = nodtot 
  IF(rankslave)THEN
     setting%node = 1
     mynum = 1
  ELSE
     setting%node = mynum     
  ENDIF
  setting%comm = comm

  dim1 = AORdefault
  dim2 = AODFdefault
  CALL LS_MPI_BUFFER(dim1,Master)        
  CALL LS_MPI_BUFFER(dim2,Master)        
  IF(mynum.NE.master)THEN !SLAVE
     call set_default_AOs(dim1,dim2)
  ENDIF


  !call lstimer('START',ts,te,6)
  IF(mynum.NE.master)THEN !SLAVE
     call typedef_init_setting(SETTING)
  ENDIF
  call mpicopy_scheme(setting%SCHEME,Slave,Master)                           !LSSETTING008

  CALL LS_MPI_BUFFER(setting%numFragments,Master)        !LSSETTING074

  CALL LS_MPI_BUFFER(setting%IntegralTransformGC,Master) !LSSETTING001
  CALL LS_MPI_BUFFER(setting%DO_DFT,Master)              !LSSETTING002
  CALL LS_MPI_BUFFER(setting%EDISP,Master)               !LSSETTING003
  CALL LS_MPI_BUFFER(setting%nAO,Master)                 !LSSETTING004
  nAO = setting%nAO

  SLAVE = SETTING%node.ne.Master

  IF(SLAVE)THEN
     setting%molBuild = .TRUE.                                             !LSSETTING018
     setting%basBuild = .TRUE.                                             !LSSETTING019
  ENDIF

  do I = 1,nAO
     IF (.NOT.SLAVE) THEN
        nonemptyMolecule = associated(setting%Molecule(I)%p)
     ENDIF

     CALL LS_MPI_BUFFER(nonemptyMolecule,Master)

     IF(SLAVE)THEN
        nullify(setting%Molecule(I)%p)
        nullify(setting%Basis(I)%p)
     ENDIF
     IF (nonemptyMolecule) THEN
        IF(SLAVE)THEN
           allocate(setting%Molecule(I)%p)
           allocate(setting%Basis(I)%p)
           call nullifyMainBasis(setting%Basis(I)%p)
        ENDIF
        call mpicopy_molecule(setting%Molecule(I)%p,Slave,Master)              !LSSETTING005
        call LS_MPI_BUFFER(setting%basis(I)%p%WBASIS,nBasisBasParam,Master)
        DO J=1,nBasisBasParam
           IF(setting%basis(I)%p%WBASIS(J))THEN
              call mpicopy_basissetinfo(setting%basis(I)%p%BINFO(J),Slave,Master) 
           ENDIF
        ENDDO
     ELSE
        setting%molBuild(I) = .FALSE.                                            !LSSETTING018
        setting%basBuild(I) = .FALSE.                                            !LSSETTING019
     ENDIF
  enddo
  call LS_MPI_BUFFER(setting%Batchindex,nAO,Master)                           !LSSETTING011
  call LS_MPI_BUFFER(setting%Batchsize,nAO,Master)                            !LSSETTING012
  call LS_MPI_BUFFER(setting%Batchdim,nAO,Master)                             !LSSETTING013
  call LS_MPI_BUFFER(setting%molID,nAO,Master)                                !LSSETTING014

  !if (.NOT.SLAVE) call lstimer('copy1',ts,te,6)
  IF(SLAVE)THEN
     setting%fragBuild = .FALSE.                                              !LSSETTING020
     DO iAO=1,nAO
        setting%Fragment(iAO)%p => setting%Molecule(iAO)%p                     !LSSETTING007
     ENDDO
  ENDIF


  call LS_MPI_BUFFER(setting%sameMOL,nAO,nAO,Master)                         !LSSETTING015
  call LS_MPI_BUFFER(setting%sameBas,nAO,nAO,Master)                         !LSSETTING016
  call LS_MPI_BUFFER(setting%sameFrag,nAO,nAO,Master)                        !LSSETTING017

  IF(SLAVE)THEN
     SETTING%LHSdalloc = .FALSE.
     SETTING%RHSdalloc = .FALSE.
     SETTING%LHSdmatAlloc = .FALSE.
     SETTING%RHSdmatAlloc = .FALSE.
     SETTING%lstensor_attached = .FALSE.                                     !LSSETTING029
  ENDIF

  call LS_MPI_BUFFER(setting%nDmatLHS,Master)                                !LSSETTING027
  call LS_MPI_BUFFER(setting%nDmatRHS,Master)                                !LSSETTING028
  call LS_MPI_BUFFER(setting%LHSdmat,Master)                                 !LSSETTING032
  call LS_MPI_BUFFER(setting%RHSdmat,Master)                                 !LSSETTING033
  
  !if (.NOT.SLAVE) call lstimer('copy2',ts,te,6)
  !IF(SETTING%scheme%memdist.AND.matrix_type.EQ.mtype_scalapack)THEN
  !IF(matrix_type.EQ.mtype_scalapack)THEN
  IF(SETTING%LHSdmat)THEN
     ndmat = SETTING%nDmatLHS
     IF(SLAVE)THEN
        call mem_alloc(setting%DmatLHS,ndmat)
        setting%LHSdmatAlloc = .TRUE.                                         !LSSETTING034
     ENDIF
     do I=1,ndmat
        IF(SLAVE)NULLIFY(setting%DmatLHS(I)%p)
        IF(SLAVE)ALLOCATE(setting%DmatLHS(I)%p)
#if VAR_SCALAPACK
        IF(setting%scheme%memdist.AND.matrix_type.EQ.mtype_scalapack)THEN
           call PDM_MATRIXSYNC(setting%DmatLHS(I)%p)
        ELSE
           call mat_mpicopy(setting%DmatLHS(I)%p,Slave,Master)                   !LSSETTING023
        ENDIF
#else
        call mat_mpicopy(setting%DmatLHS(I)%p,Slave,Master)                   !LSSETTING023
#endif
     enddo
  ENDIF
  IF(SETTING%RHSdmat)THEN
     ndmat = SETTING%nDmatRHS
     IF(SLAVE)THEN
        call mem_alloc(setting%DmatRHS,SETTING%nDmatRHS)
        setting%RHSdmatAlloc = .TRUE.                                         !LSSETTING035
     ENDIF
     do I=1,SETTING%nDmatRHS
        IF(SLAVE)NULLIFY(setting%DmatRHS(I)%p)
        IF(SLAVE)ALLOCATE(setting%DmatRHS(I)%p)
#if VAR_SCALAPACK
        IF(setting%scheme%memdist.AND.matrix_type.EQ.mtype_scalapack)THEN
           call PDM_MATRIXSYNC(setting%DmatRHS(I)%p)
        ELSE
           call mat_mpicopy(setting%DmatRHS(I)%p,Slave,Master)                   !LSSETTING024
        ENDIF
#else
        call mat_mpicopy(setting%DmatRHS(I)%p,Slave,Master)                   !LSSETTING024
#endif
     enddo
  ENDIF
  !ENDIF
  !if (.NOT.SLAVE) call lstimer('copy3',ts,te,6)

  DsymLHSassociated = ASSOCIATED(setting%DsymLHS)
  call LS_MPI_BUFFER(DsymLHSassociated,Master) 
  IF(DsymLHSassociated)THEN
     IF(SLAVE)call mem_alloc(setting%DsymLHS,setting%nDmatLHS)
     call LS_MPI_BUFFER(setting%DsymLHS,setting%nDmatLHS,Master)             !LSSETTING068
  ELSE
     IF(SLAVE)nullify(setting%DsymLHS)
  ENDIF

  DsymRHSassociated = ASSOCIATED(setting%DsymRHS)
  call LS_MPI_BUFFER(DsymRHSassociated,Master) 
  IF(DsymRHSassociated)THEN
     IF(SLAVE)call mem_alloc(setting%DsymRHS,setting%nDmatRHS)
     call LS_MPI_BUFFER(setting%DsymRHS,setting%nDmatRHS,Master)             !LSSETTING069
  ELSE
     IF(SLAVE)nullify(setting%DsymRHS)
  ENDIF

  LHS_GAB = associated(setting%LST_GAB_LHS)
  RHS_GAB = associated(setting%LST_GAB_RHS)
  call LS_MPI_BUFFER(LHS_GAB,Master) 
  call LS_MPI_BUFFER(RHS_GAB,Master) 
  IF(SLAVE)THEN
     NULLIFY(SETTING%LST_GAB_LHS)                                             !LSSETTING056
     NULLIFY(SETTING%LST_GAB_RHS)                                             !LSSETTING057
     IF(LHS_GAB)ALLOCATE(SETTING%LST_GAB_LHS)                                 !LSSETTING056
     IF(RHS_GAB)ALLOCATE(SETTING%LST_GAB_RHS)                                 !LSSETTING056
  ENDIF
  !if (.NOT.SLAVE) call lstimer('copy4',ts,te,6)

  IF(LHS_GAB)call mpicopy_lstensor(setting%LST_GAB_LHS,Slave,Master)        !LSSETTING056
  IF(RHS_GAB)call mpicopy_lstensor(setting%LST_GAB_RHS,Slave,Master)        !LSSETTING057
  !if (.NOT.SLAVE) call lstimer('copy5',ts,te,6)

  call LS_MPI_BUFFER(Setting%iLST_GAB_LHS,Master)                            !LSSETTING060
  call LS_MPI_BUFFER(Setting%iLST_GAB_RHS,Master)                            !LSSETTING061
  IF(LHS_GAB)THEN
     call LS_MPI_BUFFER(Setting%CS_MAXELM_LHS,Master)                        !LSSETTING064
     call LS_MPI_BUFFER(Setting%PS_MAXELM_LHS,Master)                        !LSSETTING066
  ENDIF
  IF(RHS_GAB)THEN
     call LS_MPI_BUFFER(Setting%CS_MAXELM_RHS,Master)                        !LSSETTING065
     call LS_MPI_BUFFER(Setting%PS_MAXELM_RHS,Master)                        !LSSETTING067
  ENDIF

  call LS_MPI_BUFFER(SETTING%LHSdfull,Master)                                !LSSETTING036
  IF(SETTING%LHSdfull)THEN
     ndmat = SETTING%nDmatLHS
     IF(.NOT.SLAVE)THEN
        dim1 = SIZE(SETTING%DfullLHS, 1)  
        dim2 = SIZE(SETTING%DfullLHS, 2)  
        dim3 = SIZE(SETTING%DfullLHS, 3)  
     ELSE
        dim1 = 0
        dim2 = 0
        dim3 = 0
     ENDIF
     call LS_MPI_BUFFER(dim1,Master) 
     call LS_MPI_BUFFER(dim2,Master) 
     call LS_MPI_BUFFER(dim3,Master) 
     IF(SLAVE)THEN
        call mem_alloc(SETTING%DfullLHS,dim1,dim2,dim3)
        SETTING%LHSdalloc = .TRUE.                                           !LSSETTING038
     ELSE
     ENDIF
     call LS_MPI_BUFFER(SETTING%DfullLHS,dim1,dim2,dim3,Master)              !LSSETTING025
  ENDIF

  call LS_MPI_BUFFER(SETTING%RHSdfull,Master)                                !LSSETTING037
  IF(SETTING%RHSdfull)THEN
     ndmat = SETTING%nDmatRHS
     IF(.NOT.SLAVE)THEN
        dim1 = SIZE(SETTING%DfullRHS, 1)  
        dim2 = SIZE(SETTING%DfullRHS, 2)  
        dim3 = SIZE(SETTING%DfullRHS, 3)  
     ENDIF
     call LS_MPI_BUFFER(dim1,Master) 
     call LS_MPI_BUFFER(dim2,Master) 
     call LS_MPI_BUFFER(dim3,Master) 
     IF(SLAVE)THEN
        call mem_alloc(SETTING%DfullRHS,dim1,dim2,dim3)
        SETTING%RHSdalloc = .TRUE.
     ENDIF
     call LS_MPI_BUFFER(SETTING%DfullRHS,dim1,dim2,dim3,Master)
  ENDIF

  !if (.NOT.SLAVE) call lstimer('copy6',ts,te,6)

  !if (.NOT.SLAVE) call lstimer('copy6a',ts,te,6)
  call LS_MPI_BUFFER(SETTING%LHSdmatAOindex1,Master)                         !LSSETTING070
  call LS_MPI_BUFFER(SETTING%LHSdmatAOindex2,Master)                         !LSSETTING071
  call LS_MPI_BUFFER(SETTING%RHSdmatAOindex1,Master)                         !LSSETTING072
  call LS_MPI_BUFFER(SETTING%RHSdmatAOindex2,Master)                         !LSSETTING073
  !if (.NOT.SLAVE) call lstimer('copy6b',ts,te,6)
  !if (.NOT.SLAVE) call lstimer('copy6c',ts,te,6)

  call mpicopy_output(setting%output,Slave,Master)                           !LSSETTING078
  !if (.NOT.SLAVE) call lstimer('copy6d',ts,te,6)

  call mpicopy_GaussianGeminal(setting%GGem,slave,master)                    !LSSETTING079

  !if (.NOT.SLAVE) call lstimer('copy7',ts,te,6)

  call mpicopy_reduced_screen_info(setting%redCS,slave,master)               !LSSETTING080

  !if (.NOT.SLAVE) call lstimer('copy8',ts,te,6)
  CALL LS_MPI_BUFFER(setting%GPUMAXMEM,Master)                               !LSSETTING081

END SUBROUTINE mpicopy_setting

Subroutine TestMPIcopySetting(setting)
implicit none
type(lssetting) :: setting
!
integer(kind=ls_mpik) :: master,comm
integer :: Job
logical :: AddToBuffer2,rankslave
Job=LSMPIBROADCAST 
master=0
comm=1
!This will set the addtobuffer (lsmpiType.F90 module variable)
!and allocate the mpibuffers
call ls_mpiInitBuffer(master,Job,comm) 
!place info of setting into buffers (using setting%node = 0 => SLAVE=FALSE)
rankslave = .FALSE.
call mpicopy_setting(setting,comm,rankslave)
CALL typedef_free_setting(setting)
!This will set the addtobuffer (lsmpiType.F90 module variable) to false
AddToBuffer2 = .FALSE.
call ls_mpiInitBufferAddToBuffer(AddToBuffer2)
!place info of buffers into setting (using setting%node = 1 => SLAVE=TRUE)
rankslave = .TRUE.
call mpicopy_setting(setting,comm,rankslave)
!sets the nDP = iDP etc. in lsmpiType.F90 (otherwise an error in ls_MpiFinalizeBuffer)
call ls_mpiModbuffersizes
setting%node = 0
call ls_MpiFinalizeBuffer(master,Job,comm)
end Subroutine TestMPIcopySetting

SUBROUTINE mpicopy_screen(Slave,Master)
implicit none
logical :: slave
integer(kind=ls_mpik) :: master
!
logical :: assStart

IF(SLAVE)call screen_init()
assStart = associated(SCREENFROMLSSCREEN%start)
call LS_MPI_BUFFER(assStart,Master)
IF(assStart)THEN
   IF(SLAVE)THEN
      allocate(SCREENFROMLSSCREEN%start)
      nullify(SCREENFROMLSSCREEN%start%next)
   ENDIF
   call mpicopy_screenchain(SCREENFROMLSSCREEN,SCREENFROMLSSCREEN%start,slave,Master)
ELSE
   IF(SLAVE)nullify(SCREENFROMLSSCREEN%start)
ENDIF
end SUBROUTINE mpicopy_screen

recursive SUBROUTINE mpicopy_screenchain(screen,screenchain,Slave,Master)
implicit none
type(screenitem) :: screen
type(screenchainitem),pointer :: ScreenChain
logical :: slave
integer(kind=ls_mpik) :: master
!
type(LSTENSOR),pointer :: GAB  
logical :: assStart

call LS_MPI_BUFFER(ScreenChain%filename,80,master)
IF(SLAVE)THEN
   nullify(ScreenChain%LST%p)
   allocate(ScreenChain%LST%p)   
ENDIF
call mpicopy_lstensor(ScreenChain%LST%p,slave,master)
assStart = associated(Screenchain%next)
call LS_MPI_BUFFER(assStart,Master)
IF(assStart)THEN
   IF(SLAVE)THEN
      allocate(ScreenChain%next)
      nullify(ScreenChain%next%next)
   ENDIF
   call mpicopy_screenchain(SCREEN,ScreenChain%next,Slave,Master)
ELSE
   IF(SLAVE)THEN
      nullify(Screenchain%next)
      SCREEN%end => Screenchain
   ENDIF
ENDIF
end SUBROUTINE mpicopy_screenchain

Subroutine TestMPIcopyScreen()
implicit none
integer(kind=ls_mpik) :: master,comm
integer :: Job
logical :: AddToBuffer2,slave
Job=LSMPIBROADCAST 
master=0
comm=1
!This will set the addtobuffer (lsmpiType.F90 module variable)
!and allocate the mpibuffers
call ls_mpiInitBuffer(master,Job,comm) 
!place info of setting into buffers (using setting%node = 0 => SLAVE=FALSE)
slave = .FALSE.
call mpicopy_screen(slave,master)
call screen_free
!This will set the addtobuffer (lsmpiType.F90 module variable) to false
AddToBuffer2 = .FALSE.
call ls_mpiInitBufferAddToBuffer(AddToBuffer2)
!place info of buffers into setting (using setting%node = 1 => SLAVE=TRUE)
slave = .TRUE.
call mpicopy_screen(slave,master)
!sets the nDP = iDP etc. in lsmpiType.F90 (otherwise an error in ls_MpiFinalizeBuffer)
call ls_mpiModbuffersizes
call ls_MpiFinalizeBuffer(master,Job,comm)
end Subroutine TestMPIcopyScreen

#ifdef VAR_MPI
subroutine lsmpi_isend_lstmemrealkbuf(lstmem_index,NodeToRecv,Mynum,comm)
implicit none
integer,intent(in) :: lstmem_index,NodeToRecv
integer(kind=ls_mpik),intent(in) :: Mynum,comm
!
integer(kind=long)    :: nbuffer
real(realk),pointer   :: buffer(:)
integer(kind=ls_mpik) :: COUNT,TAG,IERR,request,RCVNODE
TAG = 155534879
RCVNODE=NodeToRecv
call set_lstmemrealkbufferpointer(lstmem_index,buffer,nbuffer)
COUNT = nbuffer
IF(nbuffer.GT.HUGE(COUNT))call lsquit('64 bit error in lsmpi_isend_lstmemrealkbuf',-1)
call MPI_ISEND(buffer,COUNT,MPI_DOUBLE_PRECISION,RCVNODE,TAG,comm,request,ierr)
IF(IERR.NE.MPI_SUCCESS)THEN
   call lsquit('MPI_ISEND ERROR',-1)
!   print*,'MPI_ERR_COMM',MPI_ERR_COMM
!   print*,'MPI_ERR_COUNT',MPI_ERR_COUNT
!   print*,'MPI_ERR_TYPE',MPI_ERR_TYPE
!   print*,'MPI_ERR_TAG',MPI_ERR_TAG
!   print*,'MPI_ERR_RANK',MPI_ERR_RANK
!   print*,'MPI_ERR_INTERN',MPI_ERR_INTERN
ENDIF
call MPI_Request_free(request,ierr)
!call MPI_WAIT(request,ierr)
end subroutine lsmpi_isend_lstmemrealkbuf

subroutine lsmpi_probe_and_irecv_add_lstmemrealkbuf(lstmem_index,Mynum,comm,&
     & inputMessageRecieved,numnodes)
implicit none
integer,intent(in) :: lstmem_index
integer(kind=ls_mpik),intent(in) :: Mynum,comm,numnodes
logical,intent(inout) :: inputMessageRecieved(numnodes)
!
integer(kind=long)  :: nbuffer
real(realk),pointer :: buffer(:)
integer(kind=ls_mpik) :: COUNT,TAG,IERR,request,COUNT2
real(realk),pointer :: buffer2(:)
integer(kind=ls_mpik) :: lsmpi_status(MPI_STATUS_SIZE),nMess,j,i
logical(kind=ls_mpik) :: MessageRecieved
logical :: ALLOC,MessageRecievedW
TAG = 155534879
call set_lstmemrealkbufferpointer(lstmem_index,buffer,nbuffer)
IF(nbuffer.GT.HUGE(COUNT))call lsquit('64 bit error in lsmpi_isend_lstmemrealkbuf',-1)
COUNT = nbuffer
MessageRecieved = .TRUE.
MessageRecievedW = .TRUE.
ALLOC=.FALSE.
DO WHILE(MessageRecievedW)
   call MPI_IPROBE(MPI_ANY_SOURCE,TAG,comm,MessageRecieved,lsmpi_status,ierr)
   MessageRecievedW = MessageRecieved 
   IF(MessageRecievedW)THEN
      IF(.NOT.ALLOC)THEN
         call mem_alloc(buffer2,COUNT)
         ALLOC=.TRUE.
      ENDIF
      !get the sender ID
      j = lsmpi_status(MPI_SOURCE)
      call MPI_GET_COUNT(lsmpi_status, MPI_DOUBLE_PRECISION, COUNT2,IERR)
      IF(COUNT2.NE.COUNT)CALL LSQUIT('COUNT WRONG IN lsmpi_probe_and_irecv_add_lstmemrealkbuf ',-1)
      inputMessageRecieved(j+1) = .TRUE.
      call MPI_IRECV(buffer2,COUNT,MPI_DOUBLE_PRECISION,J,TAG,MPI_COMM_WORLD,request,ierr)
      call MPI_Request_free(request,ierr) 
!      call MPI_WAIT(request,ierr)
      !add to buffer
      do I=1,COUNT
         buffer(I)=buffer(I)+buffer2(I) 
      enddo     
   ENDIF
ENDDO
IF(ALLOC)call mem_dealloc(buffer2)
end subroutine lsmpi_probe_and_irecv_add_lstmemrealkbuf

subroutine lsmpi_blocking_recv_add_lstmemrealkbuf(lstmem_index,Mynum,comm,&
     & inputMessageRecieved,numnodes)
implicit none
integer,intent(in) :: lstmem_index
integer(kind=ls_mpik),intent(in) :: Mynum,comm,numnodes
logical,intent(in) :: inputMessageRecieved(numnodes)
!
integer(kind=long)             :: nbuffer
real(realk),pointer :: buffer(:)
integer(kind=ls_mpik) :: COUNT,TAG,IERR,request
real(realk),pointer :: buffer2(:)
integer(kind=ls_mpik) :: lsmpi_status(MPI_STATUS_SIZE),nMess,k,i,COUNT2,j
integer :: Nmissing
logical :: ALLOC,ALLMessageRecieved
TAG = 155534879
Nmissing=0
ALLMessageRecieved = .TRUE.
DO I=1,numnodes
   IF(.NOT.inputMessageRecieved(I))THEN
      Nmissing = Nmissing + 1
      ALLMessageRecieved = .FALSE.
   ENDIF
ENDDO

call set_lstmemrealkbufferpointer(lstmem_index,buffer,nbuffer)
IF(nbuffer.GT.HUGE(COUNT))call lsquit('64 bit error in lsmpi_isend_lstmemrealkbuf',-1)
COUNT = nbuffer
ALLOC=.FALSE.
DO WHILE(.NOT.ALLMessageRecieved)
   IF(.NOT.ALLOC)THEN
      call mem_alloc(buffer2,COUNT)
      ALLOC=.TRUE.
   ENDIF
   !blocking Recv   
   call MPI_PROBE(MPI_ANY_SOURCE,TAG,comm,lsmpi_status,ierr)
   !get the sender ID
   j = lsmpi_status(MPI_SOURCE)
   call MPI_GET_COUNT(lsmpi_status, MPI_DOUBLE_PRECISION, COUNT2,IERR)
   IF(COUNT2.NE.COUNT)THEN
      print*,'COUNT',COUNT,'infpar%mynm',infpar%mynum
      print*,'COUNT2',COUNT2,'infpar%mynm',infpar%mynum
      call lsquit('COUNT WRONG lsmpi_blocking_recv_add_lstmemrealkbuf',-1)
   ENDIF
   call MPI_RECV(buffer2,COUNT,MPI_DOUBLE_PRECISION,J,TAG,comm,lsmpi_status,ierr)
   !add to buffer
   do K=1,COUNT
      buffer(K)=buffer(K)+buffer2(K) 
   enddo
   Nmissing = Nmissing - 1   
   IF(Nmissing.EQ.0)ALLMessageRecieved = .TRUE.
ENDDO
IF(ALLOC)call mem_dealloc(buffer2)
end subroutine lsmpi_blocking_recv_add_lstmemrealkbuf
#endif
!!$subroutine lsmpi_blocking_recv_add_lstmemrealkbuf2(lstmem_index,Mynum,comm,&
!!$     & inputMessageRecieved,numnodes)
!!$implicit none
!!$integer,intent(in) :: lstmem_index,Mynum,comm,numnodes
!!$logical,intent(inout) :: inputMessageRecieved(numnodes)
!!$!
!!$integer(kind=long)             :: nbuffer
!!$real(realk),pointer :: buffer(:)
!!$integer :: COUNT,TAG,IERR,request
!!$real(realk),pointer :: buffer2(:)
!!$integer :: lsmpi_status(MPI_STATUS_SIZE),nMess,k,i,COUNT2,j
!!$logical :: MessageRecieved,ALLMessageRecieved,ALLOC
!!$
!!$TAG = 155534879
!!$call set_lstmemrealkbufferpointer(lstmem_index,buffer,nbuffer)
!!$IF(nbuffer.GT.HUGE(COUNT))call lsquit('64 bit error in lsmpi_isend_lstmemrealkbuf',-1)
!!$COUNT = nbuffer
!!$ALLOC=.FALSE.
!!$ALLMessageRecieved = .TRUE.
!!$DO I=1,numnodes
!!$   IF(.NOT.inputMessageRecieved(I))THEN
!!$      ALLMessageRecieved = .FALSE.
!!$      EXIT
!!$   ENDIF
!!$ENDDO
!!$
!!$DO WHILE(.NOT.ALLMessageRecieved)
!!$   call MPI_IPROBE(MPI_ANY_SOURCE,TAG,comm,MessageRecieved,lsmpi_status,ierr)
!!$   IF(MessageRecieved)THEN
!!$      IF(.NOT.ALLOC)THEN
!!$         call mem_alloc(buffer2,COUNT)
!!$         ALLOC=.TRUE.
!!$      ENDIF
!!$      !get the sender ID
!!$      j = lsmpi_status(MPI_SOURCE)
!!$      call MPI_GET_COUNT(lsmpi_status, MPI_DOUBLE_PRECISION, COUNT2,IERR)
!!$      IF(COUNT2.NE.COUNT)CALL LSQUIT('COUNT WRONG IN lsmpi_probe_and_irecv_add_lstmemrealkbuf ',-1)
!!$      inputMessageRecieved(j+1) = .TRUE.
!!$      print*,'BLOCKING mynum',mynum,'RECV BUFFER FROM ',J,' OF SIZE',COUNT
!!$      call MPI_IRECV(buffer2,COUNT,MPI_DOUBLE_PRECISION,J,TAG,MPI_COMM_WORLD,request,ierr)
!!$      print*,'IRECV ERROR:',IERR
!!$      print*,'MPI_SUCCESS',MPI_SUCCESS
!!$      IF(IERR.NE.MPI_SUCCESS)THEN
!!$         print*,'IRECV ERROR:',IERR
!!$         print*,'MPI_SUCCESS',MPI_SUCCESS
!!$         print*,'MPI_ERR_COMM',MPI_ERR_COMM
!!$         print*,'MPI_ERR_COUNT',MPI_ERR_COUNT
!!$         print*,'MPI_ERR_TYPE',MPI_ERR_TYPE
!!$         print*,'MPI_ERR_TAG',MPI_ERR_TAG
!!$         print*,'MPI_ERR_RANK',MPI_ERR_RANK
!!$         print*,'MPI_ERR_INTERN',MPI_ERR_INTERN
!!$      ENDIF
!!$
!!$!      call MPI_Request_free(request,ierr) 
!!$      call MPI_WAIT(request,ierr)
!!$      print*,'IWAIT ERROR:',IERR
!!$      print*,'BLOCKING MYNU',mynum,'nBuffer',COUNT
!!$!      print*,'BLOCKING MYNU',mynum,'RECV BUFFER2',buffer2
!!$      !add to buffer
!!$      do I=1,COUNT
!!$         buffer(I)=buffer(I)+buffer2(I) 
!!$      enddo     
!!$!      print*,'BLOCKING MYNU',mynum,'nBuffer',COUNT
!!$!      print*,'BLOCKING MYNU',mynum,'NEW UPDATED BUFFER',buffer
!!$      ALLMessageRecieved = .TRUE.
!!$      DO I=1,numnodes
!!$         IF(.NOT.inputMessageRecieved(I))THEN
!!$            ALLMessageRecieved = .FALSE.
!!$            EXIT
!!$         ENDIF
!!$      ENDDO      
!!$   ENDIF
!!$ENDDO
!!$IF(ALLOC)call mem_dealloc(buffer2)
!!$end subroutine lsmpi_blocking_recv_add_lstmemrealkbuf2

SUBROUTINE mpicopy_lstensor(tensor,Slave,Master)
implicit none
type(lstensor) :: tensor
logical :: slave
integer(kind=ls_mpik) :: master
!
logical    :: isAssociated
INTEGER    :: I,J,K,L,offset,n1,n2,n3,n4
integer(kind=long) :: nmemsize,AllocInt,AllocRealk,AllocIntS,nsize
integer :: AllocInt4,AllocRealk4,AllocIntS4
real(realk) :: ts,te

IF(SLAVE)call lstensor_NULLIFY(TENSOR)

call LS_MPI_BUFFER(TENSOR%natom,4,Master)
call LS_MPI_BUFFER(TENSOR%nbast,4,Master)
call LS_MPI_BUFFER(TENSOR%nbatches,4,Master)
call LS_MPI_BUFFER(TENSOR%ndim5,Master)
call LS_MPI_BUFFER(TENSOR%nSLSAO,Master)
call LS_MPI_BUFFER(TENSOR%nLSAO,Master)
call LS_MPI_BUFFER(TENSOR%Gradienttensor,Master)
call LS_MPI_BUFFER(TENSOR%pChargetensor,Master)
call LS_MPI_BUFFER(TENSOR%Econtrib,Master)
call LS_MPI_BUFFER(TENSOR%Screentensor,Master)
call LS_MPI_BUFFER(TENSOR%nMBIE,Master)
call LS_MPI_BUFFER(TENSOR%LowerDiagZero,Master)
call LS_MPI_BUFFER(TENSOR%PermuteResultTensor,Master)

call LS_MPI_BUFFER(TENSOR%maxgabelm,Master)
call LS_MPI_BUFFER(TENSOR%maxprimgabelm,Master)
IF(.NOT.SLAVE)THEN
   call retrieve_lstmemval(TENSOR%lstmem_index,AllocInt,AllocRealk,AllocIntS)
ENDIF

#ifdef VAR_INT64
!64 bit integers
call LS_MPI_BUFFER(AllocInt,Master)
call LS_MPI_BUFFER(AllocIntS,Master)
call LS_MPI_BUFFER(AllocRealk,Master)
#else
!32 bit integers
IF(.NOT.SLAVE)THEN
   IF(AllocInt.LT.HUGE(offset))THEN
      AllocInt4 = AllocInt
   ELSE
      call lsquit('AllocInt too big requires 64 bit integer bcast',-1)
   ENDIF
   IF(AllocIntS.LT.HUGE(offset))THEN
      AllocIntS4 = AllocIntS
   ELSE
      call lsquit('AllocIntS too big requires 64 bit integer bcast',-1)
   ENDIF
   IF(AllocRealk.LT.HUGE(offset))THEN
      AllocRealk4 = AllocRealk
   ELSE
      call lsquit('AllocInt too big requires 64 bit integer bcast',-1)
   ENDIF
ENDIF
call LS_MPI_BUFFER(AllocInt4,Master)
call LS_MPI_BUFFER(AllocIntS4,Master)
call LS_MPI_BUFFER(AllocRealk4,Master)
AllocInt = AllocInt4
AllocIntS = AllocIntS4
Allocrealk = AllocRealk4   
#endif
IF(SLAVE)THEN
   call init_lstensorMem(AllocInt,AllocRealk,AllocIntS,TENSOR%lstmem_index)
ENDIF
!call mpicopy_lstensorMem(Master)

isAssociated = ASSOCIATED(TENSOR%maxgab)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      call Mem_alloc(TENSOR%maxgab,TENSOR%nbatches(1),TENSOR%nbatches(2))
      nsize = size(TENSOR%maxgab,KIND=long)*mem_shortintsize
      call mem_allocated_mem_lstensor(nsize)
   ENDIF
   call LS_MPI_BUFFER(TENSOR%maxgab,TENSOR%nbatches(1),TENSOR%nbatches(2),Master)
ELSE
   IF(SLAVE)NULLIFY(TENSOR%maxgab)
ENDIF

isAssociated = ASSOCIATED(TENSOR%maxprimgab)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      call mem_alloc(TENSOR%maxprimgab,TENSOR%nbatches(1),TENSOR%nbatches(2))
      nsize = size(TENSOR%maxprimgab,KIND=long)*mem_shortintsize
      call mem_allocated_mem_lstensor(nsize)
   ENDIF
   call LS_MPI_BUFFER(TENSOR%maxprimgab,TENSOR%nbatches(1),TENSOR%nbatches(2),Master)
ELSE
   IF(SLAVE)NULLIFY(TENSOR%maxprimgab)
ENDIF

isAssociated = ASSOCIATED(TENSOR%MBIE)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      call mem_alloc(TENSOR%MBIE,TENSOR%nMBIE,TENSOR%nbatches(1),TENSOR%nbatches(2))
      nsize = size(TENSOR%MBIE,KIND=long)*mem_realsize
      call mem_allocated_mem_lstensor(nsize)
   ENDIF
   call LS_MPI_BUFFER(TENSOR%MBIE,TENSOR%nMBIE,TENSOR%nbatches(1),&
        & TENSOR%nbatches(2),Master)
ELSE
   IF(SLAVE)NULLIFY(TENSOR%MBIE)
ENDIF

isAssociated = ASSOCIATED(TENSOR%nAOBATCH)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      n1 = 0
      n2 = 0
   ELSE
      n1 = SIZE(TENSOR%nAOBATCH,1)
      n2 = SIZE(TENSOR%nAOBATCH,2)
   ENDIF
   call LS_MPI_BUFFER(n1,Master)
   call LS_MPI_BUFFER(n2,Master)
   IF(SLAVE)THEN
      call mem_alloc(TENSOR%nAOBATCH,n1,n2)
      nsize = size(TENSOR%nAOBATCH,KIND=long)*mem_intsize
      call mem_allocated_mem_lstensor(nsize)
   ENDIF
   call LS_MPI_BUFFER(TENSOR%nAOBATCH,n1,n2,Master)
ELSE
   IF(SLAVE)NULLIFY(TENSOR%nAOBATCH)
ENDIF

isAssociated = ASSOCIATED(TENSOR%INDEX)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(TENSOR%Gradienttensor)THEN
      IF(SLAVE)THEN
         n1 = 0
         n2 = 0
         n3 = 0
         n4 = 0
      ELSE
         n1 = SIZE(TENSOR%INDEX,1)
         n2 = SIZE(TENSOR%INDEX,2)
         n3 = SIZE(TENSOR%INDEX,3)
         n4 = SIZE(TENSOR%INDEX,4)
      ENDIF
      call LS_MPI_BUFFER(n1,Master)
      call LS_MPI_BUFFER(n2,Master)
      call LS_MPI_BUFFER(n3,Master)
      call LS_MPI_BUFFER(n4,Master)
      IF(SLAVE)THEN
         IF(n2.NE.1)call lsquit('error in mpicopy_lstensor A.',-1)
         IF(n3.NE.1)call lsquit('error in mpicopy_lstensor B.',-1)
         call mem_alloc(TENSOR%INDEX,n1,n2,n3,n4)
         nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
         call mem_allocated_mem_lstensor(nsize)
         DO L=1,n1
            DO I=1,n4
               TENSOR%INDEX(I,1,1,L) = I
            ENDDO
         ENDDO
      ENDIF
   ELSE
      IF(SLAVE)THEN
         n1 = 0
         n2 = 0
         n3 = 0
         n4 = 0
      ELSE
         n1 = SIZE(TENSOR%INDEX,1)
         n2 = SIZE(TENSOR%INDEX,2)
         n3 = SIZE(TENSOR%INDEX,3)
         n4 = SIZE(TENSOR%INDEX,4)
      ENDIF
      call LS_MPI_BUFFER(n1,Master)
      call LS_MPI_BUFFER(n2,Master)
      call LS_MPI_BUFFER(n3,Master)
      call LS_MPI_BUFFER(n4,Master)
      IF(SLAVE)THEN
         call mem_alloc(TENSOR%INDEX,n1,n2,n3,n4)
         nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
         call mem_allocated_mem_lstensor(nsize)
      ENDIF
      call LS_MPI_BUFFER(TENSOR%INDEX,n1,n2,n3,n4,Master)
   ENDIF
ELSE
   IF(SLAVE)NULLIFY(TENSOR%INDEX)
ENDIF

isAssociated = ASSOCIATED(TENSOR%LSAO)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      call mem_alloc(TENSOR%LSAO,TENSOR%nLSAO)
   ENDIF
   DO I = 1,TENSOR%nLSAO
      call mpicopy_lsaotensor(TENSOR%LSAO(I),Slave,Master)
   ENDDO
ELSE
   IF(SLAVE)NULLIFY(TENSOR%LSAO)
ENDIF

isAssociated = ASSOCIATED(TENSOR%SLSAO)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      call mem_alloc(TENSOR%SLSAO,TENSOR%nSLSAO)
   ENDIF
   DO I = 1,TENSOR%nSLSAO
      call mpicopy_slsaotensor(TENSOR%SLSAO(I),Slave,Master)
   ENDDO
ELSE
   IF(SLAVE)NULLIFY(TENSOR%SLSAO)
ENDIF

!IF(SLAVE)THEN
!   Call Determine_slstensor_memory(tensor,nmemsize)
!   call add_mem_to_global(nmemsize)
!ENDIF
end SUBROUTINE mpicopy_lstensor

subroutine mpicopy_slsaotensor(LSAO,Slave,Master)
implicit none
type(slsaotensor) :: lsao
logical :: slave
integer(kind=ls_mpik) :: master
!
logical    :: isAssociated
INTEGER    :: I,Ibat,Jbat,dim,ielms,maxgab
INTEGER    :: nAngmomA,nAngmomB,maxgabelm,n1
integer(kind=long) :: nmemsize
integer,pointer     :: elms(:)
real(realk) :: ts,te

IF(SLAVE)call slsaotensor_nullify(lsao)
call LS_MPI_BUFFER(lsao%nelms,Master)
call LS_MPI_BUFFER(lsao%nLocal,2,Master)
call LS_MPI_BUFFER(lsao%ATOM,2,Master)
call LS_MPI_BUFFER(lsao%AOBATCH,2,Master)
call LS_MPI_BUFFER(lsao%maxBat,Master)
isAssociated = ASSOCIATED(lsao%Selms)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      call mem_LSTpointer_alloc(lsao%selms,lsao%nelms)
   ENDIF
! we have copied the memory buffer so this should not be needed
   call LS_MPI_BUFFER(lsao%selms,lsao%nelms,Master)
ELSE
   IF(SLAVE)NULLIFY(lsao%selms)
ENDIF

isAssociated = ASSOCIATED(lsao%nOrb)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      n1 = 0
   ELSE
      n1 = size(lsao%nOrb)
   ENDIF
   call LS_MPI_BUFFER(n1,Master)
   IF(SLAVE)THEN
      call mem_LSTpointer_alloc(lsao%nOrb,n1)
   ENDIF
   call LS_MPI_BUFFER(lsao%nOrb,n1,Master)
ELSE
   IF(SLAVE)NULLIFY(lsao%nOrb)
ENDIF

isAssociated = ASSOCIATED(lsao%StartLocalOrb)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      n1 = 0
   ELSE
      n1 = size(lsao%StartLocalOrb,1)
   ENDIF
   call LS_MPI_BUFFER(n1,Master)
   IF(SLAVE)THEN
      call mem_LSTpointer_alloc(lsao%StartLocalOrb,n1)
   ENDIF
   call LS_MPI_BUFFER(lsao%StartLocalOrb,n1,Master)
ELSE
   IF(SLAVE)NULLIFY(lsao%StartLocalOrb)
ENDIF
end subroutine mpicopy_slsaotensor

subroutine mpicopy_lsaotensor(LSAO,Slave,Master)
implicit none
type(lsaotensor) :: lsao
logical :: slave
integer(kind=ls_mpik) :: master
!
logical    :: isAssociated
INTEGER    :: I,Ibat,Jbat,dim,ielms,maxgab
INTEGER    :: nAngmomA,nAngmomB,maxgabelm,n1
integer(kind=long) :: nmemsize
integer,pointer     :: elms(:)

IF(SLAVE)call lsaotensor_nullify(lsao)
call LS_MPI_BUFFER(lsao%nelms,Master)
call LS_MPI_BUFFER(lsao%nLocal,4,Master)
call LS_MPI_BUFFER(lsao%ATOM,4,Master)
call LS_MPI_BUFFER(lsao%AOBATCH,4,Master)
call LS_MPI_BUFFER(lsao%maxBat,Master)
call LS_MPI_BUFFER(lsao%maxAng,Master)
   
isAssociated = ASSOCIATED(lsao%elms)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      call mem_LSTpointer_alloc(lsao%elms,lsao%nelms)
   ENDIF
   call LS_MPI_BUFFER(lsao%elms,lsao%nelms,Master)
ELSE
   IF(SLAVE)NULLIFY(lsao%elms)
ENDIF

isAssociated = ASSOCIATED(lsao%nOrb)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      n1 = 0
   ELSE
      n1 = size(lsao%nOrb,1)
   ENDIF
   call LS_MPI_BUFFER(n1,Master)
   IF(SLAVE)THEN
      call mem_LSTpointer_alloc(lsao%nOrb,n1)
   ENDIF
   call LS_MPI_BUFFER(lsao%nOrb,n1,Master)
ELSE
   IF(SLAVE)NULLIFY(lsao%nOrb)
ENDIF

isAssociated = ASSOCIATED(lsao%StartLocalOrb)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      n1 = 0
   ELSE
      n1 = size(lsao%StartLocalOrb,1)
   ENDIF
   call LS_MPI_BUFFER(n1,Master)
   IF(SLAVE)THEN
      call mem_LSTpointer_alloc(lsao%StartLocalOrb,n1)
   ENDIF
   call LS_MPI_BUFFER(lsao%StartLocalOrb,n1,Master)
ELSE
   IF(SLAVE)NULLIFY(lsao%StartLocalOrb)
ENDIF

isAssociated = ASSOCIATED(lsao%StartGlobalOrb)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      n1 = 0
   ELSE
      n1 = size(lsao%StartGlobalOrb,1)
   ENDIF
   call LS_MPI_BUFFER(n1,Master)
   IF(SLAVE)THEN
      call mem_LSTpointer_alloc(lsao%StartGlobalOrb,n1)
   ENDIF
   call LS_MPI_BUFFER(lsao%StartGlobalOrb,n1,Master)
ELSE
   IF(SLAVE)NULLIFY(lsao%StartGlobalOrb)
ENDIF

isAssociated = ASSOCIATED(lsao%nAngmom)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      n1 = 0
   ELSE
      n1 = size(lsao%nAngmom,1)
   ENDIF
   call LS_MPI_BUFFER(n1,Master)
   IF(SLAVE)THEN
      call mem_LSTpointer_alloc(lsao%nAngmom,n1)
   ENDIF
   call LS_MPI_BUFFER(lsao%nAngmom,n1,Master)
ELSE
   IF(SLAVE)NULLIFY(lsao%nAngmom)
ENDIF
end subroutine mpicopy_lsaotensor

SUBROUTINE mpicopy_output(output,Slave,Master)
implicit none
type(integralOutput) :: output
logical :: slave,isAssociated
integer(kind=ls_mpik) :: master
integer :: n1
call LS_MPI_BUFFER(output%ndim,5,Master)
!call LS_MPI_BUFFER(output%ndim3D,3,Master)
call LS_MPI_BUFFER(output%doGrad,Master)
call LS_MPI_BUFFER(output%USEBUFMM,Master)
call LS_MPI_BUFFER(output%MMBUFLEN,Master)
call LS_MPI_BUFFER(output%MAXBUFI,Master)
call LS_MPI_BUFFER(output%MAXBUFR,Master)
call LS_MPI_BUFFER(output%MAXBUFN,Master)
call LS_MPI_BUFFER(output%IBUFI,Master)
call LS_MPI_BUFFER(output%IBUFN,Master)
call LS_MPI_BUFFER(output%LUITNM,Master)
call LS_MPI_BUFFER(output%LUITNMR,Master)
call LS_MPI_BUFFER(output%decpacked,Master)
call LS_MPI_BUFFER(output%decpacked2,Master)
call LS_MPI_BUFFER(output%decpackedK,Master)
call LS_MPI_BUFFER(output%FullAlphaCD,Master)
call LS_MPI_BUFFER(output%exchangeFactor,Master)

isAssociated = ASSOCIATED(output%postprocess)
call LS_MPI_BUFFER(isAssociated,Master)
IF(isAssociated)THEN
   IF(SLAVE)THEN
      n1 = 0
   ELSE
      n1 = size(output%postprocess,1)
   ENDIF
   call LS_MPI_BUFFER(n1,Master)
   IF(SLAVE)THEN
      call mem_alloc(output%postprocess,n1)
   ENDIF
   call LS_MPI_BUFFER(output%postprocess,n1,Master)
ELSE
   IF(SLAVE)NULLIFY(output%postprocess)
ENDIF

END SUBROUTINE mpicopy_output

#ifdef VAR_MPI
SUBROUTINE mpicopy_integralconfig(dalton,Slave,Master)
implicit none
type(integralconfig) :: dalton
logical :: slave
integer(kind=ls_mpik) :: master

call LS_MPI_BUFFER(dalton%contang,Master)
call LS_MPI_BUFFER(dalton%NOGCINTEGRALTRANSFORM,Master)
call LS_MPI_BUFFER(dalton%FORCEGCBASIS,Master)
call LS_MPI_BUFFER(dalton%noOMP,Master)
call LS_MPI_BUFFER(dalton%IchorForceCPU,Master)
call LS_MPI_BUFFER(dalton%IchorForceGPU,Master)
call LS_MPI_BUFFER(dalton%UNRES,Master)
call LS_MPI_BUFFER(dalton%CFG_LSDALTON,Master)
call LS_MPI_BUFFER(dalton%TRILEVEL,Master)
call LS_MPI_BUFFER(dalton%DOPASS,Master)
call LS_MPI_BUFFER(dalton%DENSFIT,Master)
call LS_MPI_BUFFER(dalton%DF_K,Master)
call LS_MPI_BUFFER(dalton%INTEREST,Master)
call LS_MPI_BUFFER(dalton%LINSCA,Master)
call LS_MPI_BUFFER(dalton%MATRICESINMEMORY,Master)
call LS_MPI_BUFFER(dalton%MEMDIST,Master)
call LS_MPI_BUFFER(dalton%LOW_ACCURACY_START,Master)
call LS_MPI_BUFFER(dalton%LINSCAPRINT,Master)
call LS_MPI_BUFFER(dalton%AOPRINT,Master)
call LS_MPI_BUFFER(dalton%MOLPRINT,Master)
call LS_MPI_BUFFER(dalton%INTPRINT,Master)
call LS_MPI_BUFFER(dalton%BASPRINT,Master)
call LS_MPI_BUFFER(dalton%PRINTATOMCOORD,Master)
call LS_MPI_BUFFER(dalton%NOBQBQ,Master)
call LS_MPI_BUFFER(dalton%JENGINE,Master)
call LS_MPI_BUFFER(dalton%LOCALLINK,Master)
call LS_MPI_BUFFER(dalton%LOCALLINKmulthr,Master)
call LS_MPI_BUFFER(dalton%LOCALLINKDcont,Master)
call LS_MPI_BUFFER(dalton%LOCALLINKDthr,Master)
call LS_MPI_BUFFER(dalton%LOCALLINKsimmul,Master)
call LS_MPI_BUFFER(dalton%LOCALLINKoption,Master)
call LS_MPI_BUFFER(dalton%LOCALLINKincrem,Master)

call LS_MPI_BUFFER(dalton%FTUVmaxprim,Master)
call LS_MPI_BUFFER(dalton%maxpasses,Master)
call LS_MPI_BUFFER(dalton%FMM,Master)
call LS_MPI_BUFFER(dalton%LINK,Master)

call LS_MPI_BUFFER(dalton%LSDASCREEN,Master)
call LS_MPI_BUFFER(dalton%LSDAJENGINE,Master)
call LS_MPI_BUFFER(dalton%LSDACOULOMB,Master)
call LS_MPI_BUFFER(dalton%LSDALINK,Master)
call LS_MPI_BUFFER(dalton%LSDASCREEN_THRLOG,Master)
call LS_MPI_BUFFER(dalton%DAJENGINE,Master)
call LS_MPI_BUFFER(dalton%DACOULOMB,Master)
call LS_MPI_BUFFER(dalton%DALINK,Master)
call LS_MPI_BUFFER(dalton%DASCREEN_THRLOG,Master)

call LS_MPI_BUFFER(dalton%DEBUGOVERLAP,Master)
call LS_MPI_BUFFER(dalton%DEBUG4CENTER,Master)
call LS_MPI_BUFFER(dalton%DEBUG4CENTER_ERI,Master)
call LS_MPI_BUFFER(dalton%DEBUGPROP,Master)
call LS_MPI_BUFFER(dalton%DEBUGICHOR,Master)
call LS_MPI_BUFFER(dalton%DEBUGICHORLINK,Master)
call LS_MPI_BUFFER(dalton%DEBUGICHORLINKFULL,Master)
call LS_MPI_BUFFER(dalton%DEBUGICHOROPTION,Master)
call LS_MPI_BUFFER(dalton%DEBUGGEN1INT,Master)
call LS_MPI_BUFFER(dalton%DEBUGCGTODIFF,Master)
call LS_MPI_BUFFER(dalton%DEBUGEP,Master)
call LS_MPI_BUFFER(dalton%DEBUGscreen,Master)
call LS_MPI_BUFFER(dalton%DEBUGGEODERIVOVERLAP,Master)
call LS_MPI_BUFFER(dalton%DEBUGGEODERIVKINETIC,Master)
call LS_MPI_BUFFER(dalton%DEBUGGEODERIVEXCHANGE,Master)
call LS_MPI_BUFFER(dalton%DEBUGGEODERIVCOULOMB,Master)
call LS_MPI_BUFFER(dalton%DEBUGMAGDERIV,Master)
call LS_MPI_BUFFER(dalton%DEBUGMAGDERIVOVERLAP,Master)
call LS_MPI_BUFFER(dalton%DEBUGCCFRAGMENT,Master)
call LS_MPI_BUFFER(dalton%DEBUGKINETIC,Master)
call LS_MPI_BUFFER(dalton%DEBUGNUCPOT,Master)
call LS_MPI_BUFFER(dalton%DEBUGGGEM,Master)
call LS_MPI_BUFFER(dalton%DEBUGLSlib,Master)
call LS_MPI_BUFFER(dalton%DEBUGuncontAObatch,Master)
call LS_MPI_BUFFER(dalton%DEBUGDECPACKED,Master)

call LS_MPI_BUFFER(dalton%DO4CENTERERI,Master)
call LS_MPI_BUFFER(dalton%DUMP4CENTERERI,Master)
call LS_MPI_BUFFER(dalton%OVERLAP_DF_J,Master)
call LS_MPI_BUFFER(dalton%PARI_J,Master)
call LS_MPI_BUFFER(dalton%PARI_K,Master)
call LS_MPI_BUFFER(dalton%MOPARI_K,Master)
call LS_MPI_BUFFER(dalton%SIMPLE_PARI,Master)
call LS_MPI_BUFFER(dalton%NON_ROBUST_PARI,Master)

call LS_MPI_BUFFER(dalton%PARI_CHARGE,Master)
call LS_MPI_BUFFER(dalton%PARI_DIPOLE,Master)
call LS_MPI_BUFFER(dalton%TIMINGS,Master)
call LS_MPI_BUFFER(dalton%nonSphericalETUV,Master)
call LS_MPI_BUFFER(dalton%HIGH_RJ000_ACCURACY,Master)
!*FMM PARAMETERS
call LS_MPI_BUFFER(dalton%MM_LMAX,Master)
call LS_MPI_BUFFER(dalton%MM_TLMAX,Master)
call LS_MPI_BUFFER(dalton%MM_SCREEN,Master)
call LS_MPI_BUFFER(dalton%NO_MMFILES,Master)
call LS_MPI_BUFFER(dalton%MM_NO_ONE,Master)
call LS_MPI_BUFFER(dalton%CREATED_MMFILES,Master)
call LS_MPI_BUFFER(dalton%USEBUFMM,Master)
call LS_MPI_BUFFER(dalton%DO_MMGRD,Master)
call LS_MPI_BUFFER(dalton%MM_NOSCREEN,Master)
call LS_MPI_BUFFER(dalton%MMunique_ID1,Master)
!*BASIS PARAMETERS
call LS_MPI_BUFFER(dalton%ATOMBASIS,Master)
call LS_MPI_BUFFER(dalton%BASIS,nBasisBasParam,Master)
call LS_MPI_BUFFER(dalton%NOFAMILY,Master)
call LS_MPI_BUFFER(dalton%Hermiteecoeff,Master)
call LS_MPI_BUFFER(dalton%DoSpherical,Master)
call LS_MPI_BUFFER(dalton%UNCONT,Master) 
call LS_MPI_BUFFER(dalton%NOSEGMENT,Master)
!* JOB REQUESTS
call LS_MPI_BUFFER(dalton%DO3CENTEROVL,Master)
call LS_MPI_BUFFER(dalton%DO2CENTERERI,Master)
call LS_MPI_BUFFER(dalton%MIXEDOVERLAP,Master)

!*CAUCHY-SCHWARZ INTEGRAL PARAMETERS
!THE ONE THRESHOLD TO RULE THEM ALL
call LS_MPI_BUFFER(dalton%THRESHOLD,Master)
!THESE THRESHOLDS TELL HOW THEY SHOULD BE SET COMPARED TO THE ONE THRESHOLD
call LS_MPI_BUFFER(dalton%CS_THRESHOLD,Master)
call LS_MPI_BUFFER(dalton%OE_THRESHOLD,Master)
call LS_MPI_BUFFER(dalton%PS_THRESHOLD,Master)
call LS_MPI_BUFFER(dalton%OD_THRESHOLD,Master)
call LS_MPI_BUFFER(dalton%PARI_THRESHOLD,Master)
call LS_MPI_BUFFER(dalton%J_THR,Master)
call LS_MPI_BUFFER(dalton%K_THR,Master)
call LS_MPI_BUFFER(dalton%ONEEL_THR,Master)
!OTHER CAUCHY-SCHWARZ INTEGRAL PARAMETERS
call LS_MPI_BUFFER(dalton%CS_SCREEN,Master)
call LS_MPI_BUFFER(dalton%PARI_SCREEN,Master)
call LS_MPI_BUFFER(dalton%OE_SCREEN,Master)
call LS_MPI_BUFFER(dalton%savegabtomem,Master)
!PRIMITIVE INTEGRAL PARAMETERS
call LS_MPI_BUFFER(dalton%PS_SCREEN,Master)
call LS_MPI_BUFFER(dalton%PS_DEBUG,Master)
!Screen OD-batches by AO-batch extent
call LS_MPI_BUFFER(dalton%OD_SCREEN,Master)
call LS_MPI_BUFFER(dalton%MBIE_SCREEN,Master)
!Fragment molecule into to distinct parts, and construct matrices block by block
call LS_MPI_BUFFER(dalton%FRAGMENT,Master)
!Approximate number of atoms per fragment
call LS_MPI_BUFFER(dalton%numAtomsPerFragment,Master)
!FMM
call LS_MPI_BUFFER(dalton%LU_LUINTM,Master)
call LS_MPI_BUFFER(dalton%LU_LUINTR,Master)
call LS_MPI_BUFFER(dalton%LU_LUINDM,Master)
call LS_MPI_BUFFER(dalton%LU_LUINDR,Master)
call LS_MPI_BUFFER(dalton%LR_EXCHANGE_DF,Master)
call LS_MPI_BUFFER(dalton%LR_EXCHANGE_PARI,Master)
call LS_MPI_BUFFER(dalton%LR_EXCHANGE,Master)
call LS_MPI_BUFFER(dalton%ADMM_EXCHANGE,Master)
call LS_MPI_BUFFER(dalton%ADMM1,Master)
call LS_MPI_BUFFER(dalton%ADMM_2ERI,Master)
call LS_MPI_BUFFER(dalton%ADMMQ,Master)
call LS_MPI_BUFFER(dalton%ADMM_FUNC,len(dalton%ADMM_FUNC),Master)
call LS_MPI_BUFFER(dalton%ADMMS,Master)
call LS_MPI_BUFFER(dalton%ADMMP,Master)
call LS_MPI_BUFFER(dalton%ADMM_separateX,Master)
call LS_MPI_BUFFER(dalton%ADMMexchangeMetric,Master)
call LS_MPI_BUFFER(dalton%PRINT_EK3,Master)
call LS_MPI_BUFFER(dalton%ADMMBASISFILE,Master)
call LS_MPI_BUFFER(dalton%SR_EXCHANGE,Master)

!Coulomb attenuated method CAM parameters
call LS_MPI_BUFFER(dalton%CAM,Master)
call LS_MPI_BUFFER(dalton%CAMalpha,Master)
call LS_MPI_BUFFER(dalton%CAMbeta,Master)
call LS_MPI_BUFFER(dalton%CAMmu,Master)
!DFT PARAMETERS
call mpicopy_DFTparam(dalton%DFT,master)
call LS_MPI_BUFFER(dalton%exchangeFactor,Master)
call LS_MPI_BUFFER(dalton%nelectrons,Master)
call LS_MPI_BUFFER(dalton%molcharge,Master)
call LS_MPI_BUFFER(dalton%run_dec_gradient_test,Master)

call LS_MPI_BUFFER(dalton%ForceRIMP2memReduced,Master)
call LS_MPI_BUFFER(dalton%SolveNMRResponseSimultan,Master)
call LS_MPI_BUFFER(dalton%ResponseMatNormConvTest,Master)
call LS_MPI_BUFFER(dalton%PreCalcDFscreening,Master)
call LS_MPI_BUFFER(dalton%PreCalcF12screening,Master)

END SUBROUTINE MPICOPY_INTEGRALCONFIG
#endif

SUBROUTINE mpicopy_scheme(scheme,slave,master)
implicit none
type(lsintscheme) :: scheme
logical :: slave
integer(kind=ls_mpik) :: master

!PARAMETERS FROM **INTEGRALS   DECLERATION
call LS_MPI_BUFFER(scheme%NOBQBQ,Master)
call LS_MPI_BUFFER(scheme%doMPI,Master)
call LS_MPI_BUFFER(scheme%MasterWakeSlaves,Master)
call LS_MPI_BUFFER(scheme%noOMP,Master)
call LS_MPI_BUFFER(scheme%IchorForceCPU,Master)
call LS_MPI_BUFFER(scheme%IchorForceGPU,Master)
call LS_MPI_BUFFER(scheme%CFG_LSDALTON,Master)
call LS_MPI_BUFFER(scheme%DOPASS,Master)
call LS_MPI_BUFFER(scheme%DENSFIT,Master)
call LS_MPI_BUFFER(scheme%DF_K,Master)
call LS_MPI_BUFFER(scheme%INTEREST,Master)
call LS_MPI_BUFFER(scheme%MATRICESINMEMORY,Master)
call LS_MPI_BUFFER(scheme%MEMDIST,Master)
call LS_MPI_BUFFER(scheme%AOPRINT,Master)
call LS_MPI_BUFFER(scheme%INTPRINT,Master)
call LS_MPI_BUFFER(scheme%JENGINE,Master)
call LS_MPI_BUFFER(scheme%FTUVmaxprim,Master)
call LS_MPI_BUFFER(scheme%maxpasses,Master)
call LS_MPI_BUFFER(scheme%FMM,Master)
call LS_MPI_BUFFER(scheme%LINK,Master)
call LS_MPI_BUFFER(scheme%LSDASCREEN,Master)
call LS_MPI_BUFFER(scheme%LSDAJENGINE,Master)
call LS_MPI_BUFFER(scheme%LSDACOULOMB,Master)
call LS_MPI_BUFFER(scheme%LSDALINK,Master)
call LS_MPI_BUFFER(scheme%LSDASCREEN_THRLOG,Master)
call LS_MPI_BUFFER(scheme%DAJENGINE,Master)
call LS_MPI_BUFFER(scheme%DACOULOMB,Master)
call LS_MPI_BUFFER(scheme%DALINK,Master)
call LS_MPI_BUFFER(scheme%DASCREEN_THRLOG,Master)
call LS_MPI_BUFFER(scheme%DEBUGOVERLAP,Master)
call LS_MPI_BUFFER(scheme%DEBUG4CENTER,Master)
call LS_MPI_BUFFER(scheme%DEBUG4CENTER_ERI,Master)
call LS_MPI_BUFFER(scheme%DEBUGCCFRAGMENT,Master)
call LS_MPI_BUFFER(scheme%DEBUGKINETIC,Master)
call LS_MPI_BUFFER(scheme%DEBUGNUCPOT,Master)
call LS_MPI_BUFFER(scheme%DO4CENTERERI,Master)
call LS_MPI_BUFFER(scheme%OVERLAP_DF_J,Master)
call LS_MPI_BUFFER(scheme%PARI_J,Master)
call LS_MPI_BUFFER(scheme%PARI_K,Master)
call LS_MPI_BUFFER(scheme%MOPARI_K,Master)
call LS_MPI_BUFFER(scheme%SIMPLE_PARI,Master)
call LS_MPI_BUFFER(scheme%NON_ROBUST_PARI,Master)
call LS_MPI_BUFFER(scheme%PARI_CHARGE,Master)
call LS_MPI_BUFFER(scheme%PARI_DIPOLE,Master)
call LS_MPI_BUFFER(scheme%TIMINGS,Master)
call LS_MPI_BUFFER(scheme%nonSphericalETUV,Master)
call LS_MPI_BUFFER(scheme%HIGH_RJ000_ACCURACY,Master)
!*FMM PARAMETERS
call LS_MPI_BUFFER(scheme%MM_LMAX,Master)
call LS_MPI_BUFFER(scheme%MM_TLMAX,Master)
call LS_MPI_BUFFER(scheme%MM_SCREEN,Master)
call LS_MPI_BUFFER(scheme%NO_MMFILES,Master)
call LS_MPI_BUFFER(scheme%MM_NO_ONE,Master)
call LS_MPI_BUFFER(scheme%CREATED_MMFILES,Master)
call LS_MPI_BUFFER(scheme%USEBUFMM,Master)
call LS_MPI_BUFFER(scheme%DO_MMGRD,Master)
call LS_MPI_BUFFER(scheme%MM_NOSCREEN,Master)
call LS_MPI_BUFFER(scheme%MMunique_ID1,Master)
!FMM
call LS_MPI_BUFFER(scheme%LU_LUINTM,Master)
call LS_MPI_BUFFER(scheme%LU_LUINTR,Master)
call LS_MPI_BUFFER(scheme%LU_LUINDM,Master)
call LS_MPI_BUFFER(scheme%LU_LUINDR,Master)
!*BASIS PARAMETERS
call LS_MPI_BUFFER(scheme%BASIS,nBasisBasParam,Master)
call LS_MPI_BUFFER(scheme%NOFAMILY,Master)
call LS_MPI_BUFFER(scheme%Hermiteecoeff,Master)
call LS_MPI_BUFFER(scheme%DoSpherical,Master)
call LS_MPI_BUFFER(scheme%UNCONT,Master) !FORCE UNCONTRACTED BASIS
call LS_MPI_BUFFER(scheme%NOSEGMENT,Master) !DISABLE SEGMENTS 
call LS_MPI_BUFFER(scheme%CONTANG,Master) !Used contracted-angular AO ordering rather than angular-contracted
!* JOB REQUESTS
call LS_MPI_BUFFER(scheme%DO3CENTEROVL,Master)
call LS_MPI_BUFFER(scheme%DO2CENTERERI,Master)
call LS_MPI_BUFFER(scheme%CMORDER,Master)
call LS_MPI_BUFFER(scheme%CMIMAT,Master)
call LS_MPI_BUFFER(scheme%MIXEDOVERLAP,Master)
call LS_MPI_BUFFER(scheme%OD_MOM,Master)
call LS_MPI_BUFFER(scheme%MOM_CENTER,3,Master)
!*CAUCHY-SCHWARZ INTEGRAL PARAMETERS
!THE ONE THRESHOLD TO RULE THEM ALL
call LS_MPI_BUFFER(scheme%THRESHOLD,Master)
!THESE THRESHOLDS TELL HOW THEY SHOULD BE SET COMPARED TO THE ONE THRESHOLD
call LS_MPI_BUFFER(scheme%CS_THRESHOLD,Master)
call LS_MPI_BUFFER(scheme%OE_THRESHOLD,Master)
call LS_MPI_BUFFER(scheme%PS_THRESHOLD,Master)
call LS_MPI_BUFFER(scheme%OD_THRESHOLD,Master)
call LS_MPI_BUFFER(scheme%PARI_THRESHOLD,Master)
call LS_MPI_BUFFER(scheme%J_THR,Master)
call LS_MPI_BUFFER(scheme%K_THR,Master)
call LS_MPI_BUFFER(scheme%ONEEL_THR,Master)
call LS_MPI_BUFFER(scheme%IntThreshold,Master)
!OTHER CAUCHY-SCHWARZ INTEGRAL PARAMETERS
call LS_MPI_BUFFER(scheme%CS_SCREEN,Master)
call LS_MPI_BUFFER(scheme%PARI_SCREEN,Master)
call LS_MPI_BUFFER(scheme%OE_SCREEN,Master)
call LS_MPI_BUFFER(scheme%savegabtomem,Master)
call LS_MPI_BUFFER(scheme%ReCalcGab,Master)
call LS_MPI_BUFFER(scheme%CS_int,Master)
call LS_MPI_BUFFER(scheme%PS_int,Master)
!PRIMITIVE INTEGRAL PARAMETERS
call LS_MPI_BUFFER(scheme%PS_SCREEN,Master)
call LS_MPI_BUFFER(scheme%PS_DEBUG,Master)
!Screen OD-batches by AO-batch extent
call LS_MPI_BUFFER(scheme%OD_SCREEN,Master)
call LS_MPI_BUFFER(scheme%MBIE_SCREEN,Master)
!Fragment molecule into to distinct parts, and construct matrices block by block
call LS_MPI_BUFFER(scheme%FRAGMENT,Master)
!Approximate number of atoms per fragment
call LS_MPI_BUFFER(scheme%numAtomsPerFragment,Master)


!Coulomb attenuated method CAM parameters
call LS_MPI_BUFFER(scheme%LR_EXCHANGE_DF,Master)
call LS_MPI_BUFFER(scheme%LR_EXCHANGE_PARI,Master)
call LS_MPI_BUFFER(scheme%LR_EXCHANGE,Master)
call LS_MPI_BUFFER(scheme%SR_EXCHANGE,Master)

call LS_MPI_BUFFER(scheme%ADMM_EXCHANGE,Master)
call LS_MPI_BUFFER(scheme%ADMM1,Master)
call LS_MPI_BUFFER(scheme%ADMM_2ERI,Master)
call LS_MPI_BUFFER(scheme%ADMMQ,Master)
call LS_MPI_BUFFER(scheme%ADMMS,Master)
call LS_MPI_BUFFER(scheme%ADMMP,Master)
call LS_MPI_BUFFER(scheme%ADMM_separateX,Master)
call LS_MPI_BUFFER(scheme%ADMMexchangeMetric,Master)
call LS_MPI_BUFFER(scheme%PRINT_EK3,Master)
call LS_MPI_BUFFER(scheme%ADMM_CONSTRAIN_FACTOR,Master)
call LS_MPI_BUFFER(scheme%ADMM_LARGE_LAMBDA,Master)

call LS_MPI_BUFFER(scheme%CAM,Master)
call LS_MPI_BUFFER(scheme%CAMalpha,Master)
call LS_MPI_BUFFER(scheme%CAMbeta,Master)
call LS_MPI_BUFFER(scheme%CAMmu,Master)
call LS_MPI_BUFFER(scheme%exchangeFactor,Master)
!DFT PARAMETERS
call mpicopy_DFTparam(scheme%DFT,master)

call LS_MPI_BUFFER(scheme%INCREMENTAL,Master)
call LS_MPI_BUFFER(scheme%DO_PROP,Master)
call LS_MPI_BUFFER(scheme%PropOper,Master)

call LS_MPI_BUFFER(scheme%ForceRIMP2memReduced,Master)
call LS_MPI_BUFFER(scheme%AONuclearSpecID,Master)
call LS_MPI_BUFFER(scheme%PreCalcDFscreening,Master)
call LS_MPI_BUFFER(scheme%PreCalcF12screening,Master)

END SUBROUTINE mpicopy_scheme

#ifdef VAR_MPI
subroutine lsmpi_lstensor_reduction(Tensor,master,mynum,comm)
implicit none
type(lstensor) :: tensor
integer(kind=ls_mpik) :: master,mynum,comm

real(realk) :: t1,t2,t3,t4
CALL LS_GETTIM(t1,t2)
call ls_mpiInitBuffer(master,LSMPIREDUCTION,comm)
!plug lstensor into buffer
call lstensor_buffer(Tensor,master)
!reduction of buffer and free of slave buffer
CALL LS_GETTIM(t3,t4)
!IF(infpar%mynum.EQ.infpar%master)write(*,*) 'debug:master-reduction0',t3 - t1
CALL LS_GETTIM(t1,t2)
call ls_mpiFinalizeBuffer(master,LSMPIREDUCTION,comm)
IF(mynum.EQ.infpar%master)THEN
CALL LS_GETTIM(t3,t4)
!write(*,*) 'debug:master-reduction1',t3 - t1
CALL LS_GETTIM(t1,t2)
!build lstensor from buffer
call lstensor_buffer(Tensor,master)
!free the masters buffers
call ls_mpiFinalizeBuffer(master,LSMPIREDUCTIONmaster,comm)
CALL LS_GETTIM(t3,t4)
!write(*,*) 'debug:master-reduction2',t3 - t1
ENDIF
end subroutine lsmpi_lstensor_reduction
#endif

!> \brief fill the lstensor into buffer or extract from buffer
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR1 the original lstensor
SUBROUTINE lstensor_buffer(TENSOR1,master)
  implicit none
  TYPE(LSTENSOR)     :: TENSOR1
  integer(kind=ls_mpik) :: master
!
  INTEGER    :: I,nmat,dim,dim1,dim2,nMBIE
  !All logicals and integers should be the same for the lstensor on all nodes
  nmat = TENSOR1%ndim5
  IF(TENSOR1%gradienttensor)THEN
     dim = nmat*3
     DO I = 1,TENSOR1%nLSAO
        call LS_MPI_BUFFER(TENSOR1%LSAO(I)%elms(1:dim),dim,Master)
     ENDDO
  ELSE
     IF(ASSOCIATED(TENSOR1%LSAO))THEN
        DO I = 1,TENSOR1%nLSAO
           dim = nmat*TENSOR1%LSAO(I)%nelms
           IF(dim.GT.0)THEN
              call LS_MPI_BUFFER(TENSOR1%LSAO(I)%elms(1:dim),dim,Master)
           ENDIF
        ENDDO
     ENDIF
     IF(ASSOCIATED(TENSOR1%SLSAO))THEN
        DO I = 1,TENSOR1%nSLSAO
           dim = TENSOR1%SLSAO(I)%nelms
           IF(dim.GT.0)THEN
              call LS_MPI_BUFFER(TENSOR1%SLSAO(I)%Selms(1:dim),dim,Master)
           ENDIF
        ENDDO
     ENDIF
     IF(ASSOCIATED(TENSOR1%maxgab))THEN
        dim1 = TENSOR1%nbatches(1)
        dim2 = TENSOR1%nbatches(2)
        call LS_MPI_BUFFER(TENSOR1%maxgab,dim1,dim2,Master)
     ENDIF
     IF(ASSOCIATED(TENSOR1%maxprimgab))THEN
        dim1 = TENSOR1%nbatches(1)
        dim2 = TENSOR1%nbatches(2)
        call LS_MPI_BUFFER(TENSOR1%maxprimgab,dim1,dim2,Master)
     ENDIF
     IF(ASSOCIATED(TENSOR1%MBIE))THEN
        dim1 = TENSOR1%nbatches(1)
        dim2 = TENSOR1%nbatches(2)
        nMBIE = TENSOR1%nMBIE
        call LS_MPI_BUFFER(TENSOR1%MBIE,nMBIE,dim1,dim2,Master)
     ENDIF
  ENDIF
  call LS_MPI_BUFFER(TENSOR1%maxgabelm,Master)
  call LS_MPI_BUFFER(TENSOR1%maxprimgabelm,Master)
end SUBROUTINE lstensor_buffer

!> \brief MPI Copies(Broadcasts) a MOLECULE
!> \author T. Kjaergaard
!> \date 2010-05-31
!> \param Molecule Contains the information about the molecule
!> \param Slave if this processor is a slave 
!> \param Master the integer for the master process
subroutine mpicopy_molecule(Molecule,SLAVE,Master)
implicit none
TYPE(MOLECULEINFO),intent(INOUT) :: MOLECULE
logical :: Slave
integer(kind=ls_mpik) :: master
!
integer :: I

call LS_MPI_BUFFER(MOLECULE%nAtoms,Master)
call LS_MPI_BUFFER(MOLECULE%nAtomsNPC,Master)
call LS_MPI_BUFFER(MOLECULE%nelectrons,Master)
call LS_MPI_BUFFER(MOLECULE%charge,Master)

IF(SLAVE)THEN
   MOLECULE%nbastREG = 0
   MOLECULE%nbastAUX = 0
   MOLECULE%nbastCABS = 0
   MOLECULE%nbastJK = 0
   MOLECULE%nbastADMM = 0
   MOLECULE%nbastVAL = 0
   MOLECULE%nprimbastREG = 0
   MOLECULE%nprimbastAUX = 0
   MOLECULE%nprimbastCABS = 0
   MOLECULE%nprimbastJK = 0
   MOLECULE%nprimbastADMM = 0
   MOLECULE%nprimbastVAL = 0
   call mem_alloc(MOLECULE%ATOM,MOLECULE%nAtoms)
ENDIF
do I = 1, MOLECULE%nAtoms
   call mpicopy_atom(MOLECULE,I,Slave,Master)
enddo

call LS_MPI_BUFFER(MOLECULE%pointMolecule,Master)
call LS_MPI_BUFFER(MOLECULE%nSubSystems,Master)
call LS_MPI_BUFFER(MOLECULE%label,22,Master)
IF(Molecule%nSubSystems.NE.0)THEN
   IF(SLAVE)THEN
      call mem_alloc(Molecule%SubSystemLabel,Molecule%nSubSystems)
   ENDIF
   IF(len(Molecule%SubSystemLabel(1)).NE.80)THEN
      CALL LSQUIT('Dim mismatch in mpicopy_molecule',-1)
   ENDIF
   do I = 1,Molecule%nSubSystems       
      call LS_MPI_BUFFER(Molecule%SubSystemLabel(I),80,Master)
   enddo
ELSE
   IF(SLAVE)THEN
      NULLIFY(Molecule%SubSystemLabel)
   ENDIF
ENDIF

end subroutine mpicopy_molecule

!> \brief MPI Copies(Broadcasts) an atom from MOLECULE
!> \author S. Reine
!> \date 2010-02-21
!> \param Molecule Contains the information about the original molecule
!> \param I Atomic number of atom to be copied
!> \param Slave if this processor is a slave 
!> \param Master the integer for the master process
SUBROUTINE mpicopy_atom(MOLECULE,I,Slave,Master)
IMPLICIT NONE
INTEGER,intent(IN)               :: I
INTEGER(kind=ls_mpik),intent(IN) :: Master
TYPE(MOLECULEINFO),intent(INOUT)    :: MOLECULE
Logical :: Slave
!
integer :: K

call LS_MPI_BUFFER(MOLECULE%ATOM(I)%Isotope,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%Name,len(MOLECULE%ATOM(I)%Name),Master)

call LS_MPI_BUFFER(MOLECULE%ATOM(I)%MASS,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%CovRad,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%Frag,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%CENTER,3,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%Atomic_number,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%molecularIndex,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%SubsystemIndex,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%Charge,Master)
!call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nbasis,Master)
do K = 1,nBasisBasParam!MOLECULE%ATOM(I)%nbasis
   call LS_MPI_BUFFER(MOLECULE%ATOM(I)%basislabel(K),len(MOLECULE%ATOM(I)%basislabel(K)),Master)
   call LS_MPI_BUFFER(MOLECULE%ATOM(I)%basisindex(K),Master)
   call LS_MPI_BUFFER(MOLECULE%ATOM(I)%IDtype(K),Master)
enddo
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%phantom,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%pointcharge,Master)

call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nContOrbREG,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nPrimOrbREG,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nContOrbAUX,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nPrimOrbAUX,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nContOrbCABS,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nPrimOrbCABS,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nContOrbJK,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nPrimOrbJK,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nContOrbADMM,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nPrimOrbADMM,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nContOrbVAL,Master)
call LS_MPI_BUFFER(MOLECULE%ATOM(I)%nPrimOrbVAL,Master)
IF(SLAVE)THEN
   MOLECULE%nbastREG     = MOLECULE%nbastREG     + MOLECULE%ATOM(I)%nContOrbREG
   MOLECULE%nPrimbastREG = MOLECULE%nPrimbastREG + MOLECULE%ATOM(I)%nPrimOrbREG
   MOLECULE%nbastAUX     = MOLECULE%nbastAUX     + MOLECULE%ATOM(I)%nContOrbAUX
   MOLECULE%nPrimbastAUX = MOLECULE%nPrimbastAUX + MOLECULE%ATOM(I)%nPrimOrbAUX
   MOLECULE%nbastCABS     = MOLECULE%nbastCABS   + MOLECULE%ATOM(I)%nContOrbCABS
   MOLECULE%nPrimbastCABS= MOLECULE%nPrimbastCABS+ MOLECULE%ATOM(I)%nPrimOrbCABS
   MOLECULE%nbastJK      = MOLECULE%nbastJK     + MOLECULE%ATOM(I)%nContOrbJK
   MOLECULE%nPrimbastJK = MOLECULE%nPrimbastJK + MOLECULE%ATOM(I)%nPrimOrbJK
   MOLECULE%nbastADMM     = MOLECULE%nbastADMM   + MOLECULE%ATOM(I)%nContOrbADMM
   MOLECULE%nPrimbastADMM= MOLECULE%nPrimbastADMM+ MOLECULE%ATOM(I)%nPrimOrbADMM
   MOLECULE%nbastVAL     = MOLECULE%nbastVAL     + MOLECULE%ATOM(I)%nContOrbVAL
   MOLECULE%nPrimbastVAL = MOLECULE%nPrimbastVAL + MOLECULE%ATOM(I)%nPrimOrbVAL
ENDIF
END SUBROUTINE mpicopy_atom

!> \brief MPI braodcasts basissetinfo
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param bas is the BASISSETINFO
!> \param slave logical to tell if processor is slave
!> \param master integer of master process
SUBROUTINE mpicopy_basissetinfo(BAS,slave,master)
  implicit none
  TYPE(BASISSETINFO),intent(inout)    :: BAS
  Integer(kind=ls_mpik)               :: master
  logical            :: slave
!
  INTEGER            :: I,J,K,nrow,ncol,nsize1,nsize2

  IF(slave)THEN
     call nullifyBasisset(BAS)
  ENDIF
  call LS_MPI_BUFFER(BAS%natomtypes,Master)
  call LS_MPI_BUFFER(BAS%label,len(BAS%label),Master)
  call LS_MPI_BUFFER(BAS%labelindex,Master)
  call LS_MPI_BUFFER(BAS%nChargeindex,Master)
  call LS_MPI_BUFFER(BAS%nbast,Master)
  call LS_MPI_BUFFER(BAS%nprimbast,Master)
  call LS_MPI_BUFFER(BAS%DunningsBasis,Master)
  call LS_MPI_BUFFER(BAS%GeminalScalingFactor,Master)
  call LS_MPI_BUFFER(BAS%GCbasis,Master)
  call LS_MPI_BUFFER(BAS%Spherical,Master)
  call LS_MPI_BUFFER(BAS%Gcont,Master)     
  IF(BAS%natomtypes.NE. 0)THEN
     IF(slave)THEN
        call mem_alloc(BAS%ATOMTYPE,BAS%natomtypes)
     ENDIF
     DO I = 1,BAS%natomtypes
        call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%nAngmom,Master)
        call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%family,Master)
        call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%ToTnorb,Master)
        call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%ToTnprim,Master)
        call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%Charge,Master)
        call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%NAME,LEN(BAS%ATOMTYPE(I)%NAME),Master)
        DO J=1,BAS%ATOMTYPE(I)%nAngmom
           call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%SHELL(J)%nprim,Master)
           call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%SHELL(J)%norb,Master)
           call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%SHELL(J)%nsegments,Master)
           DO K=1,BAS%ATOMTYPE(I)%SHELL(J)%nsegments
              call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow,Master)
              call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol,Master)
              !in the case of trilevel the elms maybe be allocd larger 
              !and then some contracted are removed 
              !(without dealloc and reallocation) so it is safest to 
              !use nsize instead of nrow,ncol.
              IF(slave)THEN
                 nsize1 = 0
              ELSE
                 nsize1 = SIZE(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms)
              ENDIF
              call LS_MPI_BUFFER(nsize1,Master)
              IF(SLAVE)THEN
                 call mem_alloc(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms,nsize1)
              ENDIF
              call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms,nsize1,Master)
              IF(slave)THEN
                 nsize1 = 0
              ELSE
                 nsize1 = SIZE(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms)
              ENDIF
              call LS_MPI_BUFFER(nsize1,Master)
              IF(SLAVE)THEN
                 call mem_alloc(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms,nsize1)
              ENDIF
              call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms,nsize1,Master)
              IF(slave)THEN
                 nsize2 = 0
              ELSE
                 nsize2=SIZE(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents)
              ENDIF
              call LS_MPI_BUFFER(nsize2,Master)
              IF(SLAVE)THEN
                 call mem_alloc(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents,nsize2)
              ENDIF
              call LS_MPI_BUFFER(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents,nsize2,Master)
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  IF(BAS%nChargeindex .NE. 0)THEN
     IF(SLAVE)THEN
        call mem_alloc(BAS%Chargeindex,BAS%nChargeindex,.TRUE.)
     ENDIF
     DO I = 0,BAS%nChargeindex
        call LS_MPI_BUFFER(BAS%Chargeindex(I),Master)
     ENDDO
  ENDIF

END SUBROUTINE mpicopy_basissetinfo

SUBROUTINE mpicopy_GaussianGeminal(GGem,slave,master)
implicit none
TYPE(GaussianGeminal),intent(inout) :: GGem
Integer(kind=ls_mpik)               :: master
logical                             :: slave
!
INTEGER            :: I,J,K,nrow,ncol,nsize1,nsize2
LOGICAL  :: isassociated

call LS_MPI_BUFFER(GGem%is_set,Master)
IF (GGem%is_set) THEN
  call LS_MPI_BUFFER(GGem%N,Master)
  IF (GGem%N.GT.0) THEN
    call LS_MPI_BUFFER(GGem%coeff,GGem%N,Master)
    call LS_MPI_BUFFER(GGem%exponent,GGem%N,Master)
    call LS_MPI_BUFFER(GGem%expProd,GGem%N,Master)
  ENDIF
ENDIF
END SUBROUTINE mpicopy_GaussianGeminal

SUBROUTINE mpicopy_reduced_screen_info(redCS,slave,master)
implicit none
type(reducedScreeningInfo),intent(INOUT) :: redCS
Integer(kind=ls_mpik)                    :: master
logical                                  :: slave

!Does not copy. The setup of redCS is done later for the slave
IF (SLAVE) call init_reduced_screen_info(redCS)

END SUBROUTINE mpicopy_reduced_screen_info

#ifdef VAR_MPI
  subroutine init_slave_timers(times,comm)
     implicit none
     integer(kind=ls_mpik),intent(in) :: comm
     real(realk), pointer,intent(inout) :: times(:)
     integer(kind=ls_mpik) :: nnod,me

     call time_start_phase(PHASE_WORK)
     if(associated(times)) call lsquit("ERROR(init_slave_timers):pointer already associtated",-1)

     call get_rank_for_comm( comm, me   )
     call get_size_for_comm( comm, nnod )


     if(me==0_ls_mpik) then
        call time_start_phase(PHASE_COMM)
        call ls_mpibcast(INITSLAVETIME,me,comm)
        call time_start_phase(PHASE_WORK)
     endif

     call mem_alloc(times,nphases*nnod)
     times = 0.0E0_realk
     call time_start_phase(PHASE_WORK, &
        &swinit = times(me*nphases+PHASE_INIT_IDX) ,&
        &swwork = times(me*nphases+PHASE_WORK_IDX) ,&
        &swcomm = times(me*nphases+PHASE_COMM_IDX) ,&
        &swidle = times(me*nphases+PHASE_IDLE_IDX) )


     call time_start_phase(PHASE_COMM)
     call lsmpi_reduction(times,nphases*nnod,infpar%master,comm)
     call time_start_phase(PHASE_WORK)


  end subroutine init_slave_timers

  subroutine get_slave_timers(times,comm)
     implicit none
     integer(kind=ls_mpik),intent(in) :: comm
     real(realk), pointer,intent(inout) :: times(:)
     integer(kind=ls_mpik) :: nnod,me

     call time_start_phase(PHASE_WORK)
     if(.not.associated(times)) call lsquit("ERROR(get_slave_timers):pointer not associtated",-1)

     call get_rank_for_comm( comm, me   )
     call get_size_for_comm( comm, nnod )


     if(me==infpar%master) then
        call time_start_phase(PHASE_COMM)
        call ls_mpibcast(GETSLAVETIME,me,comm)
        call time_start_phase(PHASE_WORK)
        times(nphases+1:) = -1.E0_realk * times(nphases+1:)
     endif

     call time_start_phase(PHASE_WORK, &
        &dwinit = times(me*nphases+PHASE_INIT_IDX) ,&
        &dwwork = times(me*nphases+PHASE_WORK_IDX) ,&
        &dwcomm = times(me*nphases+PHASE_COMM_IDX) ,&
        &dwidle = times(me*nphases+PHASE_IDLE_IDX) )


     call time_start_phase(PHASE_COMM)
     call lsmpi_reduction(times,nphases*nnod,infpar%master,comm)
     call time_start_phase(PHASE_WORK)


  end subroutine get_slave_timers

  subroutine mem_init_background_alloc_all_nodes(comm,bytes)
     implicit none
     integer(kind=8),intent(in) :: bytes
     integer(kind=ls_mpik),intent(in) :: comm
     integer(kind=ls_mpik) :: nnod,me
     integer(kind=8) :: bytes_int

     call time_start_phase(PHASE_WORK)
     
     bytes_int  = bytes

     call get_rank_for_comm( comm, me   )
     call get_size_for_comm( comm, nnod )


     call time_start_phase(PHASE_COMM)
     if(me==infpar%master) then
        call ls_mpibcast(INIT_BG_BUF,infpar%master,comm)
     endif
     call ls_mpibcast(bytes_int,infpar%master,comm)
     call time_start_phase(PHASE_WORK)

     call mem_init_background_alloc(bytes_int)

  end subroutine mem_init_background_alloc_all_nodes
  subroutine mem_free_background_alloc_all_nodes(comm)
     implicit none
     integer(kind=ls_mpik),intent(in) :: comm
     integer(kind=ls_mpik) :: nnod,me
     real(realk) :: bytes_int

     call time_start_phase(PHASE_WORK)
     
     call get_rank_for_comm( comm, me   )
     call get_size_for_comm( comm, nnod )

     call time_start_phase(PHASE_COMM)
     if(me==infpar%master) then
        call ls_mpibcast(FREE_BG_BUF,me,comm)
     endif
     call time_start_phase(PHASE_WORK)

     call mem_free_background_alloc()

  end subroutine mem_free_background_alloc_all_nodes

  subroutine mem_change_background_alloc_all_nodes(comm,bytes)
     implicit none
     integer(kind=8),intent(in) :: bytes
     integer(kind=ls_mpik),intent(in) :: comm
     integer(kind=ls_mpik) :: nnod,me
     integer(kind=8) :: bytes_int
     call time_start_phase(PHASE_WORK)
     
     bytes_int  = bytes

     call get_rank_for_comm( comm, me   )
     call get_size_for_comm( comm, nnod )


     call time_start_phase(PHASE_COMM)
     if(me==infpar%master) then
        call ls_mpibcast(CHANGE_BG_BUF,infpar%master,comm)
     endif
     call ls_mpibcast(bytes_int,infpar%master,comm)
     call time_start_phase(PHASE_WORK)

     call mem_change_background_alloc(bytes_int)

   end subroutine mem_change_background_alloc_all_nodes
#endif

end module lsmpi_op

#ifdef VAR_MPI
subroutine init_slave_timers_slave(comm)
   use precision, only: realk
   use lstiming
   use memory_handling, only: mem_dealloc
   use lsmpi_op, only: init_slave_timers
   implicit none
   integer(kind=ls_mpik) :: comm,nnod
   real(realk), pointer :: times(:)

   times => null()

   call init_slave_timers(times,comm)

   call mem_dealloc(times)

end subroutine init_slave_timers_slave

subroutine get_slave_timers_slave(comm)
   use precision, only: realk, ls_mpik
   use lstiming, only: nphases
   use lsmpi_type, only: get_size_for_comm
   use memory_handling, only: mem_alloc,mem_dealloc
   use lsmpi_op, only: get_slave_timers
   implicit none
   integer(kind=ls_mpik) :: comm,nnod
   real(realk), pointer :: times(:)

   call get_size_for_comm(comm, nnod)
   call mem_alloc(times,nphases*nnod) 
   times = 0.0E0_realk

   call get_slave_timers(times,comm)

   call mem_dealloc(times)

end subroutine get_slave_timers_slave

subroutine mem_init_background_alloc_slave(comm)
   use precision, only: realk, ls_mpik
   use lsmpi_op, only: mem_init_background_alloc_all_nodes
   implicit none
   integer(kind=ls_mpik),intent(in) :: comm
   integer(kind=8):: bytes

   bytes=8
   call mem_init_background_alloc_all_nodes(comm,bytes)

end subroutine mem_init_background_alloc_slave

subroutine mem_free_background_alloc_slave(comm)
   use precision, only: realk, ls_mpik
   use lsmpi_op, only: mem_free_background_alloc_all_nodes
   implicit none
   integer(kind=ls_mpik),intent(in) :: comm

   call mem_free_background_alloc_all_nodes(comm)

end subroutine mem_free_background_alloc_slave
#endif
