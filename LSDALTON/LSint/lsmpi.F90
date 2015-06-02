module lsmpi_mod
  use precision
  use lsmpi_type
  use lsmpi_op
  use LSparameters
  use memory_handling
  use typedefType
  use typedef
  use molecule_module
#ifdef VAR_MPI
  use infpar_module
#endif
contains
    SUBROUTINE lsmpi_set_task_manager(tasks,Routine,AO1,AO2,AO3,AO4,Spec,intType,setting,BothPartitioning,lupri,luerr)
    implicit none
    Integer               :: LUPRI,LUERR
    Type(LSSETTING)       :: SETTING
    integer               :: AO1,AO2,AO3,AO4,intType,Spec
    Character*(*)         :: Routine
    type(ls_task_manager) :: tasks
    logical               :: BothPartitioning
!
    logical             :: noA,noB,noC,noD
    integer             :: nAtomsA,nAtomsB,nAtomsC,nAtomsD
    integer             :: nAtomsA_full,nAtomsB_full,nAtomsC_full,nAtomsD_full
    integer             :: iAO,natomsLHS,natomsRHS
    real(realk),pointer :: timeMatrixLHS(:,:)
    real(realk),pointer :: timeMatrixRHS(:,:)
    logical             :: emptyMol(4),no(4)
    integer             :: target_ntasks
    integer             :: nAtoms(4),nAtoms_full(4)

    real(realk) :: t1,t2,t3,t4

#ifdef VAR_MPI
CALL LS_GETTIM(t1,t2)
    DO iAO=1,4 
       ! Consistency testing
       IF(.NOT.associated(setting%molecule(iAO)%p)) THEN
          emptyMol(iAO) = .TRUE.
          ALLOCATE(setting%molecule(iAO)%p)
          call build_empty_molecule(setting%molecule(iAO)%p)
       ELSE
          emptyMol(iAO) = .FALSE.
       ENDIF
       ! Set up the information about the orbitals
       call setMolecularOrbitalInfo(setting%molecule(iAO)%p,tasks%orbInfo(iAO))
    ENDDO

    nAtomsA_full = setting%molecule(1)%p%nAtomsNPC
    nAtomsB_full = setting%molecule(2)%p%nAtomsNPC
    nAtomsC_full = setting%molecule(3)%p%nAtomsNPC
    nAtomsD_full = setting%molecule(4)%p%nAtomsNPC
    IF (AO1.EQ.AOEmpty) nAtomsA_full = 1
    IF (AO2.EQ.AOEmpty) nAtomsB_full = 1
    IF (AO3.EQ.AOEmpty) nAtomsC_full = 1
    IF (AO4.EQ.AOEmpty) nAtomsD_full = 1

    IF (AO1.EQ.AONuclear) nAtomsA_full = setting%molecule(1)%p%nAtoms
    IF (AO2.EQ.AONuclear) nAtomsB_full = setting%molecule(2)%p%nAtoms
    IF (AO3.EQ.AONuclear) nAtomsC_full = setting%molecule(3)%p%nAtoms
    IF (AO4.EQ.AONuclear) nAtomsD_full = setting%molecule(4)%p%nAtoms

    IF (AO1.EQ.AOPcharge) nAtomsA_full = setting%molecule(1)%p%nAtoms
    IF (AO2.EQ.AOPcharge) nAtomsB_full = setting%molecule(2)%p%nAtoms
    IF (AO3.EQ.AOPcharge) nAtomsC_full = setting%molecule(3)%p%nAtoms
    IF (AO4.EQ.AOPcharge) nAtomsD_full = setting%molecule(4)%p%nAtoms
      
    IF (AO1.EQ.AOelField) nAtomsA_full = setting%molecule(1)%p%nAtoms
    IF (AO2.EQ.AOelField) nAtomsB_full = setting%molecule(2)%p%nAtoms
    IF (AO3.EQ.AOelField) nAtomsC_full = setting%molecule(3)%p%nAtoms
    IF (AO4.EQ.AOelField) nAtomsD_full = setting%molecule(4)%p%nAtoms

    IF (AO1.EQ.AOdfCABS) nAtomsA_full = 2*setting%molecule(1)%p%nAtoms
    IF (AO2.EQ.AOdfCABS) nAtomsB_full = 2*setting%molecule(2)%p%nAtoms
    IF (AO3.EQ.AOdfCABS) nAtomsC_full = 2*setting%molecule(3)%p%nAtoms
    IF (AO4.EQ.AOdfCABS) nAtomsD_full = 2*setting%molecule(4)%p%nAtoms

      
    !Specify if auxiliary basis sets are used for left- or right-hand side, used when building fragments
    tasks%lhs_aux = AO1.EQ.AOdfdefault.OR.AO2.EQ.AOdfdefault
    tasks%rhs_aux = AO3.EQ.AOdfdefault.OR.AO4.EQ.AOdfdefault

    CALL GetSymmetries(tasks%sameAOsLHS,tasks%sameAOsRHS,tasks%sameODs,AO1,AO2,AO3,AO4,Setting%sameMol,Setting%nAO)

    natomsLHS = getNAtomsSide(AO1,AO2,tasks%sameAOsLHS,setting%molecule(1)%p,setting%molecule(2)%p)
    natomsRHS = getNAtomsSide(AO3,AO4,tasks%sameAOsRHS,setting%molecule(3)%p,setting%molecule(4)%p)

    target_ntasks = setting%numnodes
    IF (setting%numnodes.EQ.1) target_ntasks = 1

    IF(BothPartitioning)THEN
       IF ((AO3.EQ.AORdefault).AND.(AO4.EQ.AORdefault)) THEN
          IF ((AO1.EQ.AORdefault).AND.(AO2.EQ.AORdefault)) THEN
             tasks%target_nrhs = min(target_ntasks,natomsRHS)
             tasks%target_nlhs = min(target_ntasks,natomsLHS)
          ELSE
             call lsquit('error1 in use of BothPartitioning lsmpi_set_task_m',-1)
          ENDIF
       ELSE
          call lsquit('error2 in use of BothPartitioning lsmpi_set_task_m',-1)
       ENDIF
    ELSEIF(Setting%scheme%cs_int)THEN
       tasks%target_nlhs = min(target_ntasks,natomsLHS)
       tasks%target_nrhs = 1
    ELSEIF ((AO3.EQ.AORdefault).AND.(AO4.EQ.AORdefault)) THEN
      tasks%target_nrhs = min(target_ntasks,natomsRHS)
      tasks%target_nlhs = 1
      IF (natomsRHS.LT.target_ntasks) tasks%target_nlhs = natomsLHS
      IF (natomsLHS.EQ.0) tasks%target_nlhs=1
    ELSE IF ((AO1.EQ.AORdefault).AND.(AO2.EQ.AORdefault)) THEN
      tasks%target_nlhs = min(target_ntasks,natomsLHS)
      tasks%target_nrhs = 1
      IF (natomsLHS.LT.target_ntasks) tasks%target_nrhs = natomsRHS
    ELSE
      tasks%target_nlhs = 1
      IF (natomsRHS.GT.target_ntasks) &
    &      tasks%target_nlhs = min(target_ntasks,natomsLHS)
      tasks%target_nrhs = min(target_ntasks,natomsRHS)
    ENDIF
    ! Set up LHS and RHS time-estimates
    IF (ROUTINE.EQ.'LinK') THEN
       call LinK_time_estimate(timeMatrixLHS,timeMatrixRHS,nAtomsA,nAtomsB,nAtomsC,nAtomsD,&
     &                         tasks%sameODs,noA,noB,noC,noD,AO1,AO2,AO3,AO4,&
     &                         Spec,setting,lupri)
       nAtoms(1) = nAtomsA
       nAtoms(2) = nAtomsC
       nAtoms(3) = nAtomsB
       nAtoms(4) = nAtomsD
       nAtoms_full(1) = nAtomsA_full
       nAtoms_full(2) = nAtomsC_full
       nAtoms_full(3) = nAtomsB_full
       nAtoms_full(4) = nAtomsD_full
       no(1)     = noA
       no(2)     = noC
       no(3)     = noB
       no(4)     = noD
    ELSE IF (ROUTINE.EQ.'JENGINE') THEN
       call jengine_time_estimate(timeMatrixLHS,timeMatrixRHS,nAtomsA,nAtomsB,nAtomsC,nAtomsD,&
     &                            tasks%sameAOsLHS,tasks%sameAOsRHS,noA,noB,noC,noD,AO1,AO2,AO3,AO4,&
     &                            Spec,setting,lupri)
       nAtoms(1) = nAtomsA
       nAtoms(2) = nAtomsB
       nAtoms(3) = nAtomsC
       nAtoms(4) = nAtomsD
       nAtoms_full(1) = nAtomsA_full
       nAtoms_full(2) = nAtomsB_full
       nAtoms_full(3) = nAtomsC_full
       nAtoms_full(4) = nAtomsD_full
       no(1)     = noA
       no(2)     = noB
       no(3)     = noC
       no(4)     = noD
    ELSE
       call lsquit('Programming error in lsmpi_set_task_manager. Wrong routine specification.',-1)
    ENDIF

CALL LS_GETTIM(t3,t4)
!IF (setting%node.EQ.infpar%master) write(lupri,*) 'debug:estimate',t3-t1
CALL LS_GETTIM(t1,t2)

    ! Set up LHS tasks
    call ls_init_tasks(tasks%lhs,tasks%target_nlhs,nAtoms(1),nAtoms(2),nAtoms_full(1),nAtoms_full(2),no(1),no(2))
    call ls_time_tasks(tasks%lhs,timeMatrixLHS,nAtoms(1),nAtoms(2),tasks%target_nlhs,&
     &                 lupri,setting%scheme%intprint)
CALL LS_GETTIM(t3,t4)
!IF (setting%node.EQ.infpar%master) write(lupri,*) 'debug:task1',t3-t1
CALL LS_GETTIM(t1,t2)
    IF (setting%scheme%intprint.GT.0) THEN
      if (setting%node.EQ.infpar%master) call ls_print_tasks(tasks%lhs,lupri)
    ENDIF
    call mem_dealloc(timeMatrixLHS)
    
    ! Set up RHS tasks
    call ls_init_tasks(tasks%rhs,tasks%target_nrhs,nAtoms(3),nAtoms(4),nAtoms_full(3),nAtoms_full(4),no(3),no(4))
    call ls_time_tasks(tasks%rhs,timeMatrixRHS,nAtoms(3),nAtoms(4),tasks%target_nrhs,&
     &                 lupri,setting%scheme%intprint)
    call mem_dealloc(timeMatrixRHS)
    IF (setting%scheme%intprint.GT.0) THEN
      if (setting%node.EQ.infpar%master) call ls_print_tasks(tasks%rhs,lupri)
    ENDIF

CALL LS_GETTIM(t3,t4)
!IF (setting%node.EQ.infpar%master) write(lupri,*) 'debug:task2',t3-t1
CALL LS_GETTIM(t1,t2)
    tasks%itask  = -1
    tasks%lhs_current => tasks%lhs%first
    tasks%rhs_current => tasks%rhs%first
    tasks%ntasks = tasks%lhs%ntasks * tasks%rhs%ntasks

!   Free redCS - set up in ls_create_lstensor_full and used to set up the time matrices above
    call free_reduced_screen_info(setting%redCS)

    DO iAO=1,4 
      IF (emptyMol(iAO)) THEN
        call free_moleculeinfo(setting%molecule(iAO)%p)
        DEALLOCATE(setting%molecule(iAO)%p)
        NULLIFY(setting%molecule(iAO)%p)
      ENDIF
    ENDDO
#endif

    END SUBROUTINE lsmpi_set_task_manager

#ifdef VAR_MPI
SUBROUTINE lsmpi_set_MPI_task_list(task_list,tasks,setting,lupri,luerr)
implicit none
TYPE(LSMPI_TASK_LIST),INTENT(INOUT) :: task_list
TYPE(LS_TASK_MANAGER),INTENT(IN) :: tasks
TYPE(LSSETTING),INTENT(IN)       :: setting
INTEGER,INTENT(IN)               :: lupri,luerr
!
TYPE(LSTASK),pointer       :: lhs_current,rhs_current
TYPE(LSTASK_p),allocatable :: lhs_list(:),rhs_list(:)
Integer,pointer :: node_lhs(:),node_rhs(:)
Integer :: node,i,j,itask
logical :: noLHS,noRHS

Integer :: offnode

!Offnode determine which node will get which task
offnode = 0

task_list%mynode = setting%node
IF (tasks%rhs%ntasks.GT.1) THEN
  task_list%ownRHS = .TRUE.
  task_list%ownLHS = .FALSE.
ELSE IF (tasks%lhs%ntasks.GT.1) THEN
  task_list%ownRHS = .FALSE.
  task_list%ownLHS = .TRUE.
ELSE
  task_list%ownRHS = .TRUE.
  task_list%ownLHS = .TRUE.
ENDIF

!Set up LHS information
noLHS = tasks%lhs%ntasks.EQ.1
allocate(lhs_list(tasks%lhs%ntasks))
call mem_alloc(node_lhs,tasks%lhs%ntasks)
IF (task_list%ownlHS) THEN
  task_list%nuniquelHS = 0
  lhs_current => tasks%lhs%first
  DO I=1,tasks%lhs%ntasks
     node = mod(i-offnode,setting%Numnodes)
     IF (node.EQ.setting%node)THEN
       task_list%nuniquelHS = task_list%nuniquelHS + 1
       lhs_list(task_list%nuniquelHS)%p => lhs_current
       node_lhs(task_list%nuniquelHS) = node
     ENDIF
     lhs_current => lhs_current%next
  ENDDO
ELSE
  task_list%nuniquelHS = tasks%lhs%ntasks
  lhs_current => tasks%lhs%first
  DO I=1,tasks%lhs%ntasks
     lhs_list(I)%p => lhs_current
     lhs_current => lhs_current%next
     node = mod(i-offnode,setting%Numnodes)
     IF (noLHS) node = -1 !This means no partitioning (and therefore a reduction of the result)
     node_lhs(I) = node
  ENDDO
ENDIF

!Set up RHS information
noRHS = tasks%rhs%ntasks.EQ.1
allocate(rhs_list(tasks%rhs%ntasks))
call mem_alloc(node_rhs,tasks%rhs%ntasks)
IF (task_list%ownRHS) THEN
  task_list%nuniqueRHS = 0
  rhs_current => tasks%rhs%first
  DO I=1,tasks%rhs%ntasks
     node = mod(i-offnode,setting%Numnodes)
     IF (node.EQ.setting%node)THEN
       task_list%nuniqueRHS = task_list%nuniqueRHS + 1
       rhs_list(task_list%nuniqueRHS)%p => rhs_current
       node_rhs(task_list%nuniqueRHS) = node
     ENDIF
     rhs_current => rhs_current%next
  ENDDO
ELSE
  task_list%nuniqueRHS = tasks%rhs%ntasks
  rhs_current => tasks%rhs%first
  DO I=1,tasks%rhs%ntasks
     rhs_list(I)%p => rhs_current
     rhs_current => rhs_current%next
     node = mod(i-offnode,setting%Numnodes)
     IF (noRHS) node = -1 ! This means no partitioning (and a syncronization of RHS Dmat to all nodes)
     node_rhs(I) = node
  ENDDO
ENDIF

!Set up the MPI task-list
task_list%numMPItasks = task_list%nuniqueLHS*task_list%nuniqueRHS
CALL mem_alloc(task_list%taskPair,2,task_list%numMPItasks)
CALL mem_alloc(task_list%taskNode,2,task_list%numMPItasks)
ALLOCATE(task_list%LHS(task_list%numMPItasks))
ALLOCATE(task_list%RHS(task_list%numMPItasks))
!CALL mem_alloc(task_list%taskWall,task_list%numMPItasks)
!CALL mem_alloc(task_list%taskCPU,task_list%numMPItasks)
!CALL mem_alloc(task_list%taskEstimate,task_list%numMPItasks)
itask=0
DO I=1,task_list%nuniqueLHS
  DO J=1,task_list%nuniqueRHS
    itask=itask+1
    task_list%taskPair(1,itask) = I
    task_list%taskPair(2,itask) = J
    task_list%taskNode(1,itask) = node_lhs(I)
    task_list%taskNode(2,itask) = node_rhs(J)
    task_list%LHS(itask)%p => lhs_list(I)%p 
    task_list%RHS(itask)%p => rhs_list(J)%p 
!   task_list%taskEstimate(itask) = lhs_list(I)%p%task_part*rhs_list(I)%p%task_part
  ENDDO
ENDDO

deallocate(lhs_list)
deallocate(rhs_list)
call mem_dealloc(node_lhs)
call mem_dealloc(node_rhs)

END SUBROUTINE lsmpi_set_MPI_task_list

SUBROUTINE lsmpi_free_MPI_task_list(task_list)
implicit none
TYPE(LSMPI_TASK_LIST),INTENT(INOUT) :: task_list
!
CALL mem_dealloc(task_list%taskPair)
CALL mem_dealloc(task_list%taskNode)
DEALLOCATE(task_list%LHS)
DEALLOCATE(task_list%RHS)
!CALL mem_dealloc(task_list%taskWall)
!CALL mem_dealloc(task_list%taskCPU)
!CALL mem_dealloc(task_list%taskEstimate)

END SUBROUTINE lsmpi_free_MPI_task_list

#endif


    subroutine ls_print_tasks(tasks,iunit)
    implicit none
    TYPE(LS_TASKS),INTENT(IN) :: tasks
    INTEGER,INTENT(IN)        :: iunit
    !
    TYPE(LSTASK),pointer :: current
    INTEGER              :: itask,i

    WRITE(IUNIT,'(1X,A)')    '***        MPI task object       ***'
    WRITE(IUNIT,'(1X,A,I5)') 'Target ntasks:                      ', tasks%target_ntasks
    WRITE(IUNIT,'(1X,A,I5)') 'Acutal number of tasks:             ', tasks%ntasks
    WRITE(IUNIT,'(1X,A,I5)') 'Number of row atoms:                ', tasks%natoms_row
    WRITE(IUNIT,'(1X,A,I5)') 'Number of column atoms:             ', tasks%natoms_col
    WRITE(IUNIT,'(1X,A,I5)') 'Number of originiating row atoms:   ', tasks%natoms_row_full
    WRITE(IUNIT,'(1X,A,I5)') 'Number of originiating column atoms:', tasks%natoms_col_full
    WRITE(IUNIT,'(1X,A)')    ''

    current => tasks%first
    DO itask=1,tasks%ntasks
      IF (.NOT.ASSOCIATED(current)) CALL LSQUIT('Error in ls_print_tasks. Task not associated',IUNIT)
      WRITE(IUNIT,'(3X,A,I5)') 'Task number:                  ',itask
      WRITE(IUNIT,'(3X,A,F16.8)') 'Task partition:               ',current%task_part
      IF (current%full_mol_row) THEN
        WRITE(IUNIT,'(3X,A)')  'The row partition consist of the full molecule'
        WRITE(IUNIT,'(3X,A,I5)') 'Number of row atoms in task:      ',current%nrow
        WRITE(IUNIT,'(3X,A,I5)') 'Number of originating row atoms:  ',current%nrow_full
      ELSE
        WRITE(IUNIT,'(3X,A,I5)') 'Number of row atoms in task:      ',current%nrow
        WRITE(IUNIT,'(3X,A,I5)') 'Number of originating row atoms:  ',current%nrow_full
        IF (size(current%row_atoms).NE.current%nrow) &
     &      CALL LSQUIT('Error in ls_print_tasks. row_atoms inconsistency!',iunit)
        WRITE(IUNIT,'(5X,A7,10I6 /, (12X,10I6))')   'Atoms ',current%row_atoms
      ENDIF
      IF (current%full_mol_col) THEN
        WRITE(IUNIT,'(3X,A)')  'The column partition consist of the full molecule'
        WRITE(IUNIT,'(3X,A,I5)') 'Number of column atoms in task:    ',current%ncol
        WRITE(IUNIT,'(3X,A,I5)') 'Number of originating column atoms:',current%ncol_full
      ELSE
        WRITE(IUNIT,'(3X,A,I5)') 'Number of column atoms in task:    ',current%ncol
        WRITE(IUNIT,'(3X,A,I5)') 'Number of originating column atoms:',current%ncol_full
        IF (size(current%col_atoms).NE.current%ncol) &
     &      CALL LSQUIT('Error in ls_print_tasks. col_atoms inconsistency!',iunit)
        WRITE(IUNIT,'(5X,A7,10I6 /, (12X,10I6))')   'Atoms ',current%col_atoms
      ENDIF
      current => current%next
    ENDDO
    CALL LS_FLSHFO(IUNIT)
    IF (ASSOCIATED(current)) CALL LSQUIT('Error in ls_print_tasks. Tasks not terminated correctly!',-1)
    
    end subroutine ls_print_tasks
 
    subroutine ls_init_tasks(tasks,ntasks,nrow_atoms,ncol_atoms,nrow_atoms_full,ncol_atoms_full,&
    &                        no_row_part,no_col_part)
    implicit none
    TYPE(LS_TASKS),INTENT(INOUT) :: tasks
    INTEGER,INTENT(IN)           :: ntasks,nrow_atoms,ncol_atoms,nrow_atoms_full,ncol_atoms_full
    LOGICAL,INTENT(IN)           :: no_row_part,no_col_part
    tasks%target_ntasks   = ntasks
    tasks%ntasks          = 0
    tasks%natoms_row      = nrow_atoms
    tasks%natoms_col      = ncol_atoms
    tasks%natoms_row_full = nrow_atoms_full
    tasks%natoms_col_full = ncol_atoms_full
    tasks%no_row_part     = no_row_part
    tasks%no_col_part     = no_col_part
    NULLIFY(tasks%first)
    NULLIFY(tasks%last)
    end subroutine ls_init_tasks

    subroutine ls_free_tasks(tasks,lupri)
    implicit none
    TYPE(LS_TASKS),INTENT(INOUT) :: tasks
    INTEGER,INTENT(IN)           :: lupri
    !
    INTEGER :: itask
    TYPE(LSTASK),pointer :: current,next
    next => tasks%first
    DO itask=1,tasks%ntasks
      current => next
      next    => current%next
      call ls_free_task(current,lupri)
      DEALLOCATE(current)
      NULLIFY(current)
    ENDDO
    IF (ASSOCIATED(next)) THEN
      CALL LSQUIT('Error in ls_free_tasks. Something fishy: the last task is not last.',lupri)
    ENDIF
    end subroutine ls_free_tasks

    subroutine ls_free_task(task,lupri)
    implicit none
    TYPE(LSTASK),intent(INOUT) :: task
    INTEGER,INTENT(IN)         :: lupri
    !Consistency testing
    IF (.NOT.task%full_mol_row) THEN
      IF (size(task%row_atoms).NE.task%nrow) THEN
        write(lupri,'(1X,A,I6,A,I6)') 'Error in ls_free_task: row_atoms=',size(task%row_atoms),' .NE. nrow=',task%nrow
        CALL LSQUIT('Error in ls_free_task. Something fishy: nrow do not match the size of row_atoms.',lupri)
      ENDIF
      call mem_dealloc(task%row_atoms)
    ENDIF
    IF (.NOT.task%full_mol_col) THEN
      IF (size(task%col_atoms).NE.task%ncol) THEN
        write(lupri,'(1X,A,I6,A,I6)') 'Error in ls_free_task: col_atoms=',size(task%col_atoms),' .NE. ncol=',task%ncol
        CALL LSQUIT('Error in ls_free_task. Something fishy: ncol do not match the size of col_atoms.',lupri)
      ENDIF
      call mem_dealloc(task%col_atoms)
    ENDIF
    task%nrow = 0
    task%ncol = 0
    task%nrow_full = 0
    task%ncol_full = 0
    NULLIFY(task%next)
    NULLIFY(task%previous)
    end subroutine ls_free_task

    subroutine ls_add_task(tasks,nrow,ncol,row_to_task,col_to_task,row_id,col_id,task_part,lupri)
    implicit none
    TYPE(LS_TASKS),INTENT(INOUT) :: tasks
    Integer,intent(IN)           :: nrow,ncol,row_id,col_id,lupri
    Integer,INTENT(IN)           :: row_to_task(tasks%natoms_row),col_to_task(tasks%natoms_col)
    Real(realk),intent(IN)       :: task_part
    !
    TYPE(LSTASK),pointer :: old,new
    INTEGER              :: irow,icol,iatom

    NULLIFY(new)
    ALLOCATE(new)

    tasks%ntasks = tasks%ntasks + 1

!   If there exists a first task we must have a sequence of at least one 
!   task, and should therefore add the new task to the end of that sequence
    IF (ASSOCIATED(tasks%first)) THEN
      old => tasks%last
      IF (ASSOCIATED(old%next)) &
     &  CALL LSQUIT('Error in ls_add_task: old%next should not be associated!',lupri)
      old%next => new
      new%previous => old
      tasks%last => new
!   If not, we have to start with the task we are now adding
    ELSE
      tasks%first => new
      tasks%last => new
      nullify(new%previous)
    ENDIF
    nullify(new%next)
    
    new%nrow_full = tasks%natoms_row_full
    new%ncol_full = tasks%natoms_col_full
    IF ((nrow.EQ.-1).OR.(nrow.EQ.tasks%natoms_row)) THEN
      new%full_mol_row = .TRUE.
      new%nrow = tasks%natoms_row
      NULLIFY(new%row_atoms)
    ELSE
      new%full_mol_row = .FALSE.
      new%nrow = nrow
      call mem_alloc(new%row_atoms,nrow)
      irow = 0
      DO iatom=1,tasks%natoms_row
        IF (row_to_task(iatom).EQ.row_id) THEN
          irow = irow + 1
          new%row_atoms(irow) = iatom
        ENDIF
      ENDDO
      IF (irow.NE.nrow) THEN
        write(*,*) 'Error in ls_add_task. Strange irow.NE.nrow!'
        write(*,*) 'irow        = ',irow
        write(*,*) 'nrow        = ',nrow,tasks%natoms_row
        write(*,*) 'ncol        = ',ncol,tasks%natoms_col
        write(*,*) 'nrow_full   = ',nrow,tasks%natoms_row_full
        write(*,*) 'ncol_full   = ',ncol,tasks%natoms_col_full
        write(*,*) 'row_to_task = ',row_to_task
        write(*,*) 'col_to_task = ',col_to_task
        write(*,*) 'row_id      = ',row_id
        write(*,*) 'col_id      = ',col_id
        write(*,*) 'task_part   = ',task_part
        write(lupri,*) 'Error in ls_add_task. Strange irow.NE.nrow!'
        write(lupri,*) 'irow        = ',irow
        write(lupri,*) 'nrow        = ',nrow
        write(lupri,*) 'ncol        = ',ncol
        write(lupri,*) 'nrow_full   = ',tasks%natoms_row_full
        write(lupri,*) 'ncol_full   = ',tasks%natoms_col_full
        write(lupri,*) 'row_to_task = ',row_to_task
        write(lupri,*) 'col_to_task = ',col_to_task
        write(lupri,*) 'row_id      = ',row_id
        write(lupri,*) 'col_id      = ',col_id
        write(lupri,*) 'task_part   = ',task_part
        CALL LSQUIT('Error in ls_add_task. Strange irow.NE.nrow!',lupri)
      ENDIF
    ENDIF

    new%ncol      = ncol
    new%task_part = task_part
    IF ((ncol.EQ.-1).OR.(ncol.EQ.tasks%natoms_col)) THEN
      new%full_mol_col = .TRUE.
      new%ncol = tasks%natoms_col
      NULLIFY(new%col_atoms)
    ELSE
      new%full_mol_col = .FALSE.
      call mem_alloc(new%col_atoms,ncol)
      icol = 0
      DO iatom=1,tasks%natoms_col
        IF (col_to_task(iatom).EQ.col_id) THEN
          icol = icol + 1
          new%col_atoms(icol) = iatom
        ENDIF
      ENDDO
      IF (icol.NE.ncol) CALL LSQUIT('Error in ls_add_task. Strange icol.NE.ncol!',-1)
    ENDIF

    end subroutine ls_add_task

    ! Time estimates in timeMatrix must be positive
    SUBROUTINE ls_time_tasks(tasks,timeMatrix,natoms_row,natoms_col,ntasks,lupri,iprint)
    implicit none
    Integer,intent(IN)     :: natoms_col,natoms_row,ntasks,lupri,iprint
    real(realk),intent(IN) :: timeMatrix(natoms_row,natoms_col)
    type(ls_tasks)         :: tasks
    !
    real(realk),pointer :: column_sum(:),column_factor(:),fraction_sum(:),fraction_deviation(:)
    real(realk),pointer :: row_sum(:),row_factor(:),task_sum(:),task_deviation(:)
    real(realk)         :: time_sum,time_average,part_sum,deviation,new_deviation,max_deviation,min_deviation
    real(realk)         :: new_sum,max_sum,max_task_sum,fraction_fac
    integer             :: i,j,k,nearest_integer,min_dev,max_dev,ifraction,nfractions,itask,jtask
    integer             :: start_task,end_task
    logical,pointer     :: column_assigned(:),row_assigned(:)
    integer,pointer     :: column_to_fraction(:),fraction_nint(:),fraction_row_to_task(:)
    integer,pointer     :: fraction_ncol(:),task_nrow(:)
    logical             :: reduction
    integer             :: ncombine,icombine
    
    call mem_alloc(column_sum,natoms_col)
    call mem_alloc(column_factor,natoms_col)
    call mem_alloc(column_assigned,natoms_col)
    call mem_alloc(column_to_fraction,natoms_col)
    call mem_alloc(fraction_sum,natoms_col)
    call mem_alloc(fraction_deviation,natoms_col)
    call mem_alloc(fraction_nint,natoms_col)
    
    ! We would like to split the full calculation up into ntasks different tasks,
    ! each corresponding to a subblock of the timeMatrix (non-adjacent column or 
    ! row elements allowed)
    !
    ! Our strategy is to split the timeMatrix using a two step procedure:
    !
    !   1. Combine columns info fractions first - in such a way that each
    !      fraction should take close to a multiplum of the average 
    !      task fraction
    !
    !   2. Each fraction is then split into the final subblocks, each such
    !      subblock should ideally now take (one times) the average task fraction
    
    

    
    !*************************      STEP 1 
    !   1. Combine columns info fractions first - in such a way that each
    !      fraction should take close to a multiplum of the average 
    !      task fraction

    IF (tasks%no_col_part) THEN
      nfractions = 1
      DO i=1,natoms_col
        column_to_fraction(i) = i
      ENDDO
      fraction_sum(nfractions)       = 1E0_realk*ntasks
      fraction_deviation(nfractions) = 0E0_realk
      fraction_nint(nfractions)      = ntasks
      call mem_alloc(fraction_ncol,1)
      fraction_ncol(1)               = natoms_col
    ELSE
      ! Find first the total and column-wise time estimates
      time_sum = 0E0_realk
      DO j=1,natoms_col
        column_sum(j) = 0E0_realk
        DO i=1,natoms_row
          column_sum(j) = column_sum(j) + timeMatrix(i,j)
        ENDDO
        time_sum = time_sum + column_sum(j)
IF (iprint.GT.3) THEN
  write(lupri,*) 'column_sum',j,column_sum(j)
  call ls_flshfo(lupri)
ENDIF
      ENDDO
      
      time_average = time_sum/ntasks
IF (iprint.GT.3) THEN
  write(lupri,*) ''
  write(lupri,*) 'time_average',time_average,time_sum
  write(lupri,*) ''
  call ls_flshfo(lupri)
ENDIF
      
      ! Calculate the task-fraction of each column: with the fractional 
      ! sums adding up to the total number of tasks
      DO j=1,natoms_col
        column_factor(j) = column_sum(j)/time_average
IF (iprint.GT.3) THEN
  write(lupri,*) 'column_factor',j,column_factor(j)
  call ls_flshfo(lupri)
ENDIF
      ENDDO
      
      
      ! Combine columns together into fractions, in order to construct near-integer 
      ! time-fractions (with at least a full task in each fraction)
      column_assigned = .FALSE.
      ifraction = 0
      fraction_sum = 0E0_realk
      DO j=1,natoms_col
        IF (column_assigned(j)) CYCLE
        ifraction = ifraction + 1
        column_assigned(j) = .TRUE.
        column_to_fraction(j) = ifraction
        part_sum = column_factor(j)
      ! We would like to have at least one fraction, so close to zero is not favorable
        nearest_integer = max(1,nint(part_sum))
      ! If row is empty we would like to hit right on target with the column partitioning
        IF (tasks%no_row_part) nearest_integer = 1
        deviation = part_sum - nearest_integer
IF (iprint.GT.3) THEN
  write(lupri,*) 'deviation',j,deviation
  call ls_flshfo(lupri)
ENDIF
        DO !while not finished
          min_dev = -1
          DO i=j+1,natoms_col
            IF (column_assigned(i)) CYCLE
            new_sum = part_sum + column_factor(i)
            nearest_integer = max(1,nint(new_sum))
            IF (tasks%no_row_part) nearest_integer = 1
            new_deviation   = new_sum - nearest_integer
            IF (abs(deviation).GE.abs(new_deviation)) THEN
IF (iprint.GT.3) THEN
  write(lupri,*) 'new_deviation',i,new_deviation
  call ls_flshfo(lupri)
ENDIF
              min_dev = i
              deviation = new_deviation
            ENDIF
          ENDDO
          IF (min_dev.EQ.-1) THEN
            EXIT
          ELSE
            column_assigned(min_dev) = .TRUE.
            column_to_fraction(min_dev) = ifraction
            part_sum = part_sum + column_factor(min_dev)
          ENDIF
        ENDDO
        fraction_sum(ifraction)       = part_sum
        fraction_deviation(ifraction) = deviation
        fraction_nint(ifraction)      = nint(part_sum)
      ENDDO
      
      nfractions = ifraction
IF (iprint.GT.3) THEN
  write(lupri,*) 'column_to_fraction',nfractions,column_to_fraction
  call ls_flshfo(lupri)
ENDIF
      
      ! Add unassigned columns to exisiting tasks
      DO j=1,natoms_col
        IF (column_assigned(j)) CYCLE
        column_assigned(j) = .TRUE.
        deviation     = 1E99_realk
        max_deviation = 1E99_realk
        min_dev = -1
        max_dev = -1
        DO ifraction=1,nfractions
          new_sum = fraction_sum(ifraction) + column_factor(j)
          nearest_integer = nint(new_sum)
          IF (tasks%no_row_part) nearest_integer = 1
          new_deviation = new_sum - nearest_integer
          IF (new_deviation.LT.deviation) THEN
            min_dev = ifraction
            deviation = new_deviation
          ENDIF
          new_deviation = new_sum - fraction_nint(ifraction) - 1
          IF (new_deviation.LT.max_deviation) THEN
            max_dev = ifraction
            max_deviation = new_deviation
          ENDIF
        ENDDO
      
IF (iprint.GT.3) THEN
  write(lupri,*) 'unassigned',min_dev,deviation,max_dev,max_deviation
  call ls_flshfo(lupri)
ENDIF
      ! Add column to fraction so than deviation becomes smaller. If all deviation becomes larger
      ! add in a way so that the deviation to the next integer becomes as small as possible.
      ! For example imagine we have four fractions with fractional sums 1.33, 1.32, 0.18 and 0.17. 
      ! We cannot reduce the overall deviation by addint 0.18 to one of the first fractions. Therfore 
      ! we add it to the first to reduce the overall deviation from the nearest integer pluss one.
      ! Ie. to get 1.51, 1.32 and 0.17. We can then reduce the overall deviation by combining the 
      ! first and third fractions, to get two fractions 1.68 and 1.32
        IF (min_dev.EQ.-1) THEN
          IF (max_dev.EQ.-1) THEN
            CALL LSQUIT('Error in ls_time_tasks. Cannot add column to fraction!',lupri)
          ENDIF
          deviation = max_deviation
          min_dev   = max_dev
        ENDIF
        fraction_sum(min_dev) = fraction_sum(min_dev) + column_factor(j)
        fraction_nint(min_dev) = nint(fraction_sum(min_dev))
        fraction_deviation(min_dev) = deviation
        column_to_fraction(j) = min_dev
      ENDDO
      
      ! Fractions with close to zero time-fractions are joined with other fractions
      !
      ! Simen: Most likely better to sort fractions according to decreasing 
      !        time-fractions before this step
      !
      DO ifraction=1,nfractions
        IF (nint(fraction_sum(ifraction)).EQ.0) THEN
          deviation = 1E99_realk
          max_deviation = 1E99_realk
          min_dev = -1
          max_dev = -1
          DO i=1,nfractions
            IF (fraction_nint(i).EQ.0) CYCLE
            new_sum = fraction_sum(ifraction) + fraction_sum(i)
            nearest_integer = nint(new_sum)
            IF (tasks%no_row_part) nearest_integer = 1
            new_deviation   = new_sum - nearest_integer
            IF (abs(new_deviation).LT.abs(deviation)) THEN
              min_dev = i
              part_sum = new_sum
              deviation = new_deviation
            ENDIF
            new_deviation   = new_sum - nearest_integer - 1
            IF (abs(new_deviation).LT.abs(max_deviation)) THEN
              max_dev = i
              max_sum = new_sum
              max_deviation = new_deviation
            ENDIF
          ENDDO
          IF (min_dev.EQ.-1) THEN
            IF (max_dev.EQ.1) CALL LSQUIT('Error in ls_time_tasks: could not remove zero time-fraction',lupri)
            min_dev   = max_dev
            part_sum  = max_sum
            deviation = max_deviation
          ENDIF
IF (iprint.GT.3) THEN
  write(lupri,*) 'zero',min_dev,part_sum,deviation
  call ls_flshfo(lupri)
ENDIF
      !   Move first the min_dev fraction into the ifraction (zero-)fraction
          fraction_sum(ifraction)       = part_sum
          fraction_deviation(ifraction) = deviation
          fraction_nint(ifraction)      = nint(part_sum)
          DO j=1,natoms_col
            IF (column_to_fraction(j).EQ.min_dev) column_to_fraction(j)=ifraction
          ENDDO
      !   Then move the nfractions fraction into the now empty min_dev fraction
          IF (min_dev.NE.nfractions) THEN
            fraction_sum(min_dev)       = fraction_sum(nfractions)
            fraction_deviation(min_dev) = fraction_deviation(nfractions)
            fraction_nint(min_dev)      = fraction_nint(nfractions)
            DO j=1,natoms_col
              IF (column_to_fraction(j).EQ.nfractions) column_to_fraction(j)=min_dev
            ENDDO
          ENDIF
          nfractions = nfractions - 1
        ENDIF
      ENDDO
      
      ! Get the number of columns included in each fraction
      call mem_alloc(fraction_ncol,nfractions)
      fraction_ncol = 0
      DO i=1,natoms_col
        ifraction = column_to_fraction(i)
        fraction_ncol(ifraction) = fraction_ncol(ifraction) + 1
      ENDDO
IF (iprint.GT.3) THEN
  write(lupri,*) 'fraction_ncol',fraction_ncol
  call ls_flshfo(lupri)
ENDIF
    ENDIF
    
    !Check for consistency
    itask = 0
    DO ifraction=1,nfractions
IF (iprint.GT.3) THEN
  write(lupri,*) 'fraction',ifraction,fraction_sum(ifraction),fraction_deviation(ifraction),&
     &                            fraction_deviation(ifraction),fraction_nint(ifraction)
  call ls_flshfo(lupri)
ENDIF
      itask = itask + fraction_nint(ifraction)
    ENDDO
    do while (itask.LT.ntasks)
      ncombine = ntasks - itask
      DO icombine=1,ncombine
!       Combine the two tasks with the highest positive deviation
        min_dev = -1
        max_dev = -1
        max_deviation = 0_realk
        min_deviation = 0_realk
        DO ifraction=1,nfractions
          IF (fraction_deviation(ifraction).GE.max_deviation) THEN
            min_dev = max_dev
            min_deviation = max_deviation
            max_dev = ifraction
            max_deviation = max(max_deviation,fraction_deviation(ifraction))
          ELSEIF (fraction_deviation(ifraction).GE.min_deviation) THEN
            min_dev = ifraction
            min_deviation = max(min_deviation,fraction_deviation(ifraction))
          ENDIF
        ENDDO
        IF ((max_dev.EQ.-1).OR.(min_dev.EQ.-1)) THEN
          write(LUPRI,'(1X,A,I5,A,I5)') 'Error in ls_time_tasks. Unable to combine tasks into correct number'
          write(lupri,*) 'fraction_sum      ',fraction_sum(1:nfractions)
          write(lupri,*) 'fraction_deviation',fraction_deviation(1:nfractions)
          write(lupri,*) 'timeMatrix',timeMatrix
          CALL LSQUIT('Unable to combine tasks. Needs implementation',lupri)
        ENDIF
        reduction = lstask_combine_fragments(fraction_sum,fraction_deviation,fraction_nint,&
     &                                       fraction_ncol,column_to_fraction,min_dev,max_dev,nfractions,&
     &                                       natoms_col,.FALSE.,lupri)
        IF (reduction) itask = itask + 1
      ENDDO
    enddo
    do while (itask.GT.ntasks)
      ncombine = itask - ntasks
      DO icombine=1,ncombine
!       Combine the two tasks with the highest negative deviation
        min_dev = -1
        max_dev = -1
        max_deviation = 0_realk
        min_deviation = 0_realk
        DO ifraction=1,nfractions
          IF (fraction_deviation(ifraction).LE.max_deviation) THEN
            min_dev = max_dev
            min_deviation = max_deviation
            max_dev = ifraction
            max_deviation = min(max_deviation,fraction_deviation(ifraction))
          ELSEIF (fraction_deviation(ifraction).LE.min_deviation) THEN
            min_dev = ifraction
            min_deviation = min(min_deviation,fraction_deviation(ifraction))
          ENDIF
        ENDDO
        IF ((max_dev.EQ.-1).OR.(min_dev.EQ.-1)) THEN
          write(LUPRI,'(1X,A,I5,A,I5)') 'Error in ls_time_tasks. Unable to combine tasks into correct number'
          write(lupri,*) 'fraction_sum      ',fraction_sum(1:nfractions)
          write(lupri,*) 'fraction_deviation',fraction_deviation(1:nfractions)
          write(lupri,*) 'timeMatrix',timeMatrix
          CALL LSQUIT('Unable to combine tasks. Needs implementation',lupri)
        ENDIF
        reduction = lstask_combine_fragments(fraction_sum,fraction_deviation,fraction_nint,&
     &                                  fraction_ncol,column_to_fraction,min_dev,max_dev,nfractions,&
     &                                  natoms_col,.TRUE.,lupri)
IF (iprint.GT.3) THEN
  write(lupri,*) 'reduction+',itask,ntasks,min_dev,max_dev,reduction,min_deviation,max_deviation,nfractions
  call ls_flshfo(lupri)
ENDIF
        IF (reduction) itask = itask - 1
      ENDDO
    enddo
    IF (itask.NE.ntasks) THEN
      write(LUPRI,'(1X,A,I5,A,I5)') 'Error in ls_time_tasks. Wrong number of fractional tasks. itask =',itask,' ntasks=',ntasks
      write(lupri,*) 'fraction_sum      ',fraction_sum(1:nfractions)
      write(lupri,*) 'fraction_deviation',fraction_deviation(1:nfractions)
      CALL LSQUIT('Wrong number of fractional tasks. Needs implementation',lupri)
    ! We may here possible split up one fraction, or combine two or more fractions, in order 
    ! to make things add up to the correct number of task partitions
    ENDIF
    
IF (iprint.GT.3) THEN
  write(lupri,*) 'column_to_fraction',column_to_fraction
  call ls_flshfo(lupri)
ENDIF
    DO ifraction=1,nfractions
IF (iprint.GT.3) THEN
  write(lupri,*) 'ifraction',ifraction,fraction_sum(ifraction),fraction_deviation(ifraction),&
     &                             fraction_nint(ifraction)
  call ls_flshfo(lupri)
ENDIF
    ENDDO
    
    !*************************      STEP 2
    !
    ! Split fractions into the final subblocks (each such subblock should ideally now 
    ! take (one times) the average task fraction)
    !
    call mem_alloc(fraction_row_to_task,natoms_row)
    call mem_alloc(row_sum,natoms_row)
    call mem_alloc(row_factor,natoms_row)
    call mem_alloc(row_assigned,natoms_row)
    call mem_alloc(task_deviation,ntasks)
    call mem_alloc(task_sum,ntasks)
    call mem_alloc(task_nrow,ntasks)

    max_task_sum = 0E0_realk
    
    itask = 0
    DO ifraction=1,nfractions
    
    ! Copy fraction to task if the fragment only contains one task
      IF (fraction_nint(ifraction).EQ.1) THEN
        itask = itask + 1
        IF (tasks%no_col_part) THEN
          max_task_sum = max(max_task_sum,1E0_realk)
IF (iprint.GT.3) THEN
  write(lupri,*) 'addtask1'
  call ls_flshfo(lupri)
ENDIF
          call ls_add_task(tasks,-1,-1,fraction_row_to_task,&
           &               column_to_fraction,itask,ifraction,1E0_realk,lupri)
        ELSE
          max_task_sum = max(max_task_sum,fraction_sum(ifraction))
IF (iprint.GT.3) THEN
  write(lupri,*) 'addtask2'
  call ls_flshfo(lupri)
ENDIF
          call ls_add_task(tasks,-1,fraction_ncol(ifraction),fraction_row_to_task,&
           &               column_to_fraction,itask,ifraction,fraction_sum(ifraction),lupri)
        ENDIF
        CYCLE
      ENDIF

    ! Else make a partioning of given fraction into fraction_nint tasks

    ! Start by getting the row_sums and row_factors
      row_sum = 0E0_realk
      DO j=1,natoms_col
        IF (column_to_fraction(j).NE.ifraction) CYCLE
        DO i=1,natoms_row
          row_sum(i) = row_sum(i) + timeMatrix(i,j)
        ENDDO
      ENDDO
IF (iprint.GT.3) THEN
  write(lupri,*) 'row_sum',ifraction,row_sum,natoms_row
  call ls_flshfo(lupri)
ENDIF
      part_sum = 0E0_realk
      DO i=1,natoms_row
        part_sum = part_sum + row_sum(i)
      ENDDO
      fraction_fac = fraction_sum(ifraction)/fraction_nint(ifraction)
      time_average = part_sum/ntasks
    ! row_factor = row_sum/time_average
      row_factor = fraction_nint(ifraction)*row_sum/part_sum
IF (iprint.GT.3) THEN
  write(lupri,*) 'row_sum',row_sum
  write(lupri,*) 'row_factor',row_factor
  write(lupri,*) 'fraction_nint',ifraction,fraction_nint(ifraction),natoms_row
  call ls_flshfo(lupri)
ENDIF
      
    ! Then assign the rows to (the fraction_nint) tasks based on the row_factors 
      row_assigned = .FALSE.
      start_task = itask+1
      end_task   = itask+fraction_nint(ifraction)
    ! We now create fraction_nint(ifraction) new tasks
      DO k=1,fraction_nint(ifraction)
    !   For each task combine rows so that the task fraction becomes as close to 
    !   one as possible - and then distribute the remaning rows the best possible way
        DO j=1,natoms_row
          IF (row_assigned(j)) CYCLE
    !     Special case when we have already created fraction_nint new tasks:
    !              - add row to an existing task
          IF (itask.EQ.end_task) THEN
            min_dev   = start_task
            deviation = task_sum(start_task) + row_factor(j) - 1E0_realk
            DO jtask=start_task+1,end_task
              new_sum = task_sum(jtask) + row_factor(j)
              new_deviation = new_sum - 1E0_realk
              IF (abs(deviation).GE.abs(new_deviation)) THEN
                min_dev = jtask
                deviation = new_deviation
              ENDIF
            ENDDO
    !       Add fraction to min_dev
            row_assigned(j) = .TRUE.        
            fraction_row_to_task(j) = min_dev
            task_sum(min_dev) = task_sum(min_dev) + row_factor(j)
            task_nrow(min_dev) = task_nrow(min_dev) + 1
            CYCLE
          ENDIF
          itask = itask + 1
          task_nrow(itask) = 1
          row_assigned(j) = .TRUE.
          fraction_row_to_task(j) = itask
          part_sum = row_factor(j)
    !     Each task should take the average task fraction (i.e. one)
          deviation = part_sum - 1E0_realk
          DO !while not finished
            min_dev = -1
            DO i=j+1,natoms_row
              IF (row_assigned(i)) CYCLE
              new_sum = part_sum + row_factor(i)
              new_deviation   = new_sum - 1E0_realk
              IF (abs(deviation).GE.abs(new_deviation)) THEN
                min_dev = i
                deviation = new_deviation
              ENDIF
            ENDDO
            IF (min_dev.EQ.-1) THEN
              EXIT
            ELSE
              row_assigned(min_dev) = .TRUE.
              fraction_row_to_task(min_dev) = itask
              part_sum = part_sum + row_factor(min_dev)
              task_nrow(itask) = task_nrow(itask) + 1
            ENDIF
          ENDDO
          task_sum(itask)       = part_sum
          task_deviation(itask) = deviation
        ENDDO
      ENDDO
      
      end_task=itask
      IF (itask.GT.ntasks) CALL LSQUIT('Error in ls_time_tasks. Two many tasks created',lupri)
    ! Add any unassigned row to the fragment tasks
IF (iprint.GT.10) THEN
  write(lupri,*) 'debug:assigned rows',row_assigned
ENDIF
      DO j=1,natoms_row
        IF (row_assigned(j)) CYCLE
        min_dev = -1
        deviation = 1E99_REALK
        DO k=start_task,end_task
          new_sum = task_sum(k) + row_factor(j)
          new_deviation = new_sum - 1E0_realk
          IF (abs(deviation).GE.abs(new_deviation)) THEN
            min_dev = k
            deviation = new_deviation
          ENDIF
        ENDDO
        IF (min_dev.EQ.-1) THEN
          write(lupri,*) 'Error in ls_time_tasks. Deviation too large'
          write(lupri,*) 'timeMatrix',natoms_row,natoms_col,ntasks
          call ls_output(timeMatrix,1,natoms_row,1,natoms_col,natoms_row,natoms_col,1,lupri)
          CALL LSQUIT('Error in ls_time_tasks. Deviation too large',lupri)
        ELSE
IF (iprint.GT.10) THEN
  write(lupri,*) 'debug:unassigned row',j,min_dev,task_sum(min_dev) + row_factor(j),task_nrow(itask) + 1,start_task,end_task
ENDIF
          row_assigned(j) = .TRUE.
          fraction_row_to_task(j) = min_dev
          task_sum(min_dev) = task_sum(min_dev) + row_factor(j)
          task_nrow(itask) = task_nrow(itask) + 1
        ENDIF 
      ENDDO
      DO k=start_task,end_task 
          task_sum(k) = task_sum(k)*fraction_fac
          max_task_sum = max(max_task_sum,task_sum(k))
IF (iprint.GT.10) THEN
  write(lupri,*) 'addtask3',ifraction,k,ntasks,task_sum(k),task_deviation(k),task_nrow(k),&
    &task_nrow(k),fraction_ncol(ifraction),k,ifraction
  write(lupri,*) 'addtask4',fraction_row_to_task
  write(lupri,*) 'addtask5',column_to_fraction
  write(lupri,*) 'addtask6',row_assigned
  call ls_flshfo(lupri)
ENDIF
          !Copy fraction to task
          call ls_add_task(tasks,task_nrow(k),fraction_ncol(ifraction),fraction_row_to_task,&
           &               column_to_fraction,k,ifraction,task_sum(k),lupri)
IF (iprint.GT.10) THEN
  write(lupri,*) 'addtask6'
  call ls_flshfo(lupri)
ENDIF
      ENDDO
    ENDDO !ifraction


IF (iprint.GT.1) THEN
    write(lupri,'(1X,A,I6)') 'Starting mpi with number of tasks =',tasks%ntasks
    write(lupri,'(1X,A,F5.2,A,F5.1,A)') 'Max task partition is',max_task_sum,&
         &      ' (optimal 1), giving an estimated efficiency of',1E2_realk/max_task_sum,'%'
  call ls_flshfo(lupri)
ENDIF
    
    call mem_dealloc(row_sum)
    call mem_dealloc(row_factor)
    call mem_dealloc(row_assigned)
    call mem_dealloc(fraction_row_to_task)
    call mem_dealloc(task_deviation)
    call mem_dealloc(task_sum)
    call mem_dealloc(task_nrow)
    call mem_dealloc(column_sum)
    call mem_dealloc(column_factor)
    call mem_dealloc(column_assigned)
    call mem_dealloc(column_to_fraction)
    call mem_dealloc(fraction_sum)
    call mem_dealloc(fraction_deviation)
    call mem_dealloc(fraction_nint)
    call mem_dealloc(fraction_ncol)

    END SUBROUTINE ls_time_tasks

    FUNCTION lstask_combine_fragments(fraction_sum,fraction_deviation,fraction_nint,&
     &                                fraction_ncol,column_to_fraction,ifrag,jfrag,nfrag,natoms,add,lupri)
    implicit none
    integer,intent(IN)        :: ifrag,jfrag,natoms,lupri
    real(realk),intent(inout) :: fraction_sum(natoms)
    real(realk),intent(inout) :: fraction_deviation(natoms)
    integer,intent(inout)     :: fraction_nint(natoms)
    integer,intent(inout)     :: fraction_ncol(natoms)
    integer,intent(inout)     :: column_to_fraction(natoms)
    integer,intent(INOUT)     :: nfrag
    logical,intent(IN)        :: add
    logical                   :: lstask_combine_fragments
!
    real(realk) :: new_sum
    real(realk) :: new_deviation
    integer     :: new_nint
    logical     :: reduction
    integer     :: iatom,new_ncol

    new_sum       = fraction_sum(ifrag) + fraction_sum(jfrag)
    new_deviation = fraction_deviation(ifrag) + fraction_deviation(jfrag)
    new_ncol      = fraction_ncol(ifrag) + fraction_ncol(jfrag)
    new_nint      = nint(new_sum)

    IF (add) THEN
      reduction = ((new_nint - fraction_nint(ifrag) - fraction_nint(jfrag)).EQ.-1)
      IF (reduction) new_deviation = new_deviation + 1.0E0_realk
    ELSE
      reduction = ((new_nint - fraction_nint(ifrag) - fraction_nint(jfrag)).EQ.1)
      IF (reduction) new_deviation = new_deviation - 1.0E0_realk
    ENDIF

    !Add jfrag to ifrag
    fraction_sum(ifrag)       = new_sum
    fraction_deviation(ifrag) = new_deviation
    fraction_nint(ifrag)      = new_nint
    fraction_ncol(ifrag)      = new_ncol
    DO iatom=1,natoms
      IF (column_to_fraction(iatom).EQ.jfrag) column_to_fraction(iatom) = ifrag
    ENDDO

    !Move nfrag to the (now empty) jfrag before reducing nfrag by one
    IF (jfrag.NE.nfrag) THEN
      fraction_sum(jfrag)       = fraction_sum(nfrag)
      fraction_deviation(jfrag) = fraction_deviation(nfrag)
      fraction_nint(jfrag)      = fraction_nint(nfrag)
      fraction_ncol(jfrag)      = fraction_ncol(nfrag)
      DO iatom=1,natoms
        IF (column_to_fraction(iatom).EQ.nfrag) column_to_fraction(iatom) = jfrag
      ENDDO
    ENDIF

    nfrag = nfrag - 1
    lstask_combine_fragments = reduction
    END FUNCTION lstask_combine_fragments

#ifdef VAR_MPI
    SUBROUTINE jengine_time_estimate(tmLHS,tmRHS,nA,nB,nC,nD,sameAB,sameCD,noA,noB,noC,noD,&
     &                               AO1,AO2,AO3,AO4,Spec,setting,lupri)
    implicit none
    Real(realk),pointer           :: tmLHS(:,:)
    Real(realk),pointer           :: tmRHS(:,:)
    integer                       :: AO1,AO2,AO3,AO4,Spec
    integer,intent(IN)            :: lupri
    type(lssetting),intent(inout) :: setting
    logical                       :: sameAB,sameCD,noA,noB,noC,noD
    integer                       :: nA,nB,nC,nD
    !
    Integer                      :: power
    integer                      :: AO(4)
    type(molecule_pt)            :: mol(4)
!    type(basisset_pt)            :: bas(4)
    integer                      :: nAtoms(4)
    logical                      :: noPart(4)
    integer                      :: iA,iB,iC,iD,iAO,indAB,indCD,startB,startD
    real(realk)                  :: R2,tmp,sum
    integer                      :: nbatchA,nbatchB,nbatchC,nbatchD
    integer                      :: ibatchA,ibatchB,ibatchC,ibatchD
    real(realk),parameter        :: one=1E0_realk, zero=0E0_realk
    TYPE(AOBATCHINFO),pointer    :: aoA,aoB,aoC,aoD
    integer(kind=short),pointer  :: LHSGAB(:,:),RHSGAB(:,:),LHSDENS(:,:),RHSDENS(:,:)
    real(realk)                  :: ABatomPairTE,CDatomPairTE
    real(realk)                  :: batchTE
    Logical                      :: screen,LHSDMATset,RHSDMATset,useScreen
    integer(short)               :: outerCsProduct,csProduct,CS_THRLOG,LHS_CSTHR,RHS_CSTHR
    integer(kind=short)          :: maxgab
    integer :: globalstartbatchA,globalstartbatchB,globalstartbatchC,globalstartbatchD
    integer :: globalbatchB,globalbatchD
    !Set up power - used for the time estimate of batch-wise integrals
    power = 7
    IF (AO1.EQ.AOEmpty) power = power - 1
    IF (AO2.EQ.AOEmpty) power = power - 1
    IF (AO3.EQ.AOEmpty) power = power - 1
    IF (AO4.EQ.AOEmpty) power = power - 1
    IF (Spec.EQ.GradientSpec)  power = power + 1

    useScreen = ASSOCIATED(setting%redCS%LHSGAB).AND.ASSOCIATED(setting%redCS%RHSGAB)
    IF(useScreen)THEN
       LHSGAB  => setting%redCS%LHSGAB
       RHSGAB  => setting%redCS%RHSGAB
       CS_THRLOG  = setting%redCS%CS_THRLOG
       LHS_CSTHR  = CS_THRLOG - setting%redCS%maxgabRHS
       RHS_CSTHR  = CS_THRLOG - setting%redCS%maxgabLHS

       LHSDMATset = setting%redCS%LHSDMATset
       RHSDMATset = setting%redCS%RHSDMATset
       IF (RHSDMATset)THEN
          LHS_CSTHR=LHS_CSTHR-setting%redCS%maxDmatRHS
          RHSDENS => setting%redCS%RHSDMAT
       ENDIF
       IF (LHSDMATset)THEN
          RHS_CSTHR=RHS_CSTHR-setting%redCS%maxDmatLHS
          LHSDENS => setting%redCS%LHSDMAT
       ENDIF
    ENDIF
    !Set up molecules and basis
    mol(1)%p => setting%molecule(1)%p
    mol(2)%p => setting%molecule(2)%p
    mol(3)%p => setting%molecule(3)%p
    mol(4)%p => setting%molecule(4)%p
    AO(1) = AO1
    AO(2) = AO2
    AO(3) = AO3
    AO(4) = AO4
    DO iAO=1,4
      IF (ASSOCIATED(mol(iAO)%p)) THEN
        nAtoms(iAO) = mol(iAO)%p%nAtomsNPC
      ELSE
        nAtoms(iAO) = 1
      ENDIF
      noPart(iAO) = .FALSE.
      IF (AO(iAO).EQ.AOEmpty) THEN
        nAtoms(iAO) = 1
        noPart(iAO) = .TRUE.
      ELSE IF (AO(iAO).EQ.AOdfAux) THEN
!        bas(iAO)%p => setting%basis(iAO)%p%BINFO(AuxBasParam)
      ELSE IF (AO(iAO).EQ.AORegular) THEN
!        bas(iAO)%p => setting%basis(iAO)%p%BINFO(RegBasParam)
      ELSE IF (AO(iAO).EQ.AOVAL) THEN
!        bas(iAO)%p => setting%basis(iAO)%p%BINFO(ValBasParam)
      ELSE IF (AO(iAO).EQ.AOdfCABS) THEN
!        bas(iAO)%p => setting%basis(iAO)%p%BINFO(CABBasParam)
!         will not work
         nAtoms(iAO) = 2*mol(iAO)%p%nAtomsNPC
      ELSE IF (AO(iAO).EQ.AOdfCABO) THEN
!        bas(iAO)%p => setting%basis(iAO)%p%BINFO(CABBasParam)
      ELSE IF (AO(iAO).EQ.AOdfJK) THEN
!        bas(iAO)%p => setting%basis(iAO)%p%BINFO(JKBasParam)
      ELSE IF (AO(iAO).EQ.AOadmm) THEN
!        bas(iAO)%p => setting%basis(iAO)%p%BINFO(ADMBasParam)
      ELSE IF (AO(iAO).EQ.AOpCharge) THEN
        nAtoms(iAO) = mol(iAO)%p%nAtoms
        noPart(iAO) = .TRUE.
      ELSE IF (AO(iAO).EQ.AOelField) THEN
        nAtoms(iAO) = mol(iAO)%p%nAtoms
        noPart(iAO) = .TRUE.
      ELSE IF (AO(iAO).EQ.AONuclear) THEN
        nAtoms(iAO) = 1
        noPart(iAO) = .TRUE.
      ELSE
        WRITE(LUPRI,'(1X,2A)') 'Error in jengine_time_estimate:',iAO,AO(iAO)
        CALL LSQUIT('Error in jengine_time_estimate',lupri)
      ENDIF
!     Special case for batchindex only
      if (setting%Batchindex(iao).NE.0) nAtoms(iAO) = 1
    ENDDO
    nA = nAtoms(1)
    nB = nAtoms(2)
    nC = nAtoms(3)
    nD = nAtoms(4)
    noA = noPart(1)
    noB = noPart(2)
    noC = noPart(3)
    noD = noPart(4)
    
    call mem_alloc(tmLHS,nA,nB)
    call mem_alloc(tmRHS,nC,nD)
    
!    Zero out timeMatrices
     DO iB=1,nB
       DO iA=1,nA
         tmLHS(iA,iB) = ZERO
       ENDDO
     ENDDO
     DO iD=1,nD
       DO iC=1,nC
         tmRHS(iC,iD) = ZERO
       ENDDO
     ENDDO
#if 1
!    Then add the batch-wise time-estimates to the atomic pair estimates
     IF(useScreen)THEN
      IF(size(setting%redCS%AO(1)%batch).NE.nA)THEN
       !special for EP type integrals where 
       ! sum_CD (nuclear empty|CD)D_CD
       ! is done
       DO iA=1,nA
         DO iB=1,nB
           tmLHS(iA,iB) = 1.0E0_realk
         ENDDO
       ENDDO
      ELSE
       DO iA=1,nA!setting%redCS%AO(1)%nAtoms
        aoA => setting%redCS%AO(1)%batch(iA)
        globalstartbatchA = aoA%GlobalStartBatchindex -1 
        startB = 1
        IF (sameAB) startB=iA
        DO iB=startB,nB!setting%redCS%AO(2)%nAtoms
         aoB => setting%redCS%AO(2)%batch(iB)
         globalstartbatchB = aoB%GlobalStartBatchindex - 1
         ABatomPairTE = 0E0_realk
         DO ibatchB=1,aoB%nBatches
          globalbatchB = globalstartbatchB+ibatchB
          DO ibatchA=1,aoA%nBatches
           maxgab = LHSGAB(globalstartbatchA+ibatchA,globalbatchB)
           IF (maxgab.NE.shortzero) THEN
            IF (LHSDMATset) THEN
             csProduct = maxgab + LHSDENS(globalstartbatchA+ibatchA,globalbatchB )
            ELSE
             csProduct = maxgab
            ENDIF
            screen = (csProduct.LT.LHS_CSTHR)
            IF (.not.screen) THEN
             batchTE = aoA%nPrim(ibatchA)*aoB%nPrim(iBatchB)&
                  & *(aoA%maxAng(iBatchA)+aoB%maxAng(iBatchB)+10E0_realk)**POWER
             ABatomPairTE = ABatomPairTE + batchTE
            ENDIF !screen
           ENDIF !maxgab.NE.shortzero
          ENDDO !ibatchB
         ENDDO !ibatchA
         tmLHS(iA,iB) = ABatomPairTE
        ENDDO !iB
       ENDDO !iA
      ENDIF
      DO iC=1,nC!setting%redCS%AO(3)%nAtoms
       aoC => setting%redCS%AO(3)%batch(iC)
       globalstartbatchC = aoC%GlobalStartBatchindex -1 
       startD = 1
       IF (sameCD) startD=iC
       DO iD=startD,nD!setting%redCS%AO(4)%nAtoms
        aoD => setting%redCS%AO(4)%batch(iD)
        globalstartbatchD = aoD%GlobalStartBatchindex - 1
        CDatomPairTE = 0E0_realk
        DO ibatchD=1,aoD%nBatches
         globalbatchD = globalstartbatchD+ibatchD
         DO ibatchC=1,aoC%nBatches
          maxgab = RHSGAB(globalstartbatchC+ibatchC,globalbatchD)
          IF (maxgab.NE.shortzero) THEN
           IF (RHSDMATset) THEN
              csProduct = maxgab + RHSDENS(globalstartbatchC+ibatchC,globalbatchD)
           ELSE
              csProduct = maxgab
           ENDIF
           screen = (csProduct.LE.RHS_CSTHR)
           IF (.NOT.screen) THEN
!*******************************************************************************************
!      Time estimate (in terms of # of operations, not actual time) for the 
!      four-center two-electron AObatch-quadruplet integrals
!*******************************************************************************************
! ToDo 1. refine and test
! ToDo 2. account for primitive screening, maybe some nPrimAB*nPrimCD
             batchTE = aoc%nPrim(ibatchc)*aoD%nPrim(iBatchD)&
     &                  *(aoC%maxAng(iBatchC)+aoD%maxAng(iBatchD)+10E0_realk)**POWER
             CDatomPairTE = CDatomPairTE + batchTE
           ENDIF !CS-screening
          ENDIF !RHSBATCH%maxgab.NE.shortzero
         ENDDO !ibatchC
        ENDDO !ibatchD
        tmRHS(iC,iD) = CDatomPairTE
       ENDDO !iD
      ENDDO !iC
 ELSE
! TO WHOEVER WROTE THIS CODE (Simen I guess) THIS IS A CODE FOR NOT
! USING screening - then why the fuck do you use setting%redCS
    DO iA=1,nA
!       aoA => setting%redCS%AO(1)%batch(iA)
       startB = 1
       IF (sameAB) startB=iA
       DO iB=startB,nB
!          aoB => setting%redCS%AO(2)%batch(iB)
!          ABatomPairTE = 0_realk
!          DO ibatchA=1,aoA%nBatches
!             DO ibatchB=1,aoB%nBatches
!                batchTE = aoA%nPrim(ibatchA)*aoB%nPrim(iBatchB)&
!     &               *(aoA%maxAng(iBatchA)+aoB%maxAng(iBatchB)+10E0_realk)**POWER
!                ABatomPairTE = ABatomPairTE + batchTE
!             ENDDO !ibatchB
!          ENDDO !ibatchA
!          tmLHS(iA,iB) = ABatomPairTE
          tmLHS(iA,iB) = 1.0E0_realk
       ENDDO !iB
    ENDDO !iA    
    DO iC=1,nC
!       aoC => setting%redCS%AO(3)%batch(iC)
       startD = 1
       IF (sameCD) startD=iC
       DO iD=startD,nD
!          aoD => setting%redCS%AO(4)%batch(iD)
!          CDatomPairTE = 0_realk
!          DO ibatchC=1,aoC%nBatches
!             DO ibatchD=1,aoD%nBatches
!                !*******************************************************************************************
!                !      Time estimate (in terms of # of operations, not actual time) for the 
!                !      four-center two-electron AObatch-quadruplet integrals
!                !*******************************************************************************************
!                ! ToDo 1. refine and test
!                ! ToDo 2. account for primitive screening, maybe some nPrimAB*nPrimCD
!                batchTE = aoc%nPrim(ibatchc)*aoD%nPrim(iBatchD)&
!                     &     *(aoC%maxAng(iBatchC)+aoD%maxAng(iBatchD)+10E0_realk)**POWER
!                CDatomPairTE = CDatomPairTE + batchTE
!             ENDDO !ibatchD
!          ENDDO !ibatchC
!          tmRHS(iC,iD) = CDatomPairTE
          tmRHS(iC,iD) = 1.0E0_realk
       ENDDO !iD
    ENDDO !iC
 ENDIF
#else
!    Then add the batch-wise time-estimates to the atomic pair estimates
!!$     DO iA=1,nA
!!$      aoA => setting%redCS%AO(1)%batch(iA)
!!$      startB = 1
!!$      IF (sameAB) startB=iA
!!$      DO iB=startB,nB
!!$       indAB = LHSGAB%index(iA,iB)
!!$       IF (indAB.NE.0) THEN
!!$        aoB => setting%redCS%AO(2)%batch(iB)
!!$        ABatomPairTE = 0_realk
!!$        DO ibatchA=1,aoA%nBatches
!!$         DO ibatchB=1,aoB%nBatches
!!$          IF (LHSGAB%LSAO(indAB)%alloc) THEN
!!$           LHSBATCH => LHSGAB%LSAO(indAB)%BATCH(ibatchA,iBatchB)
!!$           IF (LHSBATCH%maxgab.NE.shortzero) THEN
!!$            IF (LHSDMATset) THEN
!!$              LHSDENSB => setting%redCS%LHSDMAT%LSAO(indAB)%BATCH(ibatchA,iBatchB)
!!$              csProduct = LHSBATCH%maxgab + LHSDENSB%maxgab
!!$            ELSE
!!$              csProduct = LHSBATCH%maxgab
!!$            ENDIF
!!$            screen = (csProduct.LT.LHS_CSTHR)
!!$            IF (.not.screen) THEN
!!$             outerCSproduct = csProduct
!!$             IF (outerCSproduct.LT.shortzero) outerCSproduct=shortzero
!!$             DO iC=1,nC
!!$              aoC => setting%redCS%AO(3)%batch(iC)
!!$              startD = 1
!!$              IF (sameCD) startD=iC
!!$              DO iD=startD,nD
!!$               indCD = RHSGAB%index(iC,iD)
!!$               IF (indCD.NE.0) THEN
!!$                aoD => setting%redCS%AO(4)%batch(iD)
!!$                CDatomPairTE = 0_realk
!!$                DO ibatchC=1,aoC%nBatches
!!$                 DO ibatchD=1,aoD%nBatches
!!$                  IF (RHSGAB%LSAO(indCD)%alloc) THEN
!!$                   RHSBATCH => RHSGAB%LSAO(indCD)%BATCH(ibatchC,iBatchD)
!!$                   IF (RHSBATCH%maxgab.NE.shortzero) THEN
!!$                    IF (RHSDMATset) THEN
!!$                     RHSDENSB => setting%redCS%RHSDMAT%LSAO(indCD)%BATCH(ibatchC,iBatchD)
!!$                     csProduct = outerCSproduct + RHSBATCH%maxgab + RHSDENSB%maxgab
!!$                    ELSE
!!$                     csProduct = outerCSproduct + RHSBATCH%maxgab
!!$                    ENDIF
!!$                    screen = (csProduct.LE.CS_THRLOG)
!!$                    IF (.NOT.screen) THEN
!!$!*******************************************************************************************
!!$!      Time estimate (in terms of # of operations, not actual time) for the 
!!$!      four-center two-electron AObatch-quadruplet integrals
!!$!*******************************************************************************************
!!$! ToDo 1. refine and test
!!$! ToDo 2. account for primitive screening, maybe some nPrimAB*nPrimCD
!!$                      batchTE = aoA%nPrim(ibatchA)*aoB%nPrim(iBatchB)*&
!!$     &                          aoc%nPrim(ibatchc)*aoD%nPrim(iBatchD)*&
!!$     &                          (aoA%maxAng(iBatchA)+aoB%maxAng(iBatchB)+&
!!$     &                           aoC%maxAng(iBatchC)+aoD%maxAng(iBatchD)+4)**POWER
!!$                      ABatomPairTE = ABatomPairTE + batchTE
!!$                      CDatomPairTE = CDatomPairTE + batchTE
!!$                    ENDIF !CS-screening
!!$                   ENDIF !RHSBATCH%maxgab.NE.shortzero
!!$                  ENDIF !RHS-alloc
!!$                 ENDDO !ibatchD
!!$                ENDDO !ibatchC
!!$                tmRHS(iC,iD) = tmRHS(iC,iD) + CDatomPairTE
!!$               ENDIF !indCD.NE.0
!!$              ENDDO !iD
!!$             ENDDO !iC
!!$            ENDIF !screen
!!$           ENDIF !maxgab.NE.shortzero
!!$          ENDIF !LHS-alloc
!!$        ENDDO !ibatchB
!!$       ENDDO !ibatchA
!!$       tmLHS(iA,iB) = ABatomPairTE
!!$      ENDIF !indAB.NE.0
!!$     ENDDO !iB
!!$    ENDDO !iA
#endif

#if 0
! May be useful when there is no CS-screening matrices available
    IF (noA.AND.noB) THEN
      DO iB=1,nB
        DO iA=1,nA
          tmLHS(iA,iB) = ONE
        ENDDO
      ENDDO
    ELSE IF (noA) THEN
      DO iB=1,nB
        DO iA=1,nA
          tmLHS(iA,iB) = mol(2)%p%atom(iB)%charge
        ENDDO
      ENDDO
    ELSE IF (noB) THEN
      DO iA=1,nA
        DO iB=1,nB
          tmLHS(iA,iB) = mol(1)%p%atom(iA)%charge
        ENDDO
      ENDDO
    ELSE
      DO iA=1,nA
        DO iB=iA,nB
          tmp = mol(1)%p%atom(iA)%center(1)-mol(2)%p%atom(iB)%center(1)
          R2 = tmp*tmp
          tmp = mol(1)%p%atom(iA)%center(2)-mol(2)%p%atom(iB)%center(2)
          R2 = R2+tmp*tmp
          tmp = mol(1)%p%atom(iA)%center(3)-mol(2)%p%atom(iB)%center(3)
          R2 = R2+tmp*tmp
          IF (R2.LT. 1E3_realk) THEN
            R2 = R2+ONE
            tmLHS(iA,iB) = mol(1)%p%atom(iA)%charge*mol(2)%p%atom(iB)%charge/R2
            tmLHS(iB,iA) = mol(1)%p%atom(iA)%charge*mol(2)%p%atom(iB)%charge/R2
          ELSE
            tmLHS(iA,iB) = ZERO
            tmLHS(iB,iA) = ZERO
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    IF (noC.AND.noD) THEN
      DO iD=1,nD
        DO iC=1,nC
          tmRHS(iC,iD) = ONE
        ENDDO
      ENDDO
    ELSE IF (noC) THEN
      DO iD=1,nD
        DO iC=1,nC
          tmRHS(iC,iD) = mol(4)%p%atom(iD)%charge
        ENDDO
      ENDDO
    ELSE IF (noD) THEN
      DO iC=1,nC
        DO iD=1,nD
          tmRHS(iC,iD) = mol(3)%p%atom(iC)%charge
        ENDDO
      ENDDO
    ELSE
      DO iC=1,nC
        DO iD=iC,nD
          tmp = mol(3)%p%atom(iC)%center(1)-mol(4)%p%atom(iD)%center(1)
          R2 = tmp*tmp
          tmp = mol(3)%p%atom(iC)%center(2)-mol(4)%p%atom(iD)%center(2)
          R2 = R2+tmp*tmp
          tmp = mol(3)%p%atom(iC)%center(3)-mol(4)%p%atom(iD)%center(3)
          R2 = R2+tmp*tmp
          IF (R2.LT. 1E3_realk) THEN
            R2 = R2+ONE
            tmRHS(iC,iD) = mol(3)%p%atom(iC)%charge*mol(4)%p%atom(iD)%charge/R2
            tmRHS(iD,iC) = mol(3)%p%atom(iC)%charge*mol(4)%p%atom(iD)%charge/R2
          ELSE
            tmRHS(iC,iD) = ZERO
            tmRHS(iD,iC) = ZERO
          ENDIF
        ENDDO
      ENDDO
    ENDIF

call ls_flshfo(lupri)

DO iA=1,nA
write(*,*) 'tmLHS',AO(1)(1:6),AO(2)(1:6)
write(*,*) 'tmLHS',iA,(tmLHS(ia,ib),ib=1,nB)
ENDDO
DO iC=1,nC
write(*,*) 'tmRHS',AO(3)(1:6),AO(4)(1:6)
write(*,*) 'tmRHS',iC,(tmRHS(iC,iD),iD=1,nD)
ENDDO

#endif

! FIXME UGLY HACK TO DEAL WITH A ZERO timeMatrix 
! but a zero time matrix means that all integrals
! are screened away so we should determine this 
! before this step
sum = 0
DO iB=1,nB
   DO iA=1,nA
      sum = sum + tmLHS(iA,iB) 
   ENDDO
ENDDO
IF(sum.LT.10E-10)THEN
   DO iB=1,nB
      DO iA=1,nA
         tmLHS(iA,iB) = 1
      ENDDO
   ENDDO
ENDIF

sum = 0
DO iD=1,nD
   DO iC=1,nC
      sum = sum + tmRHS(iC,iD)
   ENDDO
ENDDO
IF(sum.LT.10E-10)THEN
   DO iD=1,nD
      DO iC=1,nC
         tmRHS(iC,iD) = 1
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE jengine_time_estimate

SUBROUTINE LinK_time_estimate(tmLHS,tmRHS,nA,nB,nC,nD,sameOD,noA,noB,noC,noD,&
 &                            AO1,AO2,AO3,AO4,Spec,setting,lupri)
use typedef
implicit none
Real(realk),pointer           :: tmLHS(:,:)
Real(realk),pointer           :: tmRHS(:,:)
integer                       :: AO1,AO2,AO3,AO4,Spec
integer,intent(IN)            :: lupri
type(lssetting),intent(inout) :: setting
logical                       :: sameOD,noA,noB,noC,noD
integer                       :: nA,nB,nC,nD
!
Integer                      :: power
integer                      :: AO(4)
type(molecule_pt)            :: mol(4)
!type(basisset_pt)            :: bas(4)
integer                      :: nAtoms(4)
logical                      :: noPart(4)
integer                      :: iA,iB,iC,iD,iAO,indAB,indCD,startC,startD
real(realk)                  :: R2,tmp,sum
integer                      :: nbatchA,nbatchB,nbatchC,nbatchD
integer                      :: ibatchA,ibatchB,ibatchC,ibatchD
real(realk),parameter        :: one=1E0_realk, zero=0E0_realk
real(realk)                  :: L2,PDIST
TYPE(AOBATCHINFO),pointer    :: aoA,aoB,aoC,aoD
integer(kind=short),pointer  :: LHSGAB(:,:),RHSGAB(:,:),LHSDENS(:,:),RHSDENS(:,:)
real(realk)                  :: ACatomPairTE,BDatomPairTE
real(realk)                  :: batchTE
Logical                      :: screen
integer(short)               :: outerCsProduct,csProduct,CS_THRLOG,LHS_CSTHR,RHS_CSTHR
integer(kind=short)          :: maxgab
integer :: globalstartbatchA,globalstartbatchB,globalstartbatchC,globalstartbatchD
integer :: globalbatchB,globalbatchD

!Set up power - used for the time estimate of batch-wise integrals
power = 7
IF (AO1.EQ.AOEmpty) power = power - 1
IF (AO2.EQ.AOEmpty) power = power - 1
IF (AO3.EQ.AOEmpty) power = power - 1
IF (AO4.EQ.AOEmpty) power = power - 1
IF (Spec.EQ.GradientSpec)  power = power + 1

!Set up distance-damping parameters - damping of time estimate according to 1/(1 + R2/L2)*PDIST
L2    = 1E2_realk
PDIST = 4E0_realk

!Set up molecules and basis
mol(1)%p => setting%molecule(1)%p
mol(2)%p => setting%molecule(2)%p
mol(3)%p => setting%molecule(3)%p
mol(4)%p => setting%molecule(4)%p
AO(1) = AO1
AO(2) = AO2
AO(3) = AO3
AO(4) = AO4
DO iAO=1,4
  IF (ASSOCIATED(mol(iAO)%p)) THEN
    nAtoms(iAO) = mol(iAO)%p%nAtomsNPC
  ELSE
    nAtoms(iAO) = 1
  ENDIF
  noPart(iAO) = .FALSE.
  IF (AO(iAO).EQ.AOEmpty) THEN
    nAtoms(iAO) = 1
    noPart(iAO) = .TRUE.
  ELSE IF (AO(iAO).EQ.AOdfAux) THEN
!    bas(iAO)%p => setting%basis(iAO)%p%BINFO(AuxBasParam)
  ELSE IF (AO(iAO).EQ.AORegular) THEN
!    bas(iAO)%p => setting%basis(iAO)%p%BINFO(RegBasParam)
  ELSE IF (AO(iAO).EQ.AOVAL) THEN
!    bas(iAO)%p => setting%basis(iAO)%p%BINFO(ValBasParam)
  ELSE IF (AO(iAO).EQ.AOdfCABS) THEN
     nAtoms(iAO) = 2*mol(iAO)%p%nAtomsNPC
  ELSE IF (AO(iAO).EQ.AOdfCABO) THEN
!    bas(iAO)%p => setting%basis(iAO)%p%BINFO(CABBasParam)
  ELSE IF (AO(iAO).EQ.AOdfJK) THEN
!    bas(iAO)%p => setting%basis(iAO)%p%BINFO(JKBasParam)
  ELSE IF (AO(iAO).EQ.AOadmm) THEN
!    bas(iAO)%p => setting%basis(iAO)%p%BINFO(ADMBasParam)
  ELSE IF (AO(iAO).EQ.AOpCharge) THEN
    nAtoms(iAO) = mol(iAO)%p%nAtoms
    noPart(iAO) = .TRUE.
  ELSE IF (AO(iAO).EQ.AOelField) THEN
    nAtoms(iAO) = mol(iAO)%p%nAtoms
    noPart(iAO) = .TRUE.
  ELSE IF (AO(iAO).EQ.AONuclear) THEN
    nAtoms(iAO) = 1
    noPart(iAO) = .TRUE.
  ELSE
    WRITE(LUPRI,'(1X,2A)') 'Error in jengine_time_estimate:',iAO,AO(iAO)
    CALL LSQUIT('Error in jengine_time_estimate',lupri)
  ENDIF
!     Special case for batchindex only
  if (setting%Batchindex(iao).NE.0) nAtoms(iAO) = 1
ENDDO
nA = nAtoms(1)
nB = nAtoms(2)
nC = nAtoms(3)
nD = nAtoms(4)
noA = noPart(1)
noB = noPart(2)
noC = noPart(3)
noD = noPart(4)

call mem_alloc(tmLHS,nA,nC)
call mem_alloc(tmRHS,nB,nD)

!    Zero out timeMatrices
DO iC=1,nC
  DO iA=1,nA
    tmLHS(iA,iC) = ZERO
  ENDDO
ENDDO
DO iD=1,nD
  DO iB=1,nB
    tmRHS(iB,iD) = ZERO
  ENDDO
ENDDO


!Set up batch-wise time estimates
DO iA=1,nA
   aoA => setting%redCS%AO(1)%batch(iA)
   startC = 1
   IF (sameOD) startC=iA
   DO iC=startC,nC
      aoC => setting%redCS%AO(3)%batch(iC)
      ACatomPairTE = 0_realk
      tmp = mol(1)%p%atom(iA)%center(1)-mol(3)%p%atom(iC)%center(1)
      R2 = tmp*tmp
      tmp = mol(1)%p%atom(iA)%center(2)-mol(3)%p%atom(iC)%center(2)
      R2 = R2+tmp*tmp
      tmp = mol(1)%p%atom(iA)%center(3)-mol(3)%p%atom(iC)%center(3)
      R2 = R2+tmp*tmp
      tmp = (R2/L2+one)**PDIST
      DO ibatchA=1,aoA%nBatches
         DO ibatchC=1,aoC%nBatches
            batchTE = aoA%nPrim(ibatchA)*aoC%nPrim(iBatchC)&
 &               *(aoA%maxAng(iBatchA)+aoC%maxAng(iBatchC)+10E0_realk)**POWER/tmp
            ACatomPairTE = ACatomPairTE + batchTE
         ENDDO !ibatchB
      ENDDO !ibatchA
      tmLHS(iA,iC) = ACatomPairTE
   ENDDO !iB
ENDDO !iA    
DO iB=1,nB
   aoB => setting%redCS%AO(2)%batch(iB)
   startD = 1
   IF (sameOD) startD=iB
   DO iD=startD,nD
      aoD => setting%redCS%AO(4)%batch(iD)
      BDatomPairTE = 0_realk
      tmp = mol(2)%p%atom(iB)%center(1)-mol(4)%p%atom(iD)%center(1)
      R2 = tmp*tmp
      tmp = mol(2)%p%atom(iB)%center(2)-mol(4)%p%atom(iD)%center(2)
      R2 = R2+tmp*tmp
      tmp = mol(2)%p%atom(iB)%center(3)-mol(4)%p%atom(iD)%center(3)
      R2 = R2+tmp*tmp
      tmp = (R2/L2+one)**PDIST
      DO ibatchB=1,aoB%nBatches
         DO ibatchD=1,aoD%nBatches
            !*******************************************************************************************
            !      Time estimate (in terms of # of operations, not actual time) for the 
            !      four-center two-electron AObatch-quadruplet integrals
            !*******************************************************************************************
            ! ToDo 1. refine and test
            ! ToDo 2. account for primitive screening, maybe some nPrimAB*nPrimCD
            batchTE = aoB%nPrim(ibatchB)*aoD%nPrim(iBatchD)&
                 &     *(aoB%maxAng(iBatchB)+aoD%maxAng(iBatchD)+10E0_realk)**POWER/tmp
            BDatomPairTE = BDatomPairTE + batchTE
         ENDDO !ibatchD
      ENDDO !ibatchC
      tmRHS(iB,iD) = BDatomPairTE
   ENDDO !iD
ENDDO !iC

! FIXME UGLY HACK TO DEAL WITH A ZERO timeMatrix 
! but a zero time matrix means that all integrals
! are screened away so we should determine this 
! before this step
sum = 0
DO iC=1,nC
   DO iA=1,nA
      sum = sum + tmLHS(iA,iC) 
   ENDDO
ENDDO
IF(sum.LT.10E-10)THEN
   DO iC=1,nC
      DO iA=1,nA
         tmLHS(iA,iC) = 1
      ENDDO
   ENDDO
ENDIF

sum = 0
DO iD=1,nD
   DO iB=1,nB
      sum = sum + tmRHS(iB,iD)
   ENDDO
ENDDO
IF(sum.LT.10E-10)THEN
   DO iD=1,nD
      DO iB=1,nB
         tmRHS(iB,iD) = 1
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE LinK_time_estimate

SUBROUTINE set_reduced_screen_info(redCS,AO,dmat_lhs,dmat_rhs,lhs_created,rhs_created,&
      &                            CS_rhs,CS_lhs,CS_THRLOG,CS_SCREEN,lupri)
implicit none
TYPE(reducedScreeningInfo),intent(INOUT) :: redCS
TYPE(AOITEMPOINTER),intent(IN)           :: AO(4)
TYPE(lstensor),pointer                   :: dmat_lhs,dmat_rhs
Logical,intent(IN)                       :: lhs_created,rhs_created,CS_SCREEN
TYPE(lstensor),pointer                   :: CS_rhs,CS_lhs
integer(short),intent(IN)                :: CS_THRLOG
Integer,intent(IN)                       :: lupri
!
Integer :: iAO,dim1,dim2

redCS%isset = .TRUE.
redCS%CS_THRLOG = CS_THRLOG
DO iAO=1,4
  call set_aoAtomInfo(redCS%AO(iAO),AO(iAO)%p)
ENDDO
!make sls-tensors
IF(CS_SCREEN)THEN
   redCS%LHSGAB => CS_lhs%maxgab
   redCS%RHSGAB => CS_rhs%maxgab
   redCS%maxgabLHS = CS_lhs%maxgabelm
   redCS%maxgabRHS = CS_rhs%maxgabelm
   redCS%nbatches(1) = CS_lhs%nbatches(1)
   redCS%nbatches(2) = CS_lhs%nbatches(2)
   redCS%nbatches(3) = CS_rhs%nbatches(3)
   redCS%nbatches(4) = CS_rhs%nbatches(4)

   IF (lhs_created) THEN
      redCS%LHSDMATset = .TRUE.
      dim1 = dmat_lhs%nbatches(1)
      dim2 = dmat_lhs%nbatches(2)
      redCS%nbatches(1) = dim1 
      redCS%nbatches(2) = dim2
      call mem_alloc(redCS%LHSDMAT,dim1,dim2)
      call Build_full_shortint2dim_from_lstensor(dmat_lhs,redCS%LHSDMAT,dim1,dim2,lupri)
      call set_lst_maxgabelms(dmat_lhs)
      redCS%maxDmatLHS = dmat_lhs%maxgabelm
   ENDIF
   IF (rhs_created) THEN
      redCS%RHSDMATset = .TRUE.
      dim1 = dmat_rhs%nbatches(1)
      dim2 = dmat_rhs%nbatches(2)
      redCS%nbatches(3) = dim1 
      redCS%nbatches(4) = dim2
      call mem_alloc(redCS%RHSDMAT,dim1,dim2)
      call Build_full_shortint2dim_from_lstensor(dmat_rhs,redCS%RHSDMAT,dim1,dim2,lupri)
      call set_lst_maxgabelms(dmat_rhs)
      redCS%maxDmatRHS = dmat_rhs%maxgabelm
   ENDIF
ENDIF
  
END SUBROUTINE set_reduced_screen_info

SUBROUTINE set_aoAtomInfo(AOinfo,AObatch)
implicit none
TYPE(AOatomInfo),intent(INOUT) :: AOinfo
TYPE(AOitem),intent(IN)        :: AObatch
!
Integer :: ibatch,iatom
Integer,pointer :: nBatches(:)

IF (AObatch%empty) THEN
  !Make empty AOinfo
  call makeEmptyAOatomInfo(AOinfo)
ELSE
 call mem_alloc(nBatches,AObatch%natoms)
 nBatches(1:AObatch%natoms) = 0
 DO ibatch=1,AObatch%nBatches
   iatom = AObatch%batch(ibatch)%atom
   nBatches(iatom) = nBatches(iatom) + 1
 ENDDO
 
 !Set first AOATOMINFO
 AOinfo%nAtoms    = AObatch%natoms
 AOinfo%nTotBatch = AObatch%nbatches
 ALLOCATE(AOinfo%batch(AOinfo%nAtoms))
 
 !Then allocate the AOBATCHINFO
 ibatch = 1
 DO iAtom=1,AOinfo%nAtoms
   AOinfo%batch(iAtom)%nBatches   = nBatches(iAtom)
   AOinfo%batch(iAtom)%GlobalStartBatchindex = ibatch
   ibatch  = ibatch + nBatches(iAtom)
   call mem_alloc(AOinfo%batch(iAtom)%nPrim,nBatches(iAtom))
   call mem_alloc(AOinfo%batch(iAtom)%maxAng,nBatches(iAtom))
 ENDDO
 
 nBatches(1:AObatch%natoms) = 0
 !Then set up the AOBATCHINFO
 DO ibatch=1,AObatch%nBatches
   iatom = AObatch%batch(ibatch)%atom
   nBatches(iatom) = nBatches(iatom) + 1
   AOinfo%batch(iAtom)%nPrim(nBatches(iatom))  = AObatch%batch(ibatch)%nPrimitives
   AOinfo%batch(iAtom)%maxAng(nBatches(iatom)) = AObatch%batch(ibatch)%maxAngmom
 ENDDO
 
  call mem_dealloc(nBatches)
ENDIF
END SUBROUTINE set_aoAtomInfo

SUBROUTINE makeEmptyAOatomInfo(AOinfo)
implicit none
TYPE(AOatomInfo), intent(INOUT) :: AOinfo
!Make empty AOinfo
AOinfo%nAtoms    = 1
AOinfo%nTotBatch = 1
ALLOCATE(AOinfo%batch(1))
CALL makeEmptyAoBatchInfo(AOinfo%batch(1))
END SUBROUTINE makeEmptyAOatomInfo

SUBROUTINE makeEmptyAoBatchInfo(AObatch)
implicit none
TYPE(AoBatchInfo),intent(INOUT) :: AObatch
AObatch%nBatches = 1 
AObatch%GlobalStartBatchindex = 1 
call mem_alloc(AObatch%nPrim,1)
call mem_alloc(AObatch%maxAng,1)
AObatch%nPrim(1)=1
AObatch%maxAng(1)=0
END SUBROUTINE makeEmptyAOBATCHINFO
#endif


!*************************************************************************
!*
!*         Routines for the task-based fragmentation
!*
!*         Simen: to replace the fragments/blocks above for MPI
!*
!*************************************************************************

!> \brief Set up fragments accoring to specifications in tasks
!> \author S. Reine
!> \date 2011-10-18
!> \param Setting The lssetting item - here the fragments are set up
!> \param task The task - defining the framentation
!> \param Side Either the left- or right-hand side ('LHS' or 'RHS')
!> \param lupri The output print unit
!> \param aux Specifies the use of an auxiliary basis
SUBROUTINE SetTaskFragments(Setting,task,Side,aux,cabs,jk,part13,lupri)
implicit none
TYPE(LSSETTING),intent(INOUT) :: Setting
TYPE(LSTASK),intent(IN)       :: task
Integer,intent(IN)            :: lupri
Character*(*)                 :: Side
Logical,intent(in)            :: aux,cabs,jk,part13
!
Integer :: aoA,aoB,AO(4)

IF (part13) THEN
  AO(1) = 1
  AO(2) = 3
  AO(3) = 2
  AO(4) = 4
ELSE
  AO(1) = 1
  AO(2) = 2
  AO(3) = 3
  AO(4) = 4
ENDIF


! Turn off symmetries between LHS and RHS
Setting%sameFrag(AO(1),AO(3)) = .FALSE.
Setting%sameFrag(AO(1),AO(4)) = .FALSE.
Setting%sameFrag(AO(2),AO(3)) = .FALSE.
Setting%sameFrag(AO(2),AO(4)) = .FALSE.
Setting%sameFrag(AO(3),AO(1)) = .FALSE.
Setting%sameFrag(AO(4),AO(1)) = .FALSE.
Setting%sameFrag(AO(3),AO(2)) = .FALSE.
Setting%sameFrag(AO(4),AO(2)) = .FALSE.

Select case(Side)
Case('LHS')
  aoA = AO(1)
  aoB = AO(2)
Case('RHS')
  aoA = AO(3)
  aoB = AO(4)
Case Default
  CALL LSQUIT('Error in SetTaskFragments. Inncorrect Side',lupri)
End select

IF (task%full_mol_row) THEN
  IF (Setting%fragBuild(aoA)) THEN
    call free_Moleculeinfo(Setting%fragment(aoA)%p)
    DEALLOCATE(Setting%fragment(aoA)%p)
    NULLIFY(Setting%fragment(aoA)%p)
    Setting%fragBuild(aoA) = .FALSE.
  ENDIF
  Setting%fragment(aoA)%p => Setting%molecule(aoA)%p
ELSE
  IF (Setting%fragBuild(aoA)) THEN
    call free_Moleculeinfo(Setting%fragment(aoA)%p)
  ELSE
    NULLIFY(Setting%fragment(aoA)%p)
    ALLOCATE(Setting%fragment(aoA)%p)
  ENDIF
  CALL BUILD_FRAGMENT(Setting%molecule(aoA)%p,Setting%fragment(aoA)%p,&
     &                Setting%basis(aoA)%p,task%row_atoms,task%nrow,lupri)
  Setting%fragBuild(aoA) = .TRUE.
ENDIF

IF (task%full_mol_col) THEN
  IF (Setting%fragBuild(aoB)) THEN
    call free_Moleculeinfo(Setting%fragment(aoB)%p)
    DEALLOCATE(Setting%fragment(aoB)%p)
    NULLIFY(Setting%fragment(aoB)%p)
    Setting%fragBuild(aoB) = .FALSE.
  ENDIF
  Setting%fragment(aoB)%p => Setting%molecule(aoB)%p
ELSE
  IF (Setting%fragBuild(aoB)) THEN
    call free_Moleculeinfo(Setting%fragment(aoB)%p)
  ELSE
    NULLIFY(Setting%fragment(aoB)%p)
    ALLOCATE(Setting%fragment(aoB)%p)
  ENDIF
  CALL BUILD_FRAGMENT(Setting%molecule(aoB)%p,Setting%fragment(aoB)%p,&
     &                Setting%basis(aoB)%p,task%col_atoms,task%ncol,lupri)
  Setting%fragBuild(aoB) = .TRUE.
ENDIF

IF (task%full_mol_row.AND.task%full_mol_col) THEN
! Symmetry
  Setting%sameFrag(aoA,aoB) = Setting%sameMol(aoA,aoB)
  Setting%sameFrag(aoB,aoA) = Setting%sameMol(aoB,aoA)
ELSE
! No symmetry
  Setting%sameFrag(aoA,aoB) = .FALSE.
  Setting%sameFrag(aoB,aoA) = .FALSE.
ENDIF

END SUBROUTINE SetTaskFragments

!> \brief Returns the number of basis funcitons for the two fragments of a task
!> \author S. Rein
!> \date 2011-10-18
!> \param orbInfo The orbital information for each molecule
!> \param task The task - defining the framentation
!> \param Side Either the left- or right-hand side ('LHS' or 'RHS')
!> \param Spec Specifying the regular or differentiated integrals ('Regular','Gradient', etc.', etc.)
!> \param lupri The output print unit
!> \param aux Specifies the use of an auxiliary basis
SUBROUTINE getTaskDimension(orbInfo,task,AOA,AOB,Side,Spec,nbastA,nbastB,lupri)
implicit none
TYPE(MolecularOrbitalInfo),intent(IN) :: orbInfo(4)
TYPE(LSTASK),intent(IN)               :: task
Integer,intent(IN)                    :: lupri
integer                               :: AOA,AOB,Spec
Character*(*)                         :: Side
Integer,intent(OUT)                   :: nbastA,nbastB
!
Integer :: iA,iB,iatom,iatomfull

Select case(Side)
Case('LHS')
  iA = 1
  iB = 2
Case('RHS')
  iA = 3
  iB = 4
Case Default
  CALL LSQUIT('Error in getTaskDimension. Inncorrect Side',lupri)
End select

IF ((AOA.EQ.AOEmpty).OR.(AOA.EQ.AONuclear)) THEN
  nbastA = 1
ELSE
  IF (Spec.EQ.GradientSpec) THEN
    nbastA = task%nrow
  ELSE
    nbastA = 0
    DO iatom=1,task%nrow
      IF (task%full_mol_row) THEN
        iatomfull = iatom
      ELSE
        iatomfull = task%row_atoms(iatom)
      ENDIF
      IF ((AOA.EQ.AOdfAux).OR.(AOA.EQ.AOdfCABO).OR.(AOA.EQ.AOdfJK).OR.(AOA.EQ.AOadmm)) THEN
        nbastA = nbastA + orbInfo(iA)%numAtomicOrbitalsAux(iatomfull)
     ELSEIF (AOA.EQ.AOdfCABS) THEN
        IF(iatomfull.GT.size(orbInfo(iA)%numAtomicOrbitalsReg))THEN
           nbastA = nbastA + &
 &orbInfo(iA)%numAtomicOrbitalsAux(iatomfull-size(orbInfo(iA)%numAtomicOrbitalsReg))
        ELSE
           nbastA = nbastA + orbInfo(iA)%numAtomicOrbitalsReg(iatomfull)
        ENDIF
      ELSE IF ((AOA.EQ.AORegular).OR.(AOA.EQ.AOVAL)) THEN
        nbastA = nbastA + orbInfo(iA)%numAtomicOrbitalsReg(iatomfull)
      ELSE IF (AOA.EQ.AOpCharge) THEN
        nbastA = nbastA + 1
      ELSE IF (AOA.EQ.AOelField) THEN
        nbastA = nbastA + 3
      ELSE
        CALL LSQUIT('Error in getTaskDimension. Not an implemented AOA option',lupri)
      ENDIF
    ENDDO
  ENDIF
ENDIF

IF ((AOB.EQ.AOEmpty).OR.(AOB.EQ.AONuclear)) THEN
  nBastB = 1
ELSE
  IF (Spec.EQ.GradientSpec) THEN
    nbastB = task%ncol
  ELSE
    nBastB = 0
    DO iatom=1,task%ncol
      IF (task%full_mol_col) THEN
        iatomfull = iatom
      ELSE
        iatomfull = task%col_atoms(iatom)
      ENDIF
      IF ((AOB.EQ.AORegular).OR.(AOB.EQ.AOVAL)) THEN
        nBastB = nBastB + orbInfo(iB)%numAtomicOrbitalsReg(iatomfull)
      ELSE IF ((AOB.EQ.AOdfAux).OR.(AOB.EQ.AOdfCABO).OR.(AOB.EQ.AOdfJK).OR.(AOB.EQ.AOadmm)) THEN
        nBastB = nBastB + orbInfo(iB)%numAtomicOrbitalsAux(iatomfull)
     ELSEIF (AOB.EQ.AOdfCABS) THEN
        IF(iatomfull.GT.size(orbInfo(iB)%numAtomicOrbitalsReg))THEN
           nbastB = nbastB + &
& orbInfo(iB)%numAtomicOrbitalsAux(iatomfull-size(orbInfo(iB)%numAtomicOrbitalsReg))
        ELSE
           nbastB = nbastB + orbInfo(iB)%numAtomicOrbitalsReg(iatomfull)
        ENDIF
      ELSE IF (AOB.EQ.AOpCharge) THEN
        nbastB = nbastB + 1
      ELSE IF (AOB.EQ.AOelField) THEN
        nbastB = nbastB + 3
      ELSE
        CALL LSQUIT('Error in getTaskDimension. Not an implemented AOB option',lupri)
      ENDIF
    ENDDO
  ENDIF
ENDIF

END SUBROUTINE getTaskDimension

!> \brief determine the symmetries: Is both the LHS AOs the same? RHS? are both Overlap distributions the same
!> \author S. Reine
!> \date 2010
!> \param sameAOsLHS if the LHS AOS are the same
!> \param sameAOsRHS if the RHS AOS are the same
!> \param sameODs if both ODs are the same
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
SUBROUTINE GetSymmetries(sameAOsLHS,sameAOsRHS,sameODs,AO1,AO2,AO3,AO4,sameMol,nAO)
implicit none
Integer,intent(IN)  :: nAO
Logical,intent(OUT) :: sameAOsLHS,sameAOsRHS,sameODs
Logical,intent(IN)  :: sameMol(nAO,nAO)
integer       :: AO1,AO2,AO3,AO4

sameAOsLHS = (AO1.EQ.AO2).AND..NOT.(AO1.EQ.AOEmpty).AND.sameMol(1,2)
sameAOsRHS = (AO3.EQ.AO4).AND..NOT.(AO3.EQ.AOEmpty).AND.sameMol(3,4)
sameODs    = (AO1.EQ.AO3) .AND. (AO2.EQ.AO4).AND.sameMol(1,3).AND.sameMol(2,4)

END SUBROUTINE GetSymmetries

FUNCTION getNAtomsSide(AOA,AOB,symmetry,molA,molB)
implicit none
integer                 :: AOA,AOB
Logical,intent(IN)            :: symmetry
Type(moleculeinfo),intent(IN) :: molA,molB
integer                       :: getNAtomsSide
!
Logical :: emptyA,emptyB
integer :: factor,factor2

emptyA = AOA.EQ.AOEmpty
emptyB = AOB.EQ.AOEmpty
IF (AOA.EQ.AOdfCABS.AND.AOB.EQ.AOdfCABS)THEN 
   factor = 4;   factor2 = 2
ELSEIF (AOA.EQ.AOdfCABS)THEN 
   factor = 2;   factor2 = 2
ELSEIF (AOB.EQ.AOdfCABS)THEN 
   factor = 2;   factor2 = 2
ELSE
   factor = 1;   factor2 = 1  
ENDIF

IF (emptyA.AND.emptyB) THEN
  getNAtomsSide = 1
ELSE IF (emptyA) THEN
  getNAtomsSide = factor*molB%nAtomsNPC
ELSE IF (emptyB) THEN
  getNAtomsSide = factor*molA%nAtomsNPC
ELSE
  getNAtomsSide = factor*molA%nAtomsNPC*molB%nAtomsNPC
  IF (symmetry) THEN
    IF (molA%nAtomsNPC.NE.molB%nAtomsNPC) &
    & CALL LSQUIT('Error in getNAtomsSide: symmetry requested but nAtoms not equal',-1)
    getNAtomsSide = factor2*molA%nAtomsNPC*(factor2*molA%nAtomsNPC+1)/2
  ENDIF
ENDIF
END FUNCTION getNAtomsSide

SUBROUTINE lsmpi_time_harvest(mynum,cpu0,wall0,cpu,wall,cpu1,wall1,cpu2,wall2,LHSpart,RHSpart,numnodes,lupri,iprint)
implicit none
Integer(kind=ls_mpik),intent(IN) :: mynum,numnodes
Integer,intent(IN)               :: lupri,iprint
Real(realk),intent(INOUT)        :: cpu,wall,LHSpart,RHSpart,cpu0,wall0,cpu1,wall1,cpu2,wall2
!
#ifdef VAR_MPI
Integer             :: inum, nodes, i,j
Real(realk),pointer :: time_harvester(:,:)
Real(realk)         :: tot_wall
Real(realk)         :: tot_cpu
Real(realk)         :: tot_LHS
Real(realk)         :: tot_RHS

IF (iprint .GT. 0) THEN
  nodes = numnodes
  call mem_alloc(time_harvester,7,nodes)
  call ls_dzero(time_harvester,7*nodes)
  inum = mynum
  IF (inum.EQ.0) inum=nodes
  time_harvester(1,inum) = cpu
  time_harvester(2,inum) = wall
  time_harvester(3,inum) = LHSpart
  time_harvester(4,inum) = RHSpart
  time_harvester(5,inum) = cpu0
  time_harvester(6,inum) = cpu1
  time_harvester(7,inum) = cpu2
  CALL lsmpi_reduction(time_harvester,7,nodes,infpar%master,MPI_COMM_LSDALTON)
  IF (mynum .EQ. infpar%master) THEN
    tot_cpu  = 0E0_realk
    tot_wall = 0E0_realk
    tot_LHS  = 0E0_realk
    tot_RHS  = 0E0_realk
  ! Retrive timings and expected time-partitions from each slave
    DO i=1,nodes
      tot_cpu  = tot_cpu  + time_harvester(1,i)
      tot_wall = tot_wall + time_harvester(2,i)
      tot_LHS  = tot_LHS  + time_harvester(3,i)
      tot_RHS  = tot_RHS  + time_harvester(4,i)
    ENDDO
    write(lupri,'(1X,A11,8A14)') '       Node', '      CPU time','  CPU fraction', ' wall fraction', &
     & '  LHS fraction', '  RHS fraction','      CPU task','     CPU final','      CPU init'
    DO i=0,nodes-1
      inum = i
      IF (inum.EQ.0) inum=nodes
      time_harvester(1,inum) = time_harvester(1,inum)/tot_cpu*nodes
      IF (tot_wall.GT.0E0_realk) &
      & time_harvester(2,inum) = time_harvester(2,inum)/tot_wall*nodes
      time_harvester(3,inum) = time_harvester(3,inum)!/tot_LHS*nodes
      time_harvester(4,inum) = time_harvester(4,inum)!/tot_RHS*nodes
      write(lupri,'(3X,I7,8F14.4)') i,time_harvester(1,inum)/nodes*tot_cpu,(time_harvester(j,inum),j=1,7)
    ENDDO
  ENDIF
  call mem_dealloc(time_harvester)
ENDIF
#endif
END SUBROUTINE lsmpi_time_harvest

  end module lsmpi_mod
! ****************      END OF MODULE 
