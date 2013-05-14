!> @file 
!> Contains info about SCF averaging and 'old' queue operations used for DIIS, and EDIIS. 

!> \brief Queue handling routines.
!> \author L. Thogersen. Documented by S. Host.
!> \date 2003
!>
!>  Currently, we only support in-memory queues but we could easily
!>  imagine an alternative file-based implementation.
!>
MODULE av_utilities
   use files
   USE Matrix_operations

   !> \brief Stores Fock/density mats for DIIS and EDIIS
   !> \author L. Thogersen
   !> \date 2003
   !>  
   !>  util_HistoryStore stores data to be preserved between iterations.
   !>  As opposed to the metric below, HistoryStore can have HUGE storage
   !>  requirements and should support an operational mode when the data is
   !>  saved to external storage (tape, disk, etc :).
   !> 
   TYPE util_HistoryStore
      TYPE(Matrix), pointer :: D(:)   !MatrixPointer, max_history_size
      TYPE(Matrix), pointer :: F(:)   !MatrixPointer
      TYPE(Matrix), pointer :: gd(:)  !MatrixPointer
      real(realk),  pointer  :: energy(:)
      real(realk),  pointer  :: metric(:,:)
      real(realk),  pointer  :: DIIS_Amat(:,:)
      INTEGER             :: current_position, used_entries, allocated
   END TYPE util_HistoryStore

!> \brief Contains parameters used for SCF averaging
!> \author L. Thogersen
!> \date 2003
type cfg_set
   !> Smallest accepted overlap when using cfg_lshift = CFG_lshift_MOchange
   real(realk) :: min_density_overlap
   !> How many old Fock/KS and density matrices should be saved/used for averaging?
   integer     :: max_history_size
   !> Largest accepted ratio when using cfg_lshift = cfg_lshift_dorth
   real(realk) :: max_dorth_ratio
end type cfg_set

!> \author S. Host
!> \date March 2010
!> \brief Contains SCF averaging parameters
type AvItem
   !CONFIGURATION:
   !==============
      !> Logical unit number for LSDALTON.OUT
      integer :: lupri
      !> Defines what type of SCF averaging should be used    
      integer :: CFG_averaging
      !AIX compiler does not allow values to be set here - instead, set them in av_set_default_config below
      !> Use no averaging
      integer :: CFG_AVG_NONE ! = 1  
#ifdef VAR_DSM
      !> Use Density Subspace Minimization
      integer :: CFG_AVG_DSM  ! = 2
#endif
      !> Use standard DIIS
      integer :: CFG_AVG_DIIS ! = 3
      !Use E-DIIS
      integer :: CFG_AVG_EDIIS != 4
      !Use Van Lenthe scheme with custom shifts etc
      integer :: CFG_AVG_van_lenthe != 5

#ifdef VAR_DSM
      !DSM SETTINGS:
      !=============
      ! If cfg_averagning = cfg_avg_DSM, several approaches are available
      integer :: cfg_dsm_app
      !> Use DSM steps until limit or convergence
      integer :: cfg_dsm_default  ! = 1
      !> Use only one DSM step
      integer :: cfg_dsm_one      ! = 2
      !> Use line search in the first dsm-step - expensive!
      integer :: cfg_dsm_search   ! = 3
      !> Add the extra more expensive term
      integer :: cfg_dsm_xtra_term != 4    
      !> How many vectors should be saved for averaging?
#endif
      integer :: max_history_size
#ifdef VAR_DSM
      !> How many vectors should be used for DSM?
      integer :: dsm_history_size 
#endif
      !> How many vectors should be used for DIIS?
      integer :: diis_history_size 
      !> How many vectors should be used for E-DIIS?
      integer :: ediis_history_size
      !Things can be run in "safe"-mode where the Fock-matrix returned is the one 
      !evaluated from Dbar = sum_i c_i D_i and not the linear combination of pre-
      !vious Fock-matrices Fbar = sum_i c_i F_i. This should be the same in HF
      !but differ in DFT
      !> Should SCF be run in safe mode?
      logical :: cfg_safe
      !> Some thresholds depend on the calculation type HF v. DFT
      integer :: CFG_THR_HartreeFock != 1
      !> Some thresholds depend on the calculation type HF v. DFT
      integer :: CFG_THR_dft != 2
      !> Some thresholds that depend on the calculation type HF v. DFT
      TYPE(CFG_set) :: cfg_settings(2)! = (/ &
!           & CFG_set( 0.98E0_realk, 10, 0.08E0_realk),  &  ! Hartree-Fock thresholds
!           & CFG_set( 0.975E0_realk, 7, 0.03E0_realk) /)   ! DFT Thresholds 
      !> Is the calculation HF or DFT?
      INTEGER :: CFG_SET_type 
      !> Should a flush of the queue be possible?
      LOGICAL :: cfg_flush_vec
      !> Counts 3 iterations before starting DIIS, after we have started gathering info for DIIS
      integer :: vanlentheCounter
      !> Should the queue be used for averaging (referenced for Van Lenthe Scheme only)?
      logical :: usequeue
      !> Should gradient be saved in addition to Fock and density matrices?
      logical :: save_gradient
      !> True when DIIS is used for Grand Canonical SCF (trilevel)
      logical :: trilevel_gcscf
      !> If DIIS fails for GCSCF (trilevel), just skip it instead of quitting
      logical :: trilevel_gcscf_skip

      !LEVEL SHIFT:
      !============
      !We have to keep track of level shift type, both here and diag structure.
      !Not pretty, but I don't want to rewrite old and semiobsolete code to make
      !it prettier! /Stinne
      !> What kind of level shift should be used?
      integer :: CFG_lshift
      !> Level shift found from min overlap of every single new occ MO on old occ MO space
!      integer :: CFG_lshift_MOchange != 1
      !> A search is made in Escf(mu) - expensive!
      integer :: CFG_lshift_search != 2
      !> Use ratio ||Dorth||/||D|| to find mu
      integer :: cfg_lshift_dorth != 3
      !> Do no level shifting
      integer :: CFG_lshift_none  != 4
      !> Van Lenthe levelshift (see J.Comput.Chem. 27, 926-932 (2005))
      integer :: CFG_lshift_vanlenthe  != 5
   !SETTINGS:
   !=========
#ifdef VAR_DSM
      !> Which entry in queue is the last for DSM?
      integer :: dsm_pos
#endif
      !> Which entry in queue is the last for DIIS?
      integer :: diis_pos
      !> Which entry in queue is the last for EDIIS?
      integer :: ediis_pos
      !> Has queue been flushed?
      logical :: flushed
   !INFO:
   !=====
      logical :: INFO_D_proj
      logical :: INFO_DIIS
#ifdef VAR_DSM
      logical :: INFO_DSM_CNORM_MU_FIG
      logical :: INFO_DSM_EIGENVAL
      logical :: INFO_DSM_ENERGY
      logical :: INFO_DSM_EQ
      logical :: INFO_DSM_EXIT
      logical :: INFO_DSM_GRAD
      logical :: INFO_DSM_METRIC
      logical :: INFO_DSM_NIT
      logical :: INFO_DSM_PROJ
      logical :: INFO_DSM_RATIO
      logical :: INFO_DSM_STEP
      logical :: INFO_DSM_STEP_BRACKET
      logical :: INFO_DSM_STEP_TOTAL
      logical :: INFO_DSM_TRUSTR
#endif
      logical :: INFO_EDIIS
      logical :: INFO_WEIGHT_FINAL
      logical :: INFO_WEIGHTS

   !DEBUGGING:
   !==========
#ifdef VAR_DSM
      logical :: DEBUG_DSM_DCHANGE
      logical :: debug_dsm_dhistory
      logical :: debug_dsm_Ecomp_fig
      logical :: DEBUG_DSM_EMODEL
      logical :: DEBUG_DSM_LINESEARCH
      logical :: DEBUG_DSM_metric
#endif
      logical :: DEBUG_EDIIS
      logical :: DEBUG_RH_MU_E
   !DATA:
   !=====
      !> For simple Van lenthe averaging, where only last F and D are needed
      type(matrix) :: Fprev
      !> For simple Van lenthe averaging, where only last F and D are needed
      type(matrix) :: Dprev
end type AvItem


CONTAINS
!> \brief Set default configuration for SCF averaging.
!> \author S. Host
!> \date March 2010
!> \param av Used to store info about SCF averaging
subroutine av_set_default_config(av)
implicit none
   type(AvItem), intent(inout) :: av

      av%CFG_AVG_NONE       = 1
#ifdef VAR_DSM
      av%CFG_AVG_DSM        = 2
#endif
      av%CFG_AVG_DIIS       = 3
      av%CFG_AVG_EDIIS      = 4
      av%CFG_AVG_van_lenthe = 5

#ifdef VAR_DSM
      av%cfg_dsm_default    = 1
      av%cfg_dsm_one        = 2
      av%cfg_dsm_search     = 3
      av%cfg_dsm_xtra_term  = 4
#endif

      av%CFG_THR_HartreeFock = 1
      av%CFG_THR_dft         = 2

!      av%CFG_lshift_MOchange  = 1
      av%CFG_lshift_search    = 2
      av%cfg_lshift_dorth     = 3
      av%CFG_lshift_none      = 4
      av%CFG_lshift_vanlenthe = 5

      av%cfg_averaging = av%cfg_avg_none        !Default is no averaging, since default optimization is ARH
#ifdef VAR_DSM
      av%cfg_dsm_app   = av%cfg_dsm_default     !Default = standard dsm
#endif
      av%CFG_lshift    = av%CFG_lshift_none     !Default = no level shift
      av%cfg_safe      = .false.
      av%CFG_SET_type  = av%CFG_THR_HartreeFock ! Default is HF
      av%cfg_flush_vec = .false.
      av%max_history_size   = 60 

      av%cfg_settings(1)%min_density_overlap = 0.98E0_realk
      av%cfg_settings(1)%max_history_size = 10
      av%cfg_settings(1)%max_dorth_ratio = 0.08E0_realk
      av%cfg_settings(2)%min_density_overlap = 0.975E0_realk
      av%cfg_settings(2)%max_history_size = 7
      av%cfg_settings(2)%max_dorth_ratio = 0.03E0_realk

#ifdef VAR_DSM
      av%dsm_history_size   = av%cfg_settings(av%CFG_SET_type)%max_history_size
      av%diis_history_size  = av%dsm_history_size
      av%ediis_history_size = av%dsm_history_size
#else
      av%diis_history_size   = av%cfg_settings(av%CFG_SET_type)%max_history_size
      av%ediis_history_size = av%diis_history_size
#endif
      av%usequeue           = .false.
      av%vanlentheCounter   = 0
      av%save_gradient      = .false.
      av%trilevel_gcscf     = .false.
      av%trilevel_gcscf_skip = .false.
      av%flushed            = .false.
   !INFO:
   !=====
      av%INFO_D_proj           = .false.
      av%INFO_DIIS             = .false.
#ifdef VAR_DSM
      av%INFO_DSM_CNORM_MU_FIG = .false.
      av%INFO_DSM_EIGENVAL     = .false.
      av%INFO_DSM_ENERGY       = .false.
      av%INFO_DSM_EQ           = .false.
      av%INFO_DSM_EXIT         = .false.
      av%INFO_DSM_GRAD         = .false.
      av%INFO_DSM_METRIC       = .false.
      av%INFO_DSM_NIT          = .false.
      av%INFO_DSM_PROJ         = .false.
      av%INFO_DSM_RATIO        = .false.
      av%INFO_DSM_STEP         = .false.
      av%INFO_DSM_STEP_BRACKET = .false.
      av%INFO_DSM_STEP_TOTAL   = .false.
      av%INFO_DSM_TRUSTR       = .false.
#endif
      av%INFO_EDIIS            = .false.
      av%INFO_WEIGHT_FINAL     = .false.
      av%INFO_WEIGHTS          = .false.
   !DEBUGGING:
   !=========
#ifdef VAR_DSM
      av%DEBUG_DSM_DCHANGE    = .false.
      av%debug_dsm_dhistory   = .false.
      av%debug_dsm_Ecomp_fig  = .false.
      av%DEBUG_DSM_EMODEL     = .false.
      av%DEBUG_DSM_LINESEARCH = .false.
      av%DEBUG_DSM_metric     = .false.
#endif
      av%DEBUG_EDIIS          = .false.
      av%DEBUG_RH_MU_E        = .false.

end subroutine av_set_default_config

!> \brief Free matrices used for SCF averaging.
!> \author S. Host
!> \date March 2010
!> \param av Used to store info about SCF averaging
subroutine av_shutdown(av)
implicit none
   type(AvItem), intent(inout) :: av

   if (av%CFG_averaging == av%CFG_AVG_van_lenthe) then
      call mat_free(av%Fprev)
      call mat_free(av%Dprev)
   endif
end subroutine av_shutdown

   !> \brief Initialize queue.
   !> \author L. Thogersen
   !> \date 2003
   SUBROUTINE queue_init(av,queue)
      implicit none
      !> Used to store info about SCF averaging
      type(avItem),intent(inout) :: av
      !> Subspace of Fock/KS and density matrices from previous SCF iterations 
      TYPE(util_HistoryStore),intent(inout) :: queue
      integer :: i
      
      nullify(queue%D)
      nullify(queue%F)
      nullify(queue%gd)
      nullify(queue%energy)
      nullify(queue%metric)
      nullify(queue%DIIS_Amat)
      allocate(queue%D(av%max_history_size))
      allocate(queue%F(av%max_history_size))
      allocate(queue%gd(av%max_history_size))
      allocate(queue%energy(av%max_history_size))
      allocate(queue%metric(av%max_history_size,av%max_history_size))
      allocate(queue%DIIS_Amat(av%max_history_size,av%max_history_size))
      
      DO I = 1, av%max_history_size
         !NULLIFY(queue%D(i)%p, queue%F(I)%p, queue%gd(I)%p)
         !to use as allocation check later - see queue_free
         queue%gd(I)%nrow = 0
      END DO
      queue%current_position = 0
      queue%metric = 0.0E0_realk !Stinne 24/10-06
!initialize specific pointers for dsm and diis
#ifdef VAR_DSM
      av%dsm_pos = 0
#endif
      av%diis_pos = 0
      av%ediis_pos = 0
      queue%used_entries = 0
      queue%allocated = 0
      if (av%CFG_averaging > 2) then 
        !CFG_AVG_DIIS = 3, CFG_AVG_EDIIS = 4
!FIXME: should this be done for EDIIS???
        av%save_gradient = .true.
      endif
   END SUBROUTINE queue_init

   !> \brief Free queue.
   !> \author L. Thogersen
   !> \date 2003
   SUBROUTINE queue_free(av,queue)
      implicit none
      !> Used to store info about SCF averaging
      type(avItem),intent(inout) :: av
      !> Subspace of Fock/KS and density matrices from previous SCF iterations 
      TYPE(util_HistoryStore),intent(inout) :: queue
      integer :: i
      DO I = 1, queue%allocated
         CALL mat_free(queue%D(i))
         CALL mat_free(queue%F(i))
         IF(queue%gd(i)%nrow > 0) CALL mat_free(queue%gd(i))
      END DO
      queue%allocated = 0
      queue%used_entries = 0
      queue%current_position = 0
#ifdef VAR_DSM
      av%dsm_pos = 0
#endif
      av%diis_pos = 0
      av%ediis_pos = 0
      deallocate(queue%D)
      deallocate(queue%F)
      deallocate(queue%gd)
      deallocate(queue%energy)
      deallocate(queue%metric)
      deallocate(queue%DIIS_Amat)

   END SUBROUTINE queue_free
   
   !> \brief Dump queue to disk.
   !> \author S. Host
   !> \date 2005
   SUBROUTINE queue_on_disk(queue,queue_lu,ndim)
      implicit none
      !> Subspace of Fock/KS and density matrices from previous SCF iterations
      TYPE(util_HistoryStore),intent(inout) :: queue
      !> Logical unit number of file to which the queue is written
      integer, intent(out) :: queue_lu
      !> Dimension of Fock/KS and density matrices
      integer, intent(out) :: ndim
      integer :: i, j
      logical :: OnMaster
      OnMaster=.TRUE.
      if (queue%used_entries > 0) then
         ndim = queue%F(1)%nrow
         queue_lu = -1
         CALL LSOPEN(queue_lu,'queue','unknown','UNFORMATTED')
         do i = 1, queue%used_entries
            call mat_write_to_disk(queue_lu,queue%D(i),OnMaster)
            CALL mat_free(queue%D(i))
         enddo
         do i = 1, queue%used_entries
            call mat_write_to_disk(queue_lu,queue%F(i),OnMaster)
            CALL mat_free(queue%F(i))
         enddo
      endif
   END SUBROUTINE queue_on_disk

   !> \brief Restore queue to disk.
   !> \author S. Host
   !> \date 2005
   SUBROUTINE queue_from_disk(queue,queue_lu,ndim)
      implicit none
      !> Subspace of Fock/KS and density matrices from previous SCF iterations
      TYPE(util_HistoryStore),intent(inout) :: queue
      !> Logical unit number of file on which the queue is stored
      integer, intent(inout) :: queue_lu
      !> Dimension of Fock/KS and density matrices
      integer, intent(in) :: ndim
      integer :: i, j
      logical :: OnMaster
      OnMaster = .TRUE.
      if (queue%used_entries > 0) then
         rewind(queue_lu)
         do i = 1, queue%used_entries
            CALL mat_init(queue%D(i), ndim, ndim)
            call mat_read_from_disk(queue_lu,queue%D(i),OnMaster)
         enddo
         do i = 1, queue%used_entries
            CALL mat_init(queue%F(i), ndim, ndim)
            call mat_read_from_disk(queue_lu,queue%F(i),OnMaster)
         enddo
      CALL LSCLOSE(queue_lu,'DELETE')
      endif
   END SUBROUTINE queue_from_disk

   !> \brief Add info to queue.
   !> \author L. Thogersen
   !> \date 2003
   SUBROUTINE add_to_queue(av, F, D, S, E, grad, queue)
      !use scf_stats
      implicit none
      !> Used to store info about SCF averaging
      type(avItem),intent(inout) :: av
      !> Fock/KS matrix
      TYPE(Matrix), INTENT(IN) :: F
      !> Density matrix
      TYPE(Matrix), INTENT(IN) :: D
      !> Overlap matrix
      TYPE(Matrix), INTENT(IN) :: S
      !> SCF energy
      real(realk), INTENT(IN)  :: E
      !> SCF gradient
      TYPE(Matrix), INTENT(IN) :: grad
      !> Subspace of Fock/KS and density matrices from previous SCF iterations
      TYPE(util_HistoryStore)  :: queue
      integer                  :: pos, i, j, posm1
      real(realk)              :: gradsqnorm, gdotprod

      !print *, "adding energy to queue", E
      !print *, 'queue%current_position', queue%current_position
      !print *, 'cfg_settings(cfg_set_type)%max_history_size', cfg_settings(cfg_set_type)%max_history_size
      !print *, 'MOD(queue%current_position, cfg_settings(cfg_set_type)%max_history_size)', &
      !         & MOD(queue%current_position, cfg_settings(cfg_set_type)%max_history_size)
      queue%current_position = &
           &MOD(queue%current_position, av%cfg_settings(av%cfg_set_type)%max_history_size)+1
      pos = queue%current_position
      IF(queue%allocated<queue%current_position) THEN
         CALL mat_init(queue%F(pos), F%nrow, F%ncol)
         CALL mat_init(queue%D(pos), D%nrow, D%ncol)
         if(av%save_gradient)&
              & CALL mat_init(queue%gd(pos),&
              &               F%nrow, F%ncol)
         queue%allocated = pos
      END IF
      if (queue%used_entries < queue%current_position) then
        queue%used_entries = pos
      endif
      call mat_assign(queue%F(pos),F)
      call mat_assign(queue%D(pos),D)
      queue%energy(pos) = E
      if (av%cfg_flush_vec.and.queue%used_entries == av%cfg_settings(av%cfg_set_type)%max_history_size) then
        if (pos /= 1) then
          posm1 = pos - 1
        else
          posm1 = av%cfg_settings(av%cfg_set_type)%max_history_size
        endif
        if (queue%energy(pos) > queue%energy(posm1)) then
          !Energy went up - step back and flush
          call flush_queue(av,S,pos,2,queue)
          av%flushed = .true.
          WRITE(av%LUPRI,*) 'WARNING - queue flushed and backstepped, it.:' !, stat_current_iteration+1
          RETURN
        else
          av%flushed = .false.
        endif
      endif
      IF(av%save_gradient) THEN
        call mat_assign(queue%gd(pos),grad)
        gradsqnorm = mat_sqnorm2(grad)
        !update diis pointer
        av%diis_pos = MOD(av%diis_pos,av%diis_history_size) + 1
        queue%DIIS_Amat(av%diis_pos,av%diis_pos) = 2.0E0_realk*gradsqnorm   
        i = av%diis_pos - 1
        if (i == 0) i = av%diis_history_size
        j = pos - 1
        if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
        do
          if (i == av%diis_pos .or. j > queue%used_entries) exit
          gdotprod = mat_dotproduct(queue%gd(j),queue%gd(pos))
          queue%DIIS_Amat(i,av%diis_pos) = 2.0E0_realk*gdotprod
          queue%DIIS_Amat(av%diis_pos,i) = queue%DIIS_Amat(i,av%diis_pos)
          i = i - 1
          if (i == 0) i = av%diis_history_size
          j = j - 1
          if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
        enddo
      ENDIF
      !update ediis pointer
      av%ediis_pos = MOD(av%ediis_pos,av%ediis_history_size) + 1
#ifdef VAR_DSM
      if(av%cfg_averaging == av%cfg_avg_dsm .or. av%cfg_lshift == av%cfg_lshift_dorth .or. av%DEBUG_RH_MU_E) then
        !update dsm pointer
        av%dsm_pos = MOD(av%dsm_pos,av%dsm_history_size) + 1
        call queue_update_TrDSDS(av,S,queue)
      endif
#endif

   END SUBROUTINE add_to_queue

   !> \brief Make TrDiSDjS array
   !> \author L. Thogersen
   !> \date 2003
   subroutine queue_update_TrDSDS(av,S,queue)
      implicit none
      !> Used to store info about SCF averaging
      type(avItem),intent(inout) :: av
      !> Overlap matrix
      type(MAtrix), INTENT(IN) :: S
      !> Subspace of Fock/KS and density matrices from previous SCF iterations
      TYPE(util_HistoryStore)  :: queue
      type(Matrix)             :: SDpos, SDS ! ,pointer
      integer :: i, j, pos

      pos = queue%current_position
      call mat_init(SDpos,queue%D(pos)%nrow,queue%D(pos)%ncol)
      call mat_init(SDS,queue%D(pos)%nrow,queue%D(pos)%ncol)
      call mat_mul(S,queue%D(pos),'n','n',1E0_realk,0E0_realk,SDpos)
      call mat_mul(SDpos,S,'n','n',1E0_realk,0E0_realk,SDS)
#ifdef VAR_DSM
      !WRITE(LUPRI,*) 'dsm_pos',dsm_pos
      queue%metric(av%dsm_pos,av%dsm_pos) = mat_dotproduct(queue%D(pos),SDS)
      !write (LUPRI,*) "metric(i,i):", queue%metric(dsm_pos,dsm_pos)
      i = av%dsm_pos - 1
      if (i == 0) i = av%dsm_history_size
      j = pos - 1
      if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
      do
        if (i == av%dsm_pos .or. j > queue%used_entries) exit
        queue%metric(i,av%dsm_pos) =  mat_dotproduct(queue%D(j),SDS)
        queue%metric(av%dsm_pos,i) = queue%metric(i,av%dsm_pos)
        i = i - 1
        if (i == 0) i = av%dsm_history_size
        j = j - 1
        if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
      enddo         
#else
     call lsquit('Queue_update_TrDSDS',-1)
#endif
      call mat_free(SDpos)
      call mat_free(SDS)
      !write(LUPRI,*) "Matrix array after update"
      !i = MIN(queue%used_entries,dsm_history_size)
      !CALL OUTPUT(queue%metric,1,i,1,i,max_history_size,max_history_size,&
      !      &      1,lupri)
   end subroutine Queue_update_TrDSDS

   !> \brief Flush the queue.
   !> \author L. Thogersen
   !> \date 2003
   !> 
   !> To flush the queue and set certain entry from queue and a certain
   !> number of older entries as the first and only entries
   !>
   subroutine flush_queue(av,S,keep_entry,nentries,queue)
     implicit none
     !> Used to store info about SCF averaging
     type(avItem),intent(inout) :: av
     !> Overlap matrix
     type(Matrix), intent(in) :: S
     !> The entry to keep
     integer, intent(in) :: keep_entry
     !> Number of older entries to keep (incl. keep_entry)
     integer, intent(in) :: nentries
     !> Subspace of Fock/KS and density matrices from previous SCF iterations
     type(util_HistoryStore),intent(inout) :: queue
     type(Matrix),allocatable :: save(:)
     real(realk), allocatable :: saveE(:)
     real(realk) :: gradsqnorm,gdotprod
     integer :: i,j,pos,k

! care is taken not to overwrite stuff we would like to keep!
! Reallocate stuff in queue
     if (nentries > 0) then
        ALLOCATE(save(nentries),saveE(nentries))
        pos = keep_entry
        do i = nentries,1,-1
          call mat_init(save(i),S%nrow,S%ncol)
          call mat_assign(save(i),queue%D(pos))
          saveE(i) = queue%energy(pos)
          pos = pos - 1
          if (pos == 0) pos = queue%used_entries
        enddo
        pos = keep_entry
        do i = nentries,1,-1
          call mat_assign(queue%D(i),save(i))
          queue%energy(i) = saveE(i)
          call mat_assign(save(i),queue%F(pos))
          pos = pos - 1
          if (pos == 0) pos = queue%used_entries
        enddo
        pos = keep_entry
        do i = nentries,1,-1
          call mat_assign(queue%F(i),save(i))
          if (av%save_gradient)then
             call mat_assign(save(i),queue%gd(pos))
          endif
          pos = pos - 1
          if (pos == 0) pos = queue%used_entries
        enddo
        do i = 1,nentries
          if (av%save_gradient) call mat_assign(queue%gd(i),save(i))
          call mat_free(save(i))
        enddo
        DEALLOCATE(save,saveE)
     endif
!
!  Handle density subspace minimization stuff
!     
#ifdef VAR_DSM
    av%dsm_pos = 0
#endif
     av%diis_pos = 0
     if (nentries > 0) then
       do k = 1,nentries
         queue%current_position = k
         queue%used_entries = k
#ifdef VAR_DSM
         av%dsm_pos = MOD(av%dsm_pos,av%dsm_history_size) + 1
#endif
         av%diis_pos = MOD(av%diis_pos,av%diis_history_size) + 1    
#ifdef VAR_DSM
         if(av%cfg_averaging == av%cfg_avg_dsm .or. av%DEBUG_RH_MU_E .or.&
            &av%cfg_lshift == av%cfg_lshift_dorth ) then
            !** Make TrDiSDjS array
            call queue_update_TrDSDS(av,S,queue)
         endif
#endif
         pos = k

         IF(av%save_gradient) THEN
            gradsqnorm = mat_sqnorm2(queue%gd(pos))
            queue%DIIS_Amat(av%diis_pos,av%diis_pos) = 2.0E0_realk*gradsqnorm   
            i = av%diis_pos - 1
            if (i == 0) i = av%diis_history_size
            j = pos - 1
            if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
            do
              if (i == av%diis_pos .or. j > queue%used_entries) exit
              gdotprod = mat_dotproduct(queue%gd(j),queue%gd(pos))
              queue%DIIS_Amat(i,av%diis_pos) = 2.0E0_realk*gdotprod
              queue%DIIS_Amat(av%diis_pos,i) = queue%DIIS_Amat(i,av%diis_pos)
              i = i - 1
              if (i == 0) i = av%diis_history_size
              j = j - 1
              if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
            enddo
         ENDIF
       enddo
       av%ediis_pos = MOD(nentries,av%ediis_history_size)
     else
       queue%used_entries = 0
       queue%current_position = 0
       av%ediis_pos = 0
     endif
   end subroutine flush_queue

   !> \brief Makes Dbar and Fbar, that is the weighted average of the density and fockmatrices
   !> \author L. Thogersen
   !> \date 2003
   SUBROUTINE get_average_arr(av,queue_to_average,queue, his_size, weights, averaged_arr)
      implicit none
      !> Used to store info about SCF averaging
      type(avItem), intent(inout)         :: av
      !> 'D' if averaging densitites, 'F' if averaging Fock/KS matrices
      character, intent(in)               :: queue_to_average
      !> Subspace of Fock/KS and density matrices from previous SCF iterations
      TYPE(util_HistoryStore),intent(in)  :: queue
      !> Number of Fock/density matrices in queue
      integer, intent(in)                 :: his_size
      !> The weights to use for weighted average
      real(realk), INTENT(IN)             :: weights(*)
      !> The average density (if queue_to_average='D') or average Fock/KS matrix (if queue_to_average='F')
      TYPE(Matrix)                        :: averaged_arr
      integer                         :: i, his_pos, pos, j

      his_pos = queue%current_position
      if (his_size == av%cfg_settings(av%cfg_set_type)%max_history_size) then
        pos = his_pos
#ifdef VAR_DSM
      elseif (his_size == av%dsm_history_size) then
        pos = av%dsm_pos
#endif
      elseif (his_size == av%diis_history_size) then
        pos = av%diis_pos
      else
        STOP 'unknown history_size in get_average_arr'
      endif
      if (pos == 0) then
        WRITE(av%LUPRI,*) 'WARNING[get_average_arr] : Tried to average from 0 matrices, exiting'
        RETURN
      endif
      if (queue_to_average == 'D') then
        call mat_copy(weights(pos),queue%D(his_pos),averaged_arr)
        i = pos-1
        if (i == 0) i = his_size
        j = his_pos-1
        if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
        do
          if (i == pos .or. j > queue%used_entries) exit
          call mat_daxpy(weights(i),queue%D(j),averaged_arr)
          i = i - 1
          if (i == 0) i = his_size
          j = j - 1
          if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
        enddo
      elseif (queue_to_average == 'F') then 
        call mat_copy(weights(pos),queue%F(his_pos),averaged_arr)
        i = pos-1
        if (i == 0) i = his_size
        j = his_pos-1
        if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
        do
          if (i == pos .or. j > queue%used_entries) exit
          call mat_daxpy(weights(i),queue%F(j),averaged_arr)
          i = i - 1
          if (i == 0) i = his_size
          j = j - 1
          if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
        enddo
      else
        WRITE(av%LUPRI,*) 'unknown queue to average from, got ',&
                    & queue_to_average,' but expected "D" or "F"'
        STOP 'unknown queue to average from, expected D or F'
      endif
   END SUBROUTINE get_average_arr
  
   !> \brief Get coefficients for part of density that can be expanded by subspace
   !> \author L. Thogersen
   !> \date 2003
   !>
   !> Din = D(orth) + D(para);   INPUT \n
   !> D(para) = sum_i(c_i D_i); \n
   !> c_i = sum_j({(S[2])^-1}_ij<D_j|Din>;   OUTPUT \n
   !>       <D_j|Din> = TrD_j S Din S; S[2]_ij = Tr D_i S D_j S \n
   !>
   subroutine util_GET_PROJ_PART(av,queue,Din,S,coef)
     implicit none
     !> Used to store info about SCF averaging
     type(avItem), intent(inout)         :: av
     !> Subspace of Fock/KS and density matrices from previous SCF iterations
     type(util_HistoryStore),intent(in)  :: queue
     !> Input density, Din = D(orth) + D(para)
     type(matrix), intent(in)            :: Din
     !> Overlap matrix
     type(matrix), intent(in)            :: S
     !> Coefficients c in D(para) = sum_i(c_i D_i)
     real(realk), intent(inout) :: coef(*)
     real(realk),allocatable :: inv_metric(:,:)
     real(realk) :: wrk(2*queue%used_entries*queue%used_entries)
     real(realk) :: vec(queue%used_entries)
     type(matrix) :: tmp, SDinS ! pointer
     real(realk) :: TrDjSDinS
     integer :: N,i1,j,k,i2
#ifdef VAR_DSM
 !** Initializations
!     N = queue%used_entries    !number of stored fock and density matrices
     N = MIN(queue%used_entries,av%dsm_history_size) !dimension of queue%metric
     coef(1:N) = 0.0E0_realk
     call mat_init(tmp,S%nrow,S%nrow)
     call mat_init(SDinS,S%nrow,S%nrow)

     !** S[2]^-1
     ALLOCATE(inv_metric(N,N))
     inv_metric = queue%metric(1:N,1:N)
     call INVERT_BY_DIAG2(inv_metric,wrk,wrk(N*N+1),vec,N)
     !** SDmuS
     call mat_mul(Din,S,'n','n',1E0_realk,0E0_realk,tmp)
     call mat_mul(S,tmp,'n','n',1E0_realk,0E0_realk,SDinS)
 
 !** c_i = sum_j({(S[2])^-1}_ij<D_j|Din>; <D_j|Din> = TrD_j S Din S
     j = queue%current_position
     i2 = av%dsm_pos
     do k = 1,N 
       TrDjSDinS = mat_dotproduct(queue%D(j),SDinS)
       do i1 = 1,N
         coef(i1) = coef(i1) + inv_metric(i1,i2)*TrDjSDinS
       enddo
       i2 = i2 - 1
       if (i2 == 0) i2 = av%dsm_history_size
       j = j - 1
       if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
     enddo
     if (av%info_d_proj) WRITE(av%LUPRI,*) 'coef in GET_PROJ_PART:',(coef(i1),i1=1,N)
     call mat_free(tmp)
     call mat_free(SDinS)
#else
     call lsquit('util_GET_PROJ_PART',-1)
#endif
   end subroutine util_GET_PROJ_PART
end MODULE av_utilities
