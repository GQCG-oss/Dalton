!> @file 
!> Contains density optimization module.

!> \brief Wrapper module for different types of density optimization.
!> \author L. Thogersen. Documented by S. Host.
!> \date 2003
MODULE density_optimization
   use arhModule
   use configurationType
   use files
   use queue_ops
   use diagonalization
   use matrix_util
   PUBLIC :: DOPT_get_density, get_density
   PRIVATE   

CONTAINS
   !> \brief Does different initializations before get_density is called.
   !> \author L. Thogersen
   !> \date 2003
   !>
   !>  For the current Fock matrix, the density is found, 
   !>  either by diagonalizing the Fock matrix or some alternative scheme.
   !>  Different kinds of dynamical levelshifting are available.
   !>  NOTE - This means that if F(i) is used, D(i+1) is found. E = Tr(h + F(i))D(i)
   !>
   SUBROUTINE DOPT_get_density(config,fifoqueue, queue, F, H1, D, SCF_it,ls)
      IMPLICIT NONE
      type(lsitem) ::ls
      !> Contains all info about configuration/settings for SCF calculation
      type(configItem),intent(inout) :: config
      !> New queue type: Contains Fock/KS and density matrices from previous SCF iterations (if ARH)
      type(modFIFO) :: fifoqueue
      !> Old queue type: Contains Fock/KS and density matrices from previous SCF iterations (if DIIS or DSM)
      TYPE(util_historyStore),target    :: queue
      !> Fock/KS matrix. Output only if cfg_scale_virt=.true.
      TYPE(Matrix), intent(inout):: F
      !> One-electron Hamiltonian
      TYPE(Matrix), INTENT(IN)   :: H1
      !> New density matrix (output)
      TYPE(Matrix), intent(inout):: D
      !> Number of current SCF iteration
      integer, intent(in) :: SCF_it
      real(realk) :: mu,ratio
      type(Matrix) :: Dnew
      integer :: ndia, dens_lun,idum,ldum, ndim
!** test nov2005
   real(realk) :: orbE,gradnorm
   type(Matrix) :: grad
   logical :: OnMaster
   OnMaster=.TRUE.
!
!** Initializations
      ndim = config%decomp%S%nrow
      ndia = 0
      SCF_iteration = SCF_it
      if (config%opt%CFG_density_method /= config%opt%CFG_F2D_arh) call mat_init(Dnew, ndim,ndim)

      found_density = .false.
      SDSexists = .false.
      mu = 0.0E0_realk 
      if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
         Drecent => fifoqueue%D_exp
      else
         Drecent => queue%D(queue%current_position) 
      endif

      if (config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens) then
         call mat_assign(Dnew,queue%D(queue%current_position))
      endif

      !** If virtual part of F should be scaled
      if (config%opt%cfg_scale_virt)then
         call lsquit('diag_scale_virt_fock removed',-1)
         ! call diag_scale_virt_fock(config%decomp%S,H1,F,D,config%decomp%nocc)
      endif
!
!** Level-shift ??
!
      IF(config%opt%cfg_density_method == config%opt%cfg_f2d_roothaan .and. &
       & config%diag%cfg_lshift /= diag_lshift_none) THEN
         !
         ! Level shifted mode - always shift F by some scalar * SDS
         !
         call mat_init(SDS,ndim,ndim)
         !** Find SDS once and for all :-)
         call mat_mul(D,config%decomp%S,'n','n',1E0_realk,0E0_realk,Dnew)
         call mat_mul(config%decomp%S,Dnew,'n','n',1E0_realk,0E0_realk,SDS)
         SDSexists = .true.
         ndia = 0   !number of fock diagonalizations
         !** Find level-shift
         call level_shift(config%av,config%diag,config%opt,F,H1,config%decomp%S,D,queue,ndia,mu,Dnew)
      ENDIF
       
!
!** Get the density
!
      if (.not. found_density .and. config%diag%cfg_no_confs_checked_in_rh) then 
        !If the density was not found during level shift search and if no configuration shift
        !was testet, then find the final density.  
        if (config%opt%CFG_density_method == config%opt%CFG_F2D_arh) then !For arh/DD it is not necessary to save the previous density
           call get_density(config,F,mu,D,fifoqueue,H1,ls)
        else
           call get_density(config,F,mu,Dnew,fifoqueue,H1,ls)
        endif
        if (config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
          & config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then 
           mu = -config%solver%current_mu
        else
           ndia = ndia + 1  !No diagonalizations if direct density optimization
        endif
      endif
       !** all kinds of print, info og debug stuff
      call print_and_stat(config%av,config%diag,F,D,config%decomp%S,H1,Dnew,queue,mu,ndia)
      !** Updating the density - observe the order!
      if (config%opt%CFG_density_method /= config%opt%CFG_F2D_arh) then
         call mat_assign(D, Dnew)  !Why the call to mat_assign here?? why not just D = Dnew ??
         call mat_free(Dnew)
      endif
!
!     Symmetrize the newly found density matrix
!     The Dmat should be symmetric but the norm of antisymmetric part 
!     can be of size 10.0E-12_realk which is enough that the 
!     subroutine mat_get_isym will denote it as nonsymmetric and then
!     the integral code will split the matrix up into a symmetric and 
!     antisymmetric contribution which is not good for performance. 
!     Alternatively the Integral code will construct an almost 
!     Symmetric Fock matrix from the almost symmetric Density which  
!     can give problems for the LAPACK, when the input is not fully symmetric. 
!     TK 

      call util_get_symm_part(D)

!
!** Write density to disk for possible restart
!      
      IF(config%decomp%cfg_DumpDensRestart)THEN
         dens_lun = -1
         IF(config%opt%optlevel.EQ.2)THEN
            call lsopen(dens_lun,'vdens.restart','UNKNOWN','UNFORMATTED')
         ELSE
            call lsopen(dens_lun,'dens.restart','UNKNOWN','UNFORMATTED')
         ENDIF
         rewind dens_lun
         call mat_write_to_disk(dens_lun,D,OnMaster)
         call mat_write_info_to_disk(dens_lun,config%decomp%cfg_gcbasis)
         call lsclose(dens_lun,'KEEP')
      ENDIF
!
!** Finalize      
      if (config%diag%cfg_lshift .ne. diag_lshift_none) then
        call mat_free(SDS)
        SDSexists = .false.
      endif
      NULLIFY(Drecent)
   END SUBROUTINE DOPT_get_density

   !> \brief This routine branches out and calls the requested routine for evaluating the density.
   !> \author L. Thogersen
   !> \date 2003
   !>
   !> If diagonalization, the level shift mu is determined before
   !> diagonalization. If ARH/TrFD, mu is determined dynamically during the
   !> solution of the linear equations. \n  
   !> Currently (April 2010) the options are: 
   !> - Diagonalization 
   !> - Direct Density Optimization (TrFD)
   !> - Purification
   !> - Augmented Roothaan-Hall (ARH)
   !>
   subroutine get_density(config,F,mu,Dnew,fifoqueue,H1,ls)
     implicit none
     type(lsitem) ::ls
     !> Contains all info about configuration/settings for SCF calculation
     type(configItem), intent(inout)     :: config
     !> Fock/KS matrix
     type(Matrix), intent(in) :: F
     !> Level shift
     real(realk), intent(in) :: mu
     !> New density matrix (output)
     type(Matrix), intent(inout) :: Dnew
     !> One-electron Hamiltonian
     TYPE(Matrix), INTENT(IN)   :: H1
     !> Contains Fock/KS and density matrices from previous SCF iterations
     type(modFIFO), optional :: fifoqueue
     type(Matrix) :: Fdamp
     logical        :: shiftedF
     integer        :: ndim
!test
     real(realk) :: orbE,gradnorm
     type(Matrix) :: grad

     ndim = config%decomp%S%nrow
     shiftedF = .false.
     if (config%diag%cfg_lshift /= diag_lshift_none .and. &
     &   config%opt%cfg_density_method /= config%opt%cfg_f2d_direct_dens .and. &
     &   ABS(mu) > 1E-8_realk) then
       !Evaluate the shifted Fock matrix
       ! Fdamp = F - muSDS
        call mat_init(Fdamp,ndim,ndim)
        CALL mat_add(1E0_realk,F,-mu, SDS, Fdamp)
        shiftedF = .true.
     endif
     if (config%opt%cfg_density_method == config%opt%cfg_f2d_roothaan) then
       ! Diagonalization
       if(shiftedF) then
         call diag_f_to_get_d(config%diag,Fdamp,config%decomp%S,Dnew)
       else
         call diag_f_to_get_d(config%diag,F,config%decomp%S,Dnew)
       endif
     else if (config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens) then
       call arh_get_density(config%solver,config%decomp,F,fifoqueue,Dnew,SCF_iteration,config%davidSCF,H1,ls)
     !Purification no longer supported! /Stinne 16-08-2010
     !else if (config%opt%cfg_density_method == config%opt%cfg_f2d_purification) then
     !  ! Purification scheme
     !  if (shiftedF) then
     !    call purification_update(config%opt,config%decomp%nocc,Fdamp,config%decomp%S,Dnew)
     !  else
     !    call purification_update(config%opt,config%decomp%nocc,F,config%decomp%S,Dnew)
     !  endif
     else if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then 
         if (present(fifoqueue)) then
            call arh_get_density(config%solver,config%decomp,F,fifoqueue,Dnew,SCF_iteration,config%davidSCF,H1,ls)
         else
            STOP 'No queue in call to arh_get_density'
         endif
     else 
     !CASE DEFAULT
        WRITE(config%LUPRI,*) 'No appropriate type of density evaluation has been chosen'
        STOP 'No appropriate type of density evaluation has been chosen'
     endif
     !END SELECT 

     if (shiftedF) then
      call mat_free(Fdamp)
     endif

   end subroutine get_density

END MODULE density_optimization
