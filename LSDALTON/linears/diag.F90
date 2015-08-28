!> @file 
!> Contains diagonalization module.

!> Diagonalization (standard Roothaan-Hall) module. Contains also the code for level shifting used in combination with RH.
!> \author L. Thogersen. Documented by S. Host.
!> \date 2003
MODULE diagonalization
   use av_utilities
   use files
   USE matrix_operations
   USE matrix_operations_aux, only:mat_column_norm,mat_density_from_orbs
   use matrix_util
   use memory_handling
   !> Smallest accepted overlap when using levelshift = MOchange
   real(realk),save :: accpr 
   !> Largest accepted Dorth ratio when using levelshift = Dorth
   real(realk),save :: limit_ratio 
   !> HOMO-LUMO gap
   real(realk),save :: Ggap 
   !> Overlap between old and new MO coefficients
   real(realk),save :: Gmochange
   !> Coefficients for idempotent density
   type(Matrix),save :: Co
   !> Saves input density for debug purposes
   type(Matrix),save :: Dinp
   !> S*D*S, D = density and S = overlap matrix
   type(Matrix),save :: SDS
   !> Pointer to current density matrix
   type(Matrix),pointer,save :: Drecent
   !> Current SCF iterationDEBU
   integer, save :: SCF_iteration
   !> States if acceptable density been found (levelshift)
   logical,save :: found_density 
   !> States if the matrix SDS been initialized and constructed
   logical,save :: SDSexists
   !> Level shift found from min overlap of every single new occ MO on old occ MO space
!   integer,parameter :: Diag_lshift_MOchange = 1
   !> A search is made in Escf(mu) - expensive!
   integer,parameter :: Diag_lshift_search = 2
   !> Use ratio ||Dorth||/||D|| to find mu
   integer,parameter :: diag_lshift_dorth = 3
   !> Do no level shifting
   integer,parameter :: Diag_lshift_none = 4
   !> Van Lenthe levelshift (see J.Comput.Chem. 27, 926-932 (2005))
   integer,parameter :: Diag_lshift_vanlenthe = 5 

!> \author S. Host
!> \date March 2010
!> \brief Contains diagonalization parameters
type DiagItem
   !CONFIGURATION:
   !==============
   !> Logical unit number for LSDALTON.OUT
   integer :: lupri
   !> Specifies wether calculation should be restarted from old density
   logical :: CFG_restart
   !> Specifies wether calculation should be purify the old density
   logical :: CFG_purifyrestart
   !> Force level 2 re-calculation on calculations restarted from old density
   logical :: CFG_redo_l2 
   !> Specifies orbital-free DFT
   logical :: CFG_OrbFree 
   !> What kind of level shift should be used with diagonalization?
   integer :: CFG_lshift
   !> A minimum level shift can be set
   REAL(REALK) :: cfg_min_lshift
   !> Run with fixed level shift
   logical     :: cfg_fixed_shift
   !> Value of fixed level shift
   real(realk) :: cfg_fixed_shift_param
   !Run for a given number of iterations with userdefined (different) levelshifts in each iterations
   LOGICAL     :: cfg_custom_shift
   !> Do custom levelshift in this number of SCF iterations
   integer     :: cfg_nshifts 
   !> Should canonical MOs be stored on disk?
   logical     :: cfg_save_CMO
   !> In search of the level-shift sometimes a configuration shift is encountered 
   !> and tested and the chosen one put on the queue. Then this is set to .false.
   LOGICAL :: cfg_no_confs_checked_in_rh
   !> The configuration shift can cause problems for large calculations and is 
   !> therefore turned off. Can be turned on in input with .CONFSHIFT in *LINSCA section. 
   LOGICAL :: cfg_no_conf_shift
   !> True if unrestricted SCF
   logical :: cfg_unres
   !> Number of occupied orbitals (if restricted)
   integer     :: nocc
   !> Number of occupied alpha orbitals (if unrestricted)
   integer     :: nocca
   !> Number of occupied beta orbitals (if unrestricted)
   integer     :: noccb
   !> I *think* this means number of unpaired electrons - why is it called spin????
   logical    :: nofinalhomolumo
   !INFO:
   !=====
   logical :: info_rh_gap
   logical :: info_rh_gradient
   logical :: info_rh_iterations
   logical :: info_rh_mu
!   logical :: info_rh_mochange
   !DEBUGGING:
   !==========
   logical :: DEBUG_DCHANGE
   logical :: DEBUG_EMODEL_CHANGE
   logical :: debug_idempotency
   logical :: DEBUG_RH_MU_E
   !DATA:
   !=====
   REAL(REALK),pointer :: cfg_levelshifts(:)
   !> energy of Homo
   REAL(REALK) :: eHOMO
end type DiagItem

CONTAINS

!> \brief Set default configuration for diagonalization.
!> \author S. Host
!> \date March 2010
!> \param diag Used to store info about diagonalization
subroutine diag_set_default_config(diag)
implicit none
   type(DiagItem),intent(inout) :: diag

   diag%cfg_unres            = .false.
   diag%cfg_restart          = .false.
   diag%cfg_purifyrestart    = .false.
   diag%cfg_redo_l2          = .false.
   diag%cfg_OrbFree          = .false.

   diag%CFG_lshift = Diag_lshift_none !Default = no level shift
   diag%cfg_min_lshift        = 0.0E0_realk
   diag%cfg_fixed_shift       = .false.
   diag%cfg_fixed_shift_param = 0.0E0_realk
   diag%cfg_custom_shift      = .false.

   diag%cfg_save_CMO          = .false.
   diag%cfg_no_confs_checked_in_rh = .true.
   diag%cfg_no_conf_shift          = .true.
   !INFO:
   !=====
   diag%info_rh_gap           = .false.
   diag%info_rh_gradient      = .false.
   diag%info_rh_iterations    = .false.
   diag%info_rh_mu            = .false.
!   diag%info_rh_mochange      = .false.
   !DEBUG:
   !======
   diag%DEBUG_DCHANGE         = .false.
   diag%DEBUG_EMODEL_CHANGE   = .false.
   diag%debug_idempotency     = .false.
   diag%DEBUG_RH_MU_E         = .false.
   diag%nofinalhomolumo       = .false.
   diag%eHOMO                 = 0.0E0_realk

end subroutine diag_set_default_config

!> \brief Initialize the diagonalization module
!> \author L. Thogersen
!> \date 2003
!> \param Ndim Number of basis functions
!> \param param1 The smallest accepted projection of the new density on the incoming one
!> \param param2 The maximum accepted ratio in get_mu_dorth
!> \param diag Used to store info about diagonalization
   subroutine dopt_get_density_INIT(Ndim,param1,param2,diag)
      implicit none
      real(realk), intent(in) :: param1, param2
      integer, intent(in) :: Ndim
      type(diagItem),intent(in) :: diag

      if (diag%CFG_lshift .ne. diag_lshift_none) then
        !The smallest accepted projection of the new density on the incoming one
        !usually 0.975 for HF and 0.98 for DFT
        ! accpr = cfg_settings(CFG_set_type)%min_density_overlap
        accpr = param1
         WRITE(diag%LUPRI,*) 'smallest accepted overlap:',accpr
        !The maximum accepted ratio in get_mu_dorth
         !limit_ratio = cfg_settings(CFG_set_type)%max_dorth_ratio
         limit_ratio = param2
         WRITE(diag%LUPRI,*) 'maximum accepted ratio in dorth:',limit_ratio
      endif
!      if (diag%cfg_lshift == diag_lshift_MOchange .or. diag%debug_rh_mu_E) then
!        call mat_init(Co,ndim,ndim)
!      endif
      if (diag%DEBUG_DCHANGE) then
        call mat_init(Dinp,ndim,ndim)
      endif
      NULLIFY(Drecent)
   end subroutine dopt_get_density_INIT
    
!> \brief Shut down the diagonalization module
!> \author L. Thogersen
!> \date 2003
!> \param diag Used to store info about diagonalization
   subroutine dopt_get_density_SHUTDOWN(diag)
     implicit none
     type(diagItem),intent(in) :: diag

!     if (diag%cfg_lshift == diag_lshift_MOchange .or. diag%debug_rh_mu_E) then
!       call mat_free(Co)
!     endif
     if (diag%DEBUG_DCHANGE) then
       call mat_free(Dinp)
     endif
   end subroutine dopt_get_density_SHUTDOWN

!> \brief This routine branches out and calls the requested routine for obtaining the level shift.
!> \author L. Thogersen
!> \date October 2005
!> \param av Used to store info about SCF averaging
!> \param diag Used to store info about diagonalization
!> \param opt Used to store info about SCF optimization type
!> \param F Fock/Kohn-Sham matrix
!> \param H1 One-electron Hamiltonian
!> \param S Overlap matrix
!> \param D Density matrix
!> \param queue Contains density and Fock/KS matrices from previous SCF iterations
!> \param ndia Number of diagonalizations in this SCF iteration
!> \param mu Levelshift
!> \param Dnew New density matrix
!>
!> This routine branches out and calls the requested routine for obtaining the level shift.
!> Currently (October 2005) the options are: No level shift, level shift obtained by a
!> line-search in the SCF energy, level shift obtained based on restrictions on new
!> information introduced to the subspace of densities and level shift obtained based on
!> restrictions on new information introduced compared to the previous density.
!>
   subroutine level_shift(av,diag,opt,F,H1,S,D,queue,ndia,mu,Dnew)
   use opttype
      implicit none
      type(avItem),intent(inout)     :: av
      type(DiagItem),intent(inout)   :: diag
      type(optItem),intent(in)       :: opt
      TYPE(util_HistoryStore),intent(in) :: queue
      type(Matrix),intent(in) :: F,H1,S,D
      real(Realk), intent(out) :: mu
      type(Matrix), intent(inout) :: Dnew !output
      integer, intent(inout) :: ndia
      real(realk) :: mu_min

      write(diag%lupri,*) 'diag%cfg_lshift:', diag%cfg_lshift
      if (diag%cfg_lshift == diag_lshift_none) then
         !do no lshift, so return 0 for mu
         mu=0E0_realk 
      else if (diag%cfg_lshift == diag_lshift_search) then
         !** make a linesearch in E_SCF(mu) using the mu giving the
         !** largest decrease in ESCF
         call get_muopt(diag,F,S,mu)
         if (diag%info_rh_mu) then
            WRITE(diag%LUPRI,*) 'using linesearch mu = ',mu,'   grep'
         endif
      else if (diag%cfg_lshift == diag_lshift_dorth) then
         !Level shift obtained based on restrictions on new information
         !introduced to the subspace of densities
         call get_mu_dorth(av,diag,queue,F,S,ndia,mu,Dnew)      
         mu = MAX(mu,diag%cfg_min_lshift)
         !if (cfg_fixed_shift) then
         !   mu = cfg_fixed_shift_param
         !endif
         if (diag%info_rh_mu) WRITE(diag%LUPRI,*) 'using dorth or mu_mindamp',mu
!!$      else if (diag%cfg_lshift == diag_lshift_MOchange) then
!!$         !Level shift obtained based on restrictions on new information introduced compared to 
!!$         !the previous density.
!!$         !Only relevant for diagonalization since MOs are used
!!$         if (opt%cfg_density_method /= opt%cfg_f2d_roothaan) then
!!$           WRITE(diag%LUPRI,*) 'The level-shift mode .MOCHANGE can only be run with diagonalization of the ',&
!!$                        &'Fock matrix - .RH . If diagonalization is unwanted choose instead the level-shift mode ',&
!!$                        &'.DORTH'
!!$           STOP 'wrong combination of level-shift mode and density evaluation - see output'
!!$         endif
!!$         call get_mu_mochange(av,diag,queue,F,D,S,ndia,mu)
!!$         mu = MAX(mu,diag%cfg_min_lshift)
!!$         if (diag%info_rh_mu) WRITE(diag%LUPRI,*) 'using mu_mochange or mu_mindamp = ',&
!!$              & mu,'  grep'
      else if (diag%cfg_lshift == diag_lshift_vanlenthe) then
         !done below - diag%cfg_fixed_shift and diag%cfg_custom_shift are true.
      else
         WRITE(diag%LUPRI,*) 'No appropriate type of level-shift has been chosen'
         STOP 'No appropriate type of level-shift has been chosen'
      endif

      if (diag%cfg_fixed_shift) then
         if (diag%cfg_custom_shift) then
            mu = diag%cfg_levelshifts(SCF_iteration)
            !write(diag%lupri,*) 'SCF it:', SCF_iteration
            !write(diag%lupri,*) 'mu set to:', mu
         else
            mu = diag%cfg_fixed_shift_param
         endif
      endif

   end subroutine level_shift

!> \brief Solve linear eq. FC = SCe and create new density from C
!> \author L. Thogersen
!> \date October 2005
!> \param diag Used to store info about diagonalization
!> \param F Fock/Kohn-Sham matrix
!> \param S Overlap matrix
!> \param Dnew Density matrix
   subroutine diag_f_to_get_d(diag,F,S,Dnew)
     implicit none
     type(DiagItem),intent(inout) :: diag !add eHOMO otherwise intent(in)
     type(Matrix), intent(in) :: F, S
     type(Matrix), intent(inout) :: Dnew  !output
     type(Matrix)             :: Cmo
     real(realk)              :: minoverlap(1)
     real(realk), pointer :: eival(:)
     integer :: ndim, cmo_lun,idum,ldum,sz !dum = dummies
     real(realk),external :: HOMO_LUMO_gap
     real(realk),external :: HOMO_energy
     logical :: OnMaster
     ndim = S%nrow
     OnMaster=.TRUE.
     call mat_init(Cmo,ndim,ndim)
     if (diag%cfg_unres) then
        sz = 2*ndim
     else
        sz = ndim
     endif

     call mem_alloc(eival,sz) ! allow for unrestricted

     call mat_diag_f(F,S,eival,Cmo)

     call mat_density_from_orbs(Cmo,Dnew,diag%nocc,diag%nocca,diag%noccb,diag%cfg_OrbFree) 

!     if (diag%cfg_lshift == diag_lshift_MOchange .or. diag%DEBUG_RH_MU_E) then
!       call util_min_MO_overlap(diag%nocc,S,Co,Cmo,1,minoverlap(1))
!       Gmochange = minoverlap(1)
!     endif
     Ggap = HOMO_LUMO_gap(diag%cfg_unres,diag%nocc,diag%nocca,diag%noccb,eival,sz)
     diag%eHOMO = HOMO_energy(diag%cfg_unres,diag%nocc,diag%nocca,diag%noccb,eival,sz)

     if (diag%cfg_save_CMO) then
       cmo_lun = -1
       call lsopen(cmo_lun,'cmo.restart','UNKNOWN','UNFORMATTED')
       rewind cmo_lun
       call mat_write_to_disk(cmo_lun,Cmo,OnMaster)
       write(cmo_lun) Ggap
       call lsclose(cmo_lun,'KEEP')
     endif
 
     call mat_free(Cmo)
     call mem_dealloc(eival)

   end subroutine diag_f_to_get_d

!!> \brief Makes trace-purification update of the Fock matrix to obtain new D
!!> \author P. Salek
!!> \date 2003
!!> \param opt Used to store info about SCF optimization type
!!> \param nocc Number of occupied orbitals
!!> \param F Fock/Kohn-Sham matrix
!!> \param S Overlap matrix
!!> \param Dnew Density matrix
!   subroutine purification_update(opt,nocc,F,S,Dnew)
!   use opttype
!     implicit none
!     type(optItem),intent(in) :: opt
!     integer,intent(in)       :: nocc
!     type(Matrix), intent(in) :: F,S
!     type(matrix), intent(inout) :: Dnew
!     integer cycles
!
!     call purify(opt%lupri,F,S,nocc,opt%cfg_purification_method,Dnew)
!
!   end subroutine purification_update

!> \brief Level shift obtained based on restrictions on new information introduced compared to the previous density.  
!> \author L. Thogersen
!> \date October 2005
!> \param av Used to store info about SCF averaging
!> \param diag Used to store info about diagonalization
!> \param queue Contains density and Fock/KS matrices from previous SCF iterations
!> \param F Fock/Kohn-Sham matrix
!> \param D Density matrix
!> \param S Overlap matrix
!> \param ndia Number of diagonalizations in this SCF iteration
!> \param mu Levelshift
!>
!> Only relevant for diagonalization since MOs are used. Theory described in J.Chem.Phys. 123, 074103 (2005) sec. III.E.1
!> 
!!$   subroutine get_mu_mochange(av,diag,queue,F,D,S,ndia,mu)
!!$      implicit none
!!$      type(avItem), intent(inout) :: av
!!$      type(diagItem), intent(inout) :: diag
!!$      type(util_HistoryStore) :: queue
!!$      type(Matrix), intent(in) :: F,D,S
!!$      integer, intent(inout) :: ndia
!!$      real(realk), intent(out) :: mu
!!$      type(Matrix) :: Dnew
!!$      integer :: k
!!$      real(realk) :: mu_l, mu_r,mochange_0,mochange_r,mochange
!!$
!!$      call mat_init(Dnew,S%nrow,S%ncol)
!!$      !** diagonalize SDS (SDS Co = S Co e) giving the corresponding
!!$      ! idempotent Dnew = <Co|Co>  and the initial coefficients Co
!!$      call GET_Didem(diag,D,S,Dnew,Co)
!!$      k = 1   !number of diagonalizations
!!$      diag%cfg_no_confs_checked_in_rh = .true.
!!$!
!!$!** find brackets for binary search - mu = 0 is left bracket
!!$!
!!$      mu_l = 0.0E0_realk
!!$      call diag_Get_density(diag,F,S,mu_l,Dnew)
!!$      mochange_0 = Gmochange
!!$      k = k+1     
!!$      if (mochange_0 >= accpr) then
!!$        mu = mu_l
!!$        if (diag%info_rh_mochange) write(diag%lupri,*) 'mu_l, mochange',mu_l,mochange_0
!!$        call mat_free(Dnew)
!!$        RETURN
!!$      endif
!!$      !right bracket starts at mu = 10, then mu = mu+5 aso.
!!$      mu_r = mu_l + 10E0_realk
!!$      do
!!$         !find new density for (F-muSDS)C = SCe
!!$         call diag_Get_density(diag,F,S,mu_r,Dnew)
!!$         mochange = Gmochange
!!$         k = k + 1
!!$         if (mochange >= accpr) then
!!$            !right bracket found
!!$            if(diag%info_rh_mochange) write(diag%lupri,*) 'mu_l,mu_r,mochange',mu_l,mu_r,mochange
!!$            mochange_r = mochange
!!$            exit
!!$         else
!!$            !update left and right bracket
!!$            mu_l = mu_r
!!$            mu_r = mu_r + 5.0E0_realk
!!$         endif
!!$         if (mu_r > 100.0E0_realk) then
!!$           if (diag%info_rh_mochange) WRITE(diag%LUPRI,*) 'mu_r = ',mu_r ,' in mochange_damp, but mochange is ',mochange
!!$           STOP 'something wrong in find_mochange_damp1'
!!$         endif
!!$      enddo
!!$!
!!$!** Binary search
!!$!
!!$      do
!!$         mu = (mu_l + mu_r)/2.0E0_realk
!!$         if ((mu_r - mu_l) < 0.2E0_realk) then
!!$            if (diag%info_rh_mochange) WRITE(diag%LUPRI,*) 'mochange_0,mochange_r',mochange_0,mochange_r
!!$            !** Check for configuration shifts:
!!$            if (mochange_0 < 0.5E0_realk .and. mochange_r > accpr+(1-accpr)/2.0E0_realk) then  
!!$              diag%cfg_no_confs_checked_in_rh = .false.
!!$              call check_conf(av, diag,F,S,mu_r,queue,Dnew,mu,mochange)
!!$              k = k + 2
!!$            else
!!$              mu = mu_r
!!$              mochange = mochange_r
!!$              if (diag%info_rh_mu) then
!!$                 WRITE(diag%LUPRI,"('mu_l,mu_r,mu,mochange',4f10.6,'  ',i3,' Ds were made in order to find this mu')") &
!!$                             & mu_l,mu_r,mu,mochange,k
!!$              endif
!!$            endif
!!$            !** Optimal mu found 
!!$            exit
!!$         endif
!!$         !** Find new density for (F-muSDS)C = SCe
!!$         call diag_Get_density(diag,F,S,mu,Dnew)
!!$         mochange = Gmochange
!!$         if (diag%info_rh_mochange) WRITE(diag%LUPRI,*) 'mu,mochange',mu,mochange
!!$         k = k + 1
!!$         if (mochange < accpr) then
!!$            mu_l = mu
!!$         else
!!$            mu_r = mu
!!$            mochange_r = mochange
!!$         endif
!!$      enddo
!!$      ndia = ndia + k
!!$      
!!$      call mat_free(Dnew)
!!$   end subroutine get_mu_mochange

!> \brief Level shift obtained based on restrictions on new information introduced to the subspace of densities
!> \author L. Thogersen
!> \date October 2005
!> \param av Used to store info about SCF averaging
!> \param diag Used to store info about diagonalization
!> \param queue Contains density and Fock/KS matrices from previous SCF iterations
!> \param F Fock/Kohn-Sham matrix
!> \param S Overlap matrix
!> \param ndia Number of diagonalizations in this SCF iteration
!> \param mu Levelshift
!> \param Dnew New density matrix
!>
!> Theory described in PhD thesis by L. Thogersen 2005
!> 
   subroutine get_mu_dorth(av,diag,queue,F,S,ndia,mu,Dnew)
      implicit none
      type(avItem), intent(inout) :: av
      type(diagItem), intent(inout) :: diag
      type(util_HistoryStore) :: queue
      type(Matrix), intent(in) :: F,S
      integer, intent(inout) :: ndia
      real(realk), intent(out) :: mu
      type(Matrix),intent(inout) :: Dnew  !output
      type(Matrix) :: D_r
      integer :: k
      real(realk) :: mu_l, mu_r,ratio, ratio_0,ratio_r,dummy

      diag%cfg_no_confs_checked_in_rh = .true.
!
!** Find brackets for binary search - mu = 0 is left bracket
!
      mu_l = 0.0E0_realk
      call diag_Get_density(diag,F,S,mu_l,Dnew)
      k = 1   !number of diagonalizations or corresponing
      call ratio_Dorth_D(av,queue,S,Dnew,ratio)
      ratio_0 = ratio
      if (diag%info_rh_mu) WRITE(diag%LUPRI,*) 'mu, ratio',mu_l,ratio
      if (ratio < limit_ratio) then
        mu = mu_l
        RETURN
      endif
      call mat_init(D_r,S%nrow,S%ncol)
      ! Right bracket starts at mu = 5, then mu = mu+5 aso.
      mu_r = mu_l + 5E0_realk
      do
         !find new density for (F-muSDS)C = SCe
         call diag_Get_density(diag,F,S,mu_r,D_r)
         k = k + 1
         call ratio_Dorth_D(av,queue,S,D_r,ratio) 
         if (diag%info_rh_mu) WRITE(diag%LUPRI,*) 'mu, ratio',mu_r,ratio
         if (ratio < limit_ratio) then
            !right bracket found
            ratio_r = ratio
            exit
         else
            !update left and right bracket
            mu_l = mu_r
            mu_r = mu_r + 5.0E0_realk
         endif
         if (mu_r > 100.0E0_realk) then
           STOP 'something wrong in find_mochange_damp2'
         endif
      enddo
!
!** Binary search
!
      do
         mu = (mu_l + mu_r)/2.0E0_realk
         !find new density for (F-muSDS)C = SCe
         call diag_get_density(diag,F,S,mu,Dnew)
         k = k + 1
         call ratio_Dorth_D(av,queue,S,Dnew,ratio)
         if (diag%info_rh_mu) WRITE(diag%LUPRI,*) 'mu, ratio',mu,ratio
         if (ratio > limit_ratio) then
            mu_l = mu
         else
            mu_r = mu
            ratio_r = ratio
            call mat_assign(D_r,Dnew)
         endif
         if (ratio_r > limit_ratio*0.9E0_realk .or. (mu_r - mu_l) < 0.3E0_realk) then 
!         if (ratio_r > limit_ratio*0.9E0_realk .or. (mu_r - mu_l) < 0.05E0_realk) then !Stinne - didn't use this
                                                                            !after all. cfg_no_conf_shift
                                                                            !instead 
            !** Check for configuration shift
            if (diag%cfg_no_conf_shift) then
               mu = mu_r
               call mat_assign(Dnew,D_r)
               found_density = .true.
            else
               if (ratio_r < limit_ratio*0.1E0_realk) then
               !if (ratio_r < limit_ratio*0.05E0_realk) then !Stinne - didn't use this after all. cfg_no_conf_shift
                                                       !instead
                  diag%cfg_no_confs_checked_in_rh = .false.
                  call check_conf(av, diag,F,S,mu_r,queue,Dnew,mu,dummy)
                  k = k + 2
               else
                  mu = mu_r
                  call mat_assign(Dnew,D_r)
                  found_density = .true.
               endif
            endif
            !** Optimal mu found
            exit
         endif
      enddo
      ndia = ndia + k
      call mat_free(D_r)
      
   end subroutine get_mu_dorth

!> \brief Calculate ||Dorth||/||D|| where ||Dorth|| is the part of the new density Dmu that 
!> cannot be expressed in the subspace of the previous densities
!> \author L. Thogersen
!> \date October 2005
!> \param av Used to store info about SCF averaging
!> \param queue Contains density and Fock/KS matrices from previous SCF iterations
!> \param S Overlap matrix
!> \param Dmu New density matrix
!> \param ratio ||Dorth||/||D||
   subroutine ratio_Dorth_D(av,queue,S,Dmu,ratio)
      implicit none
      type(avItem), intent(inout) :: av
      TYPE(util_HistoryStore),intent(in) :: queue
      type(Matrix), intent(in) :: S,Dmu
      real(realk), intent(inout) :: ratio

      !** Initializations
      call lsquit('ratio_Dorth_D',-1)

   end subroutine ratio_Dorth_D

!> \brief Testing the configuration at mu = 0 and mu = mu_r and returning
!> the mu with the lowest corresponding E_SCF.
!> \author L. Thogersen
!> \date October 2005
!> \param av Used to store info about SCF averaging
!> \param diag Used to store info about diagonalization
!> \param F Fock/Kohn-Sham matrix
!> \param S Overlap matrix
!> \param mu_r The mu That keeps the configuration
!> \param queue Contains density and Fock/KS matrices from previous SCF iterations
!> \param Dnew Density corresponding to optimal mu on exit
!> \param mu Optimal mu
!> \param mochange HOMO-LUMO gap for accepted density
!>
!> Theory described in J.Chem.Phys. 123, 074103 (2005), sec. V
!>
   subroutine check_conf(av,diag,F,S,mu_r,queue,Dnew,mu,mochange)
      use dal_interface
      use fock_evaluator

      type(avItem),intent(inout) :: av
      type(diagItem), intent(inout) :: diag
      type(util_HistoryStore) :: queue
      type(Matrix), intent(in) :: F,S
      real(realk), intent(in) :: mu_r  !the mu that keeps the configuration
      type(Matrix), intent(inout) :: Dnew  !D corresponding to opt. mu on exit
      real(realk), intent(out) :: mochange,mu
      type(Matrix) :: Fnew,Dnew2,Fnew2
      real(realk) :: gap,mochange_0,mochange_r,&
                   & ESCF_0, ESCF_mu 

      call mat_init(Dnew2,S%nrow,S%nrow)
      call mat_init(Fnew,S%nrow,S%nrow)
      call mat_init(Fnew2,S%nrow,S%nrow)
      !Abrubt configuration shift possible:
      !FIXME: maybe the next step should only be carried out if at least 5 its
      !has been made since the last time.
      !For mu = 0:
      call diag_Get_density(diag,F,S,0.0E0_realk,Dnew)
      mochange_0 = Gmochange
      call fck_get_fock(Dnew,Fnew,ESCF_0)
      !For mu = mu_r:
      call diag_Get_density(diag,F,S,mu_r,Dnew2)
      mochange_r = Gmochange
      call fck_get_fock(Dnew2,Fnew2,ESCF_mu)
      !WRITE(diag%LUPRI,*) 'THIS IS A CONFIGURATION SHIFT:',mochange_0,mochange_r
      WRITE(diag%LUPRI,*) 'WHICH ONE IS BEST? mu = 0 or mu = ',mu_r
      WRITE(diag%LUPRI,*) 'mu_r,E(0),E(mu_r)',mu_r,ESCF_0,ESCF_mu
      !compare the energies, the lowest one "wins"
      if (ESCF_mu < ESCF_0) then
         !Configuration shift not accepted
         mu = mu_r
         if (av%save_gradient) then
           call get_AO_gradient(Fnew2, Dnew2, S, Fnew)
         endif
         CALL add_to_queue(av, Fnew2, Dnew2, S, ESCF_mu, Fnew, queue) 
         mochange = mochange_r
      else
         !Configuration shift accepted
         mu = 0.0E0_realk
         !remove old vectors connected to the old configuration
         call flush_queue(av,S,queue%current_position,0,queue)
         WRITE(diag%LUPRI,*) '*** Removing all the old vectors connected to the old configuration'
         if (av%save_gradient) then
           call get_AO_gradient(Fnew,Dnew,S,Fnew2)
         endif
         call add_to_queue(av, Fnew,Dnew,S,ESCF_0, Fnew2,queue) 
         mochange = mochange_0
      endif
      call mat_free(Fnew)
      call mat_free(Fnew2)
      call mat_free(Dnew2)
   end subroutine check_conf

!> \brief Make a line search for the EHF minimum. Searches from mu found in the latest iteration and chose the first minimum.
!> \author L. Thogersen
!> \date October 2005
!> \param diag Used to store info about diagonalization
!> \param F Fock/Kohn-Sham matrix
!> \param S Overlap matrix
!> \param muopt Optimal mu giving the energy minimum
   subroutine get_muopt(diag,F,S,muopt)
      use dal_interface
      use fock_evaluator
      use scf_stats
      !make a line search for the EHF minimum
      !searches from mu found in the latest iteration and chose the first minimum
      implicit none
      type(DiagItem),intent(inout) :: diag
      type(Matrix), intent(in) :: F,S
      real(realk), intent(out) :: muopt
      type(Matrix) :: Dnew,Fnew
      real(realk) :: mu,mu0,dmu,EHFold,EHFnew,mus(50,2)
      integer :: nfock,n,i

!
!** Initializations
!
      n = 0
      nfock = 0
      mu0 = MAX(stat_tab(stat_current_iteration,7) ,0.0E0_realk) !mu in last iteration
      mus = 0.0E0_realk
      call mat_init(Dnew,S%nrow,S%ncol)
      call mat_init(Fnew,S%nrow,S%ncol)      
      call diag_get_density(diag,F,S,mu0,Dnew)
      call fck_get_fock(Dnew,Fnew,EHFold)
      nfock = nfock + 1
      if (diag%info_rh_mu) mus(nfock,1) = mu0; mus(nfock,2) = EHFold
      dmu = MAX(mu0/10.0E0_realk,0.1E0_realk)
!  
!** Establish search direction
!
      mu = mu0 + dmu
      call diag_get_density(diag,F,S,mu,Dnew)
      call fck_get_fock(Dnew,Fnew,EHFnew)
      nfock = nfock + 1
      if (diag%info_rh_mu) then
        mus(nfock,1) = mu; mus(nfock,2) = EHFnew
        write(diag%lupri,*) mus(nfock,1),mus(nfock,2)
      endif
      if (EHFnew >= EHFold) then
      !** mu is to the right of the minimum
          !mu0 is either to the right of the minimum or close to it
          if (mu0 < 1E-15_realk) then
            muopt = 0.0E0_realk
            call mat_free(Dnew)
            call mat_free(Fnew)
            RETURN
          endif           
          !Find the closest minimum from mu0 towards mu = 0
          do
            mu = mu0 - dmu
            if (mu < 0.0E0_realk) mu = 0.0E0_realk
            call diag_get_density(diag,F,S,mu,Dnew)
            call fck_get_fock(Dnew,Fnew,EHFnew)
            nfock = nfock + 1
            if (diag%info_rh_mu) then
              mus(nfock,1) = mu; mus(nfock,2) = EHFnew
              write(diag%lupri,*) mus(nfock,1),mus(nfock,2)
            endif
            if (EHFnew > EHFold) then
              !passed minimum
              muopt = mu0
              exit
            elseif (mu < 1E-7_realk) then
              !reached zero
              muopt = 0.0E0_realk
              exit
            endif
            EHFold = EHFnew
            mu0 = mu
            dmu = MAX(mu0/6.0E0_realk,0.1E0_realk)
          enddo 
      else
      !** mu is to the left of the minimum
          mu0 = mu
          EHFold = EHFnew
          do
            mu = mu0 + dmu
            call diag_get_density(diag,F,S,mu,Dnew)
            call fck_get_fock(Dnew,Fnew,EHFnew)
            nfock = nfock + 1
            if (diag%info_rh_mu) then
              mus(nfock,1) = mu; mus(nfock,2) = EHFnew
              write(diag%lupri,*) mus(nfock,1),mus(nfock,2)
            endif
            if (EHFnew > EHFold) then
              !mu is to the right of the minimum
              !mu0 is either to the left or close to the minimum
              muopt = mu0
              exit
            elseif (mu > 500.0E0_realk) then
              WRITE(diag%LUPRI,*) 'WARNING: error in get_muopt: no minimum found before mu = 500'
              muopt = mu0
              exit
            endif
            EHFold = EHFnew
            mu0 = mu
            dmu = MAX(mu0/6.0E0_realk,0.1E0_realk) 
          enddo        
      endif 
   
      if (diag%info_rh_mu) then
        WRITE(diag%LUPRI,*) 'constructed ',nfock,' Fock-matrices in the muopt process for mus:'
        do i = 1,nfock
          WRITE(diag%LUPRI,*) mus(i,1),mus(i,2)
        enddo
      endif
      call mat_free(Dnew)
      call mat_free(Fnew)

   end subroutine get_muopt

!> \brief All kinds of print, info and debug stuff
!> \author L. Thogersen
!> \date October 2005
!> \param av Used to store info about SCF averaging
!> \param diag Used to store info about diagonalization
!> \param F Fock/Kohn-Sham matrix
!> \param D Density matrix
!> \param S Overlap matrix
!> \param H1 One-electron Hamiltonian
!> \param Dnew The new found density
!> \param queue Contains density and Fock/KS matrices from previous SCF iterations
!> \param mu Level shift
!> \param ndia Number of diagonalizations in this SCF iteration
   subroutine print_and_stat(av,diag,F,D,S,H1,Dnew,queue,mu,ndia)
      !
      use scf_stats
      !Dnew is the new found density
      implicit none
      type(avItem),intent(inout) :: av
      type(diagItem),intent(inout) :: diag
      type(Matrix), intent(in) :: F,D,S,H1,Dnew
      TYPE(util_HistoryStore),intent(in) :: queue
      real(realk), intent(in) :: mu
      integer, intent(in) :: ndia
      type(matrix) :: grad,Dchange,Didem
      real(realk) :: EHFstart,orbE,TrDSDS,muopt,ratio,Eorb_idem,Epred_0
      integer :: ndim,i

      ndim = S%nrow
      if (diag%DEBUG_EMODEL_CHANGE .or. diag%debug_rh_mu_E) then
        !** diagonalize SDS (SDS C = SCe) giving the corresponding
        ! idempotent Didem = <C|C> 
        call mat_init(Didem,ndim,ndim)
        call GET_Didem(diag,D,S,Didem,Co)
      endif
      if (diag%debug_rh_mu_E) then
        !find the pseudo SCF-energy for the idempotent Dbar
         EHFstart = pEHF(Didem)
         WRITE(diag%LUPRI,*) 'SCF energy after averaging:',EHFstart
      endif
      !if (DEBUG_RH_MU_E) then
      !  if (queue%used_entries > 0) then
      !    Epred_0 = roothan_Epred(queue,H1,S,Didem)
      !  else
      !    Epred_0 = 0.0E0_realk
      !  endif
      !  Eorb_idem = 2E0_realk*mat_dotproduct(F,Didem)
      !endif
      !if (cfg_density_method == CFG_F2D_DIRECT_DENS) then
      !   Epred_0 = debug_dd_epred
      !endif

!
! Filling up the stat_tab table
!
       !** col7 for level-shift info
       stat_tab(stat_current_iteration+1,7) = mu

       !** col8 for additional RH info
!       if (diag%cfg_lshift == diag_lshift_MOchange) then
!         stat_tab(stat_current_iteration+1,8) = Gmochange
!       elseif (diag%cfg_lshift == diag_lshift_dorth) then
       if (diag%cfg_lshift == diag_lshift_dorth) then
         call ratio_Dorth_D(av,queue,S,Dnew,ratio)
         stat_tab(stat_current_iteration+1,8) = ratio
       else
         stat_tab(stat_current_iteration+1,8) = 0.0E0_realk
       endif
 
!
! Stuff written directly to output
!
      if (diag%debug_dchange) then
         !save the input D to later use in debug_dchange_module
         call mat_assign(Dinp,D)
      endif

      !** Diagonalization statistics
      if (diag%info_rh_iterations) then
        WRITE(diag%LUPRI,*) '** ',ndia,' diagonalizations were carried out **'
      endif

      if (diag%debug_idempotency) then
        TrDSDS = util_Snorm(D,S)
        WRITE(diag%LUPRI,*) 'TrDbarSDbarS = ',TrDSDS
      endif

      if (diag%INFO_RH_GAP) then   
        WRITE(diag%LUPRI,*) 'Final HOMO-LUMO gap:',Ggap,' for mu=',mu
      endif

      !if (DEBUG_RH_MU_E) then
      !   do i = 1,cfg_nits_debug  
      !     if (cfg_its_debug(i) == stat_current_iteration + 1) then
      !       !** Find connection between mu and various stuff
      !       call DEBUG_MU_DELTA_E_ORB(diag,queue,F,S,D,H1,Epred_0,Eorb_idem,EHFstart,cfg_mu_max_debug,muopt)
      !       WRITE(diag%LUPRI,*) 'Optimal level shift:',muopt
      !       exit
      !     endif
      !   enddo
      !endif
     
      if (diag%info_rh_gradient) then
         call mat_init(grad,ndim,ndim)
         !** ROOTHAN gradient before minimizing: 2SDbarF(Dbar) - 2F(Dbar)DbarS
         call get_AO_gradient(F,D,S,grad)
         !||G||^2 = TrGSGS
         WRITE(diag%LUPRI,'("Roothan gradnorm before",f18.9,29x," GRAD")')&
              & util_Snorm(grad,S)
   
         !** ROOTHAN gradient after minimizing: 2SDF(Dbar) - 2F(Dbar)DS
         call get_AO_gradient(F,Dnew,S,grad)
         !||G||^2 = TrGSGS
         WRITE(diag%LUPRI,'("Roothan gradnorm after",f18.9,29x," GRAD")')&
              & util_Snorm(grad,S)
         call mat_free(grad)
      endif

   end subroutine print_and_stat

!> \brief Calculate energy of idempotent density matrix
!> \author L. Thogersen
!> \date October 2005
!> \param Didem Idempotent density matrix
   real(realk) function pEHF(Didem)
      use dal_interface
      implicit none
      type(Matrix), intent(inout) :: Didem
      type(Matrix) :: Fcur 
      real(realk)  :: EHFcur
      
!      call mat_init(Fcur,Didem%Nrow,Didem%Ncol)
      !** Find new Fock matrix
      pEHF = 0.0E0_realk
      call lsquit('replace call di_get_fock(Didem,Fcur,EHFcur) with this quit statement',-1)      
!      pEHF = EHFcur
!      call mat_free(Fcur)

   end function pEHF

!> \brief Solves -SDSC = SCe giving the idempotent Didem = <C|C> corresponding to D.
!> \author L. Thogersen
!> \date October 2005
!> \param diag Used to store info about diagonalization
!> \param Dmu Density matrix
!> \param S Overlap matrix
!> \param Didem Idempotent density matrix
!> \param Cmo MO coefficients
   subroutine GET_Didem(diag,Dmu,S,Didem,Cmo)
      implicit none
      type(DiagItem), intent(inout) :: diag
      type(Matrix), intent(in) :: Dmu, S
      type(matrix)  :: Didem   !output
      type(matrix), optional :: Cmo ! output
      type(Matrix)  :: DS, SDmuS
      real(realk), pointer :: eig(:)
      integer     :: i,Nlin,j,k

      call mat_init(DS,S%nrow,S%ncol)
      call mat_init(SDmuS,S%nrow,S%ncol)
      call mat_mul(Dmu,S,'n','n',1E0_realk,0E0_realk,DS)
      call mat_mul(S,DS,'n','n',-1E0_realk,0E0_realk,SDmuS)
      if(PRESENT(Cmo)) THEN
         call mem_alloc(eig,S%nrow*2) ! allow for unrestricted
         call mat_diag_f(SDmuS,S,eig,Cmo)
         call mat_density_from_orbs(Cmo,Didem,diag%nocc,diag%nocca,diag%noccb)
         if (diag%debug_idempotency) then
            WRITE(diag%LUPRI,*) 'SDS eigenvalues'
            Nlin = MOD(S%nrow,7)
            j = 1
            do i = 1,Nlin
               WRITE(diag%LUPRI,'(7f10.5)') (eig(k),k=j,j+6)
               j = j+7
            enddo
            WRITE(diag%LUPRI,'(7f10.5)') (eig(k),k=j,S%nrow)
         endif
         call mem_dealloc(eig)
      ELSE
         call diag_get_density(diag,SDmuS,S,0.0E0_realk,Didem)
      END IF
      call mat_free(DS)
      call mat_free(SDmuS)
    end subroutine GET_Didem

!***************
!* DEBUG STUFF *
!***************
!!> \brief DEBUG: Print a table of values connected to a level shift mu 
!!> \author L. Thogersen
!!> \date October 2005
!!> \param diag Used to store info about diagonalization
!!> \param queue Contains density and Fock/KS matrices from previous SCF iterations
!!> \param F Fock/Kohn-Sham matrix
!!> \param S Overlap matrix
!!> \param D Density matrix
!!> \param H1 One-electron Hamiltonian
!!> \param Epred_0 Debug variable?
!!> \param Eorb_idem Debug variable?
!!> \param EHFidem Debug variable?
!!> \param mumax Debug variable?
!!> \param muopt Debug variable?
!!>
!!> (mu,dEorb(mu)) = (mu,Eorb(mu)-Eorb_start) || (mu,dEHF(mu)) = (mu,EHF(mu) - EHF_start)
!!> (mu,dEpred(mu)) = (mu,Epred(mu) - Epred_start)
!!>
!   subroutine debug_MU_DELTA_E_ORB(diag,queue,F,S,D,H1,Epred_0,Eorb_idem,EHFidem,mumax,muopt)
!      use scf_stats
!      use dal_interface
!      use fock_evaluator
!      implicit none
!      type(DiagItem),intent(in)          :: diag
!      TYPE(util_HistoryStore),intent(in) :: queue
!      type(Matrix), intent(in) :: F,S,D,H1
!      real(realk), intent(in) :: EHFidem, mumax,Epred_0,Eorb_idem
!      real(realk), intent(out) :: muopt
!      type(Matrix) :: Dnew, Fnew 
!      real(realk) :: mu,EHFslut,Eorb_slut,minpr,dEpred,gap,dmu,EHFold,dmuold,mochange,ratio, &
!                   &  EHFslut_diag, Eorb_slut_diag
!      integer :: n
!      logical :: hitlimit, look_for_opt
!
!      call mat_init(Dnew,D%nrow,D%ncol)
!      call mat_init(Fnew,F%nrow,F%ncol)
!      dmu = MAX(mumax/20E0_realk,0.05E0_realk)
!      n = 0
!      mu = mumax
!      hitlimit = .false.
!      muopt = 0.0E0_realk
!      look_for_opt = .true.
!      EHFold = 0.0E0_realk
!      debug_dd_dont_change_mu = .true.
!      debug_dd_limit_found = .false.
!      if (CFG_density_method == CFG_F2D_ROOTHAN) then
!         WRITE(diag%LUPRI,'("#   mu,        dERH,         dEpred,          dEHF,        gap,  TrDnewSDoldS,mochange,ratio_dorth, it",&
!                      &i3,"greprh")')stat_current_iteration+1 
!         do 
!            if (hitlimit) exit
!            if (mu < 0.0E0_realk) then
!              !exit
!              mu = 0.0E0_realk
!              hitlimit = .true.
!            endif
!            call diag_Get_density(F,S,mu,Dnew)
!            gap = Ggap
!            mochange = Gmochange
!            !Find ERH_mu to use for dERH = ERH_mu - ERH_0
!            Eorb_slut = 2E0_realk*mat_dotproduct(F,Dnew)
!            if (queue%used_entries > 0) then
!              !Find dEpred
!              dEpred = roothan_Epred(queue,H1,S,Dnew) - Epred_0
!              !Find ratio = ||Dorth||/||D||
!              call ratio_Dorth_D(queue,S,Dnew,ratio)
!            else
!              dEpred = 0.0E0_realk
!              ratio = 0.0E0_realk
!            endif
!            !Find minpr = TrDnewSDoldS
!            call Cchange(Dnew,D,S,minpr)
!            !** Find new Fock matrix and corresponding SCF energy
!            !To use for dEHF = EHFslut - EHF_0
!            call fck_get_fock(Dnew,Fnew,EHFslut)
!            !find the optimal mu
!            if (look_for_opt) then
!              if (EHFslut < EHFold) then
!                EHFold = EHFslut            
!              else
!                muopt = mu + dmuold
!                look_for_opt = .false.
!              endif
!            endif
!            WRITE(diag%LUPRI,'("#",f8.3,3e15.5,2f11.6,2f9.4,"          greprh")') &
!                 & mu,Eorb_slut-Eorb_idem,dEpred,EHFslut-EHFidem,gap,minpr,mochange,ratio
!
!            !** Update mu
!            mu = mu - dmu
!            dmuold = dmu
!            n = n + 1
!            if (n*dmu >= mumax/4E0_realk) then
!              n = 0
!              dmu = MAX(dmu/2.0E0_realk,0.15E0_realk)
!            endif
!         enddo
!         WRITE(diag%LUPRI,*) 'muopt = ',muopt,'   greprh'
!      else if (CFG_density_method == CFG_F2D_DIRECT_DENS) then
!         WRITE(diag%LUPRI,'("#   mu,  dERH(direcdens), dEHF(directdens), dERH(diag),    dEHF(diag), &
!                    & xmax(dd),    xnorm,   Epred,   it  ",i3,"  greprh")') stat_current_iteration+1 
!         do 
!            CFG_lshift = Diag_lshift_none
!            CFG_density_method = CFG_F2D_DIRECT_DENS
!            if (debug_dd_limit_found) exit
!            !if (mu < 0.6E0_realk*abs(cfg_mu_dd)) exit !cfg_mu_dd is the necessary damping found, 
!                                                !it doesn't make sense to go much below it 
!            if (hitlimit) exit
!            if (mu < 0.0E0_realk) then
!              !exit
!              mu = 0.0E0_realk
!              hitlimit = .true.
!            endif
!!##### First do direct dens.....######################################################
!            write (diag%lupri,*) 'mu before direcdens', mu
!            call diag_Get_density(F,S,mu,Dnew)
!            write (diag%lupri,*) 'mu after direcdens', mu
!            !Find ERH_mu to use for dERH = ERH_mu - ERH_0
!            Eorb_slut = 2E0_realk*mat_dotproduct(F,Dnew)
!            !if (queue%used_entries > 0) then
!            !  !Find dEpred
!            !  dEpred = roothan_Epred(queue,H1,S,Dnew) - Epred_0
!            !  !Find ratio = ||Dorth||/||D||
!            !  call ratio_Dorth_D(queue,S,Dnew,ratio)
!            !else
!            !  dEpred = 0.0E0_realk
!            !  ratio = 0.0E0_realk
!            !endif
!            !Find minpr = TrDnewSDoldS
!            !call Cchange(Dnew,D,S,minpr)
!            !** Find new Fock matrix and corresponding SCF energy
!            !To use for dEHF = EHFslut - EHF_0
!            call fck_get_fock(Dnew,Fnew,EHFslut)
!            !find the optimal mu
!            !if (look_for_opt) then
!            !  if (EHFslut < EHFold) then
!            !    EHFold = EHFslut            
!            !  else
!            !    muopt = mu + dmuold
!            !    look_for_opt = .false.
!            !  endif
!            !endif
!!##### ..... then do diagonalization ######################################################
!            CFG_density_method = CFG_F2D_ROOTHAN
!            CFG_lshift = diag_lshift_dorth
!            write (diag%lupri,*) 'mu before diag', mu
!            call diag_Get_density(F,S,mu,Dnew)
!            write (diag%lupri,*) 'mu after diag', mu
!            gap = Ggap
!            mochange = Gmochange
!            !Find ERH_mu to use for dERH = ERH_mu - ERH_0
!            Eorb_slut_diag = 2E0_realk*mat_dotproduct(F,Dnew)
!            !if (queue%used_entries > 0) then
!            !  !Find dEpred
!            !  dEpred = roothan_Epred(queue,H1,S,Dnew) - Epred_0
!            !  !Find ratio = ||Dorth||/||D||
!            !  call ratio_Dorth_D(queue,S,Dnew,ratio)
!            !else
!            !  dEpred = 0.0E0_realk
!            !  ratio = 0.0E0_realk
!            !endif
!            !Find minpr = TrDnewSDoldS
!            call Cchange(Dnew,D,S,minpr)
!            !** Find new Fock matrix and corresponding SCF energy
!            !To use for dEHF = EHFslut - EHF_0
!            call fck_get_fock(Dnew,Fnew,EHFslut_diag)
!            !find the optimal mu
!            !if (look_for_opt) then
!            !  if (EHFslut_diag < EHFold) then
!            !    EHFold = EHFslut_diag            
!            !  else
!            !    muopt = mu + dmuold
!            !    look_for_opt = .false.
!            !  endif
!            !endif
!            print *, 'Epred 0:',  Epred_0
!            print *, 'Epred ', debug_dd_epred
!            WRITE(diag%LUPRI,'("#",f8.3,4e15.5,2e12.4,e18.8, "          greprh")') &
!                 & mu,Eorb_slut-Eorb_idem,EHFslut-EHFidem,Eorb_slut_diag-Eorb_idem,EHFslut_diag-EHFidem,&
!                 & debug_dd_maxelm, debug_dd_xnorm, Epred_0 - debug_dd_epred
!            !** Update mu
!            mu = mu - dmu
!            dmuold = dmu
!            n = n + 1
!            if (n*dmu >= mumax/4E0_realk) then
!              n = 0
!              dmu = MAX(dmu/2.0E0_realk,0.05E0_realk)
!            endif
!         enddo
!      endif
!
!      debug_dd_dont_change_mu = .false.
!      call mat_free(Dnew)
!      call mat_free(Fnew)
!   end subroutine debug_MU_DELTA_E_ORB

!!$!> \brief Calculate minpr = TrDoldSDnewS/sqrt(TrDoldSDoldS*TrDnewSDnewS)
!!$!> \author L. Thogersen
!!$!> \date October 2005
!!$!> \param Didem Idempotent density matrix
!!$!> \param Dnonidem Non-idempotent density matrix
!!$!> \param S Overlap matrix
!!$!> \param minpr = TrDoldSDnewS/sqrt(TrDoldSDoldS*TrDnewSDnewS)
!!$   SUBROUTINE Cchange(Didem, Dnonidem, S, minpr)
!!$     !     implicit none
!!$     type(Matrix), intent(in) :: Didem, Dnonidem, S
!!$     real(realk), intent(out) :: minpr
!!$     type(Matrix) :: SDnewS,SDnew !, pointer
!!$     real(realk) :: nrm1, nrm2
!!$ 
!!$     call mat_init(SDnew, S%nrow,S%ncol)
!!$     call mat_init(SDnewS,S%nrow,S%ncol)
!!$     call mat_mul(S,Dnonidem,'n','n',1E0_realk,0E0_realk,SDnew)
!!$     call mat_mul(SDnew,S,'n','n',1E0_realk,0E0_realk,SDnewS)
!!$     !TrDnewSDoldS
!!$     minpr = mat_dotproduct(Didem,SDnewS)
!!$     !TrDoldSDoldS
!!$     nrm1  = mat_dotproduct(Dnonidem,SDnewS)
!!$     call mat_mul(S,Didem,'n','n',1E0_realk,0E0_realk,SDnew)
!!$     call mat_mul(SDnew,S,'n','n',1E0_realk,0E0_realk,SDnewS)
!!$     !TrDnewSDnewS
!!$     nrm2 = mat_dotproduct(Didem,SDnewS)
!!$     minpr = minpr / sqrt(nrm1*nrm2)
!!$     call mat_free(SDnew)
!!$     call mat_free(SDnewS)
!!$   end subroutine Cchange

!> \brief Diagonalize Fock/KS matrix to get density.
!> \author S. Host
!> \date January 2010
!> \param diag Used to store info about diagonalization
!> \param F Fock/Kohn-Sham matrix
!> \param S Overlap matrix
!> \param mu Level shift
!> \param Dnew New density matrix
!> \param queue Contains density and Fock/KS matrices from previous SCF iterations
   subroutine diag_get_density(diag,F,S,mu,Dnew,queue)
     !Made a local reduced copy of get_density to avoid circular dependency. /Stinne
     implicit none
     type(diagItem),intent(inout) :: diag !eHOMO is added in diag_f_to_get_d but otherwise intent(in)
     type(Matrix), intent(in)  :: F, S
     real(realk), intent(in)   :: mu
     type(Matrix), intent(inout) :: Dnew  !output
     type(Matrix) :: Fdamp
     logical        :: shiftedF
!test
     real(realk) :: orbE,gradnorm
     type(Matrix) :: grad
     TYPE(util_historyStore), optional :: queue

     shiftedF = .false.
     if (diag%cfg_lshift /= diag_lshift_none .and. ABS(mu) > 1E-8_realk) then
       !Evaluate the shifted Fock matrix
       ! Fdamp = F - muSDS
        call mat_init(Fdamp,S%nrow,S%nrow)
        CALL mat_add(1E0_realk,F,-mu, SDS, Fdamp)
        shiftedF = .true.
     endif
     ! Diagonalization
     if(shiftedF) then
       call diag_f_to_get_d(diag,Fdamp,S,Dnew)
     else
       call diag_f_to_get_d(diag,F,S,Dnew)
     endif
     if (shiftedF) then
      call mat_free(Fdamp)
     endif
   end subroutine diag_get_density

!!$!> \brief Scale virtual orbitals.
!!$!> \author L. Thogersen
!!$!> \date 2003
!!$!> \param S Overlap matrix
!!$!> \param H1 One-electron Hamiltonian
!!$!> \param F Fock/Kohn-Sham matrix
!!$!> \param D Density matrix
!!$!>
!!$!> Create a N-1 potential for the virtual orbitals
!!$!> allowing for more 'occupied-like' virtual orbitals
!!$!> good when configuration shifts are needed!!!
!!$!> ** F_mod = F + (n-1/n - 1)Q^T (F - H1) Q
!!$!>
!!$   subroutine diag_scale_virt_fock(S,H1,F,D,cfg_nocc)
!!$     implicit none
!!$     type(Matrix), intent(in) :: S,H1,D
!!$     type(Matrix), intent(inout) :: F
!!$     type(matrix) :: Q,wrk1,wrk2
!!$     real(realk) :: konst
!!$     integer :: ndim,cfg_nocc
!!$     
!!$     ndim = S%nrow
!!$     call mat_init(Q,ndim,ndim)
!!$     call mat_init(wrk1,ndim,ndim)
!!$     call mat_init(wrk2,ndim,ndim)
!!$
!!$     !** wrk2 = F - h1
!!$     call mat_add(1E0_realk,F,-1E0_realk,H1,wrk2)
!!$     !** Q = 1 - DS
!!$     call mat_mul(D,S,'n','n',1E0_realk,0E0_realk,wrk1)
!!$     call mat_add_identity(1E0_realk,-1E0_realk,wrk1,Q)
!!$     !** wrk1 = (F - h1) Q = wrk2 Q
!!$     call mat_mul(wrk2,Q,'n','n',1E0_realk,0E0_realk,wrk1)
!!$     !** konst = -1/n
!!$     konst = -1E0_realk/(2E0_realk*cfg_nocc)
!!$     !** F_mod = F + konst Q^T (F - H1) Q = F + konst Q^T wrk1
!!$     call mat_mul(Q,wrk1,'T','n',konst,1E0_realk,F)
!!$
!!$     call mat_free(Q)
!!$     call mat_free(wrk1)
!!$     call mat_free(wrk2)
!!$   end subroutine diag_scale_virt_fock

!Returns the overlap of the new MO with the smallest overlap to the old MOs
!!$   subroutine util_min_MO_overlap(cfg_nocc,S,Cold,Cnew,noverlap,overlap)
!!$      implicit none
!!$      integer, intent(in)      :: cfg_nocc
!!$      type(Matrix), intent(in) :: S,Cold,Cnew
!!$      integer, intent(in) :: noverlap
!!$      real(realk), intent(out) :: overlap(noverlap)
!!$      type(Matrix) :: ColdTS,ColdTSCnew
!!$      real(realk) :: scr(cfg_nocc),xxx
!!$      real(realk) :: overlap_x
!!$      integer :: i,j
!!$
!!$      call mat_init(ColdTS,S%nrow,S%ncol)
!!$      call mat_init(ColdTSCnew,S%nrow,S%ncol)
!!$      call mat_mul(Cold,S,'t','n',1E0_realk,0E0_realk,ColdTS)
!!$      call mat_mul(ColdTS,Cnew,'n','n',1E0_realk,0E0_realk,ColdTSCnew)
!!$      if (noverlap < 1) STOP ' wrong input for util_min_MO_overlap'
!!$      if (noverlap == 1) then
!!$         overlap = mat_column_norm(ColdTSCnew,1,1,CFG_NOCC)
!!$         do i = 2,CFG_NOCC
!!$           !Returns the sum of the elements in th ith column squared
!!$           !with the row limits 1:CFG_NOCC
!!$           overlap_x = mat_column_norm(ColdTSCnew,i,1,CFG_NOCC)
!!$           overlap = MIN(overlap,overlap_x)
!!$         enddo
!!$      else
!!$         do i = 1,CFG_NOCC
!!$           !Returns the sum of the elements in th ith column squared
!!$           !with the row limits 1:CFG_NOCC
!!$           scr(i) = mat_column_norm(ColdTSCnew,i,1,CFG_NOCC)
!!$         enddo
!!$         !sort them such that the lowest is first
!!$         !FIXME: use better sorting algorithm
!!$         do i = 1,cfg_nocc
!!$           do j = 1,cfg_nocc-1
!!$             if (scr(j) > scr(j+1)) then
!!$               xxx = scr(j+1)
!!$               scr(j+1) = scr(j)
!!$               scr(j) = xxx
!!$             endif
!!$           enddo
!!$         enddo
!!$         overlap(1:noverlap) = scr(1:noverlap)
!!$      endif
!!$      call mat_free(ColdTS)
!!$      call mat_free(ColdTSCnew)
!!$   end subroutine util_min_MO_overlap

END MODULE diagonalization
