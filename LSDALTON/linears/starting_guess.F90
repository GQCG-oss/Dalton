!> @file 
!> Contains initial guess module.

!> \brief Get the initial guess for SCF optimization.
!> \author L. Thogersen. Documented by S. Host.
!> \date 2003
!>
!> Find a starting guess for the density matrix. Options are:
!> - eigenvectors of H1
!> - Hueckel guess
!> - from atomic densities
!>
MODULE initial_guess
  use configurationType
  use files
  use matrix_module
  use matrix_operations
  use matrix_operations_aux
  use linsca_debug
  use matrix_util
  use typedefTYPE, only: LSSETTING, lsitem
  use typedef
  use decompMod
  use dal_interface
  use memory_handling
  use precision
  implicit none
  private
  public :: starting_guess_entry, get_initial_dens!, asymmetrize_starting_guess

CONTAINS

!> \brief Wrapper for starting guess for density matrix.
!> \author L. Thogersen
!> \date 2002
subroutine get_initial_dens(H1,S,D,ls,config)
   implicit none
   !> One-electron Hamiltonian
   type(Matrix),intent(inout)     :: H1
   !> Overlap matrix
   type(Matrix),intent(inout),target :: S
   !> Density matrix (output)
   type(Matrix),intent(inout)     :: D(1)
   !> Contains setting for integrals
   type(lsitem) :: ls
   !> Contains all info about configuration/settings for SCF calculation
   type(ConfigItem) :: config
   type(Matrix) :: C
   real(realk) :: trace, Nelectrons
   integer :: restart_lun, idum, ldum, i, j
   logical :: restart_from_dens, restart_from_cmo, purify_failed, gcbasis,dens_exsist,cmo_exsist, start_from_guess, restart_l2
   logical lcv_exsist,lcv_on_file,OnMaster 
   real(realk),parameter :: THRNEL = 1E-3_realk
   data lcv_on_file /.false./
   config%decomp%S => S !For ATOMS to work
   restart_from_dens = .false.
   restart_from_cmo = .false.
   OnMaster = .TRUE.
  !The if-statment below breaks linsca_uhf_stability and fck3_df_molgrad_linsca
  !testcases, therefore commented out. 31-01-2009, /filip.
  !IF (cfg_start_guess(1:5).eq.'FILE') THEN
   INQUIRE(file='dens.restart',EXIST=dens_exsist) 
   INQUIRE(file='cmo.restart',EXIST=cmo_exsist)

   restart_from_dens = dens_exsist.AND.config%diag%cfg_restart
   restart_from_cmo  = cmo_exsist.AND.config%diag%cfg_restart
   restart_l2        = config%diag%cfg_redo_l2.and.restart_from_dens    

   start_from_guess  = (.not.(restart_from_dens.or.restart_from_cmo)).or. &
                     & (restart_l2)

   !If we want to redo L2, we need to move dens.restart so that
   !it is not overwritten by L2
   if (restart_l2) then 
#ifdef SYS_AIX
    call rename('dens.restart\0','fdens.restart\0')
#else
    call rename('dens.restart','fdens.restart')
#endif 
   endif

   if (start_from_guess) call starting_guess_entry(H1,D,S,ls,config) 

   if (restart_l2)  then 
#ifdef SYS_AIX
    call rename('fdens.restart\0','dens.restart\0')
#else
    call rename('fdens.restart','dens.restart')
#endif
   endif
   
   INQUIRE(file='lcv_basis.restart',EXIST=lcv_exsist)

  !ENDIF
   !IF(cfg_run_Exgrad .AND. .NOT.cfg_rsp_grad_purify) restart_from_dens = .false.
   !if(cfg_rsp_run_grad) & 
   !           & restart_from_dens = cfg_rsp_grad_drestart.and.restart_from_dens
   if (restart_from_dens) then
     restart_lun = -1  !initialization
     call lsopen(restart_lun,'dens.restart','OLD','UNFORMATTED')
     rewind restart_lun
     call mat_read_from_disk(restart_lun,D(1),OnMaster)
     call mat_read_info_from_disk(restart_lun,gcbasis)
     call lsclose(restart_lun,'KEEP')
     IF (.NOT.(config%optinfo%optimize.AND.(config%optinfo%ItrNmr.GT. 0))) THEN
       WRITE(config%lupri,*)
       WRITE(config%lupri,*) '*** RESTART FROM DENSITY ON DISK - READ FROM dens.restart  ***'
       WRITE(config%lupri,*)
       WRITE(*,*)
       WRITE(*,*) '*** RESTART FROM DENSITY ON DISK - READ FROM dens.restart  ***'
       WRITE(*,*)
     ENDIF

     if (lcv_exsist) then
        ! Print header
        WRITE(config%lupri,*)
        WRITE(config%lupri,*) '*** LOADING LEVEL 2 BASIS FROM DISK - READ FROM lcv_basis.restart  ***'
        WRITE(config%lupri,*)

        WRITE(*,*)
        WRITE(*,*) '*** LOADING LEVEL 2 BASIS FROM DISK - READ FROM lcv_basis.restart  ***'
        WRITE(*,*)

        !Load in saved lcv basis
        !Allocate space and make decomposition aware of the basis
        call mat_init(config%decomp%lcv_CMO,D(1)%nrow,D(1)%ncol)
        config%decomp%decompMatInit_lcv_CMO = .TRUE.
        config%decomp%lcv_basis = .true.
 
        ! Load the file, check if that is plain L2 basis or localized LCV basis
        restart_lun = -1;
        call lsopen(restart_lun,'lcv_basis.restart','unknown','UNFORMATTED')
        call mat_read_from_disk(restart_lun,config%decomp%lcv_Cmo,OnMaster)
        call mat_read_info_from_disk(restart_lun,lcv_on_file)
        call lsclose(restart_lun,'KEEP')
     endif

     ! Quit in case orbitals on file are not localized lcv basis and we are
     ! trying to restart lcm calculation for which lcv is needed.
     if (config%opt%cfg_start_guess=='TRILEVEL') then
        if (config%decomp%cfg_lcm .and. (.not.lcv_on_file))then
           call lsquit("CAN NOT RESTART TRILEVEL .LCM calculation, LCV basis NOT ON FILE",config%lupri)
        endif
     endif

     if (gcbasis .and. .not. config%decomp%cfg_gcbasis) then
        WRITE(config%lupri,*) 'Your dens.restart was constructed using the grand-canonical (GC) basis,'
        WRITE(config%lupri,*) 'while your LSDALTON.INP uses the standard basis.'
        WRITE(config%lupri,*) 'The GC basis is default unless you use a dunnings basis set,'
        WRITE(config%lupri,*) 'or you specify .NOGCBASIS under *GENERAL'
        WRITE(config%lupri,*) 'Either contruct a new dens.restart or '
        WRITE(config%lupri,*) 'remove .NOGCBASIS or add .FORCEGCBASIS under *GENERAL'
        call lsquit('Calculation in standard basis, dens.restart in GC basis!',config%lupri)
     else if (config%decomp%cfg_gcbasis .and. .not. gcbasis) then
        IF(config%decomp%cfg_transformrestart)THEN
           WRITE(config%lupri,*) 'Your dens.restart was constructed using the standard basis'
           WRITE(config%lupri,*) 'We transform to the grand-canonical (GC) basis.'
           call AO2GCAO_transform_matrixD(D(1),ls%setting,config%lupri)
        ELSE
           WRITE(config%lupri,*) 'Your dens.restart was constructed using the standard basis, while your'
           WRITE(config%lupri,*) 'LSDALTON.INP uses the grand-canonical (GC) basis.'
           WRITE(config%lupri,*) 'The GC basis is default unless you use a dunnings basis set,'
           WRITE(config%lupri,*) 'or you specify .NOGCBASIS under *GENERAL'
           WRITE(config%lupri,*) 'Either contruct a new dens.restart or use GC basis!'
           WRITE(config%lupri,*) 'Alternativly you can add the keyword '
           WRITE(config%lupri,*) '.TRANSFORMRESTART'
           call lsquit('Calculation in GC basis, dens.restart in standard basis!',config%lupri)
        ENDIF
     else if (config%decomp%cfg_gcbasis .and. gcbasis) then
        WRITE(config%lupri,*) 'Basis check ok: Using GC basis consistently'
     else if (.not. config%decomp%cfg_gcbasis .and. .not. gcbasis) then
        WRITE(config%lupri,*) 'Basis check ok: Using standard basis consistently'
     else
        call lsquit('Basis check is messed up!!',config%lupri)
     endif
     IF (config%diag%CFG_restart .AND. &
          & ((config%optinfo%optimize .OR. config%dynamics%do_dynamics)&
          & .OR.config%diag%CFG_purifyrestart)) THEN
        !Using the density matrix from the previous geometry iteration in the current
        !geometry iteration speeds up geometry optimization significantly. 
        !Note that in some (early) geometry iterations this might lead to problems in SCF 
        !convergence (or to a loss of electrons), because the overlap matrix has
        !changed substantially and, therefore, the trace and idempotency conditions are 
        !not satisfied.  
        !Below this problem is solved by McWeeny purification. 
        !FIXME:It might require some changes for the open-shell systems!!!
        write(config%lupri,*)'Performing McWeeny purification for molecular gradient run.'

        call McWeeney_purify(S,D(1),purify_failed)
        !If purification failed, fall back to a default initial guess!
        if (purify_failed) then
         write(config%lupri,*) 'Warning:McWeeny purification failed. We fall back to an initial guess!'
         call starting_guess_entry(H1,D,S,ls,config)
        else
           ! In some cases the purification might change the occupation drastically
           ! enough to lead to a loss/gain of electrons. Therefore, we calculate here 
           ! Tr(DS) and compare with the number of electrons. If the difference is large 
           ! we instead start from a new initial guess
           Nelectrons = FLOAT(2*config%decomp%nocc + config%decomp%nactive)
           trace =  mat_dotproduct(D(1),S)
           if ( abs(2.0*trace-Nelectrons) .gt. THRNEL*Nelectrons) then
             write(config%lupri,*)'Warning: Descrepancy between 2.0*Tr(DS) after purification and Nelectrons '
             write(config%lupri,*)'2.0*Tr(DS) = ',2.0*trace,' Nelectrons = ', Nelectrons
             write(config%lupri,*)'We fall back to an initial guess!!!'
             call starting_guess_entry(H1,D,S,ls,config)
           endif !2*Tr(DS)-Nelectrons > THRNEL
        endif !cfg_dd_purify_failed
     endif !purify
   elseif (restart_from_cmo) then
     restart_lun = -1  !initialization
     call lsopen(restart_lun,'cmo.restart','OLD','UNFORMATTED')
     rewind restart_lun
     call mat_init(C,D(1)%nrow,D(1)%ncol)
     call mat_read_from_disk(restart_lun,C,OnMaster)
     call mat_density_from_orbs(C,D(1),config%decomp%nocc,config%decomp%nocca,config%decomp%noccb)
     call mat_free(C)
     call lsclose(restart_lun,'KEEP')
     WRITE(config%lupri,*)
     WRITE(config%lupri,*) '*** RESTART FROM ORBITALS ON DISK - READ FROM cmo.restart  ***'
     WRITE(config%lupri,*)
   endif


   config%decomp%S => S !For TRILEVEL to work

   if (config%opt%DEBUG_CONVERT) then
      call debug_convert_density(config%opt,D(1))
   endif

   IF (ls%setting%scheme%intprint.GE.2) write(config%lupri,'(A46,F18.8)') &
     & 'Initial-density-matrix dot product:',mat_dotproduct(D(1),D(1))

end subroutine get_initial_dens

!> \brief Branches out and chooses the proper initial guess (as requested by default or input)
!> \author L. Thogersen
!> \date 2003
  subroutine starting_guess_entry(H1,D,S,ls,config)
    implicit none
    !> One-electron Hamiltonian
    type(matrix),intent(inout)     :: H1,S
    !> Initial density matrix (output)
    type(matrix),intent(inout)     :: D(1)
    !> Contains info about integrals (?)
    type(lsitem),intent(inout)     :: ls
    !> Contains all info about configuration/settings for SCF calculation
    type(ConfigItem),intent(inout) :: config
    integer, parameter :: sguess_h1=1, sguess_hueck=2, sguess_atden=3
    integer :: ndmat
    logical :: do_huckel
    interface
       subroutine trilevel_start(D,ls,config)
         use typedeftype
         use matrix_module
         use configurationType
         TYPE(lsitem),target :: ls
         Type(Matrix) :: D(1)
         type(ConfigItem) :: config
       end subroutine trilevel_start
    end interface
    interface
       subroutine atoms_start(config,D,H1,S,ls,ndmatalloc)
         use typedeftype
         use matrix_module
         use configurationType
         type(ConfigItem),intent(in) :: config
         TYPE(lsitem),intent(inout) :: ls
         Type(Matrix),target        :: H1
         Type(Matrix),intent(inout) :: D(ndmatalloc),S
         integer,intent(in)         :: ndmatalloc
       end subroutine atoms_start
    end interface
    !Huckel doesn't work for unrestricted - not sure why! /Stinne
    if (config%decomp%cfg_unres) then
       if (config%opt%cfg_start_guess=='HUCKEL') then
          write(config%lupri,*) 'FALLBACK: Huckel guess does not work with open shell systems'
          write(config%lupri,*) '- will do H1DIAG instead (if you are the Huckel author, please fix it!)'
          config%opt%cfg_start_guess='H1DIAG'
       endif
    endif
    write(config%lupri,*)
    if (config%opt%cfg_start_guess=='H1DIAG') then
       write(*,*) 'Optimize first density for H1 operator'
       write(config%lupri,*) 'Optimize first density for H1 operator'
       write(config%lupri,*)
       call starting_guess_h1(config,H1,D(1),S)
    else if (config%opt%cfg_start_guess=='HUCKEL') then
        call lsquit('Huckel guess not implemented',config%lupri)
    !   write(*,*) 'Take first density from Hueckel guess'
    !   write(config%lupri,*) 'Take first density from Hueckel guess'
    !   call starting_guess_hueckel(config%decomp,H1,D)
    !case(sguess_atden)
    !   write(*,*) 'Take first density from atomic densities'
    !   write(lupri,*) 'Take first density from atomic densities'
    !   call starting_guess_atden(D)
    else if (config%opt%cfg_start_guess=='ATOMS') then
       write(config%lupri,*) 'First density: Atoms in molecule guess'
       write(config%lupri,*)
       ndmat = 1
       call atoms_start(config,D,H1,S,ls,ndmat)
    else if (config%opt%cfg_start_guess=='TRILEVEL') then
       write(config%lupri,*) 'First density: Trilevel procedure '
       write(config%lupri,*)
       call trilevel_start(D,ls,config)
    else if (config%opt%cfg_start_guess=='LINCOMB') then
       write(config%lupri,*) 'Starting guess is linear combination of densities on disk'
       write(config%lupri,*)
       call lsquit('starting_guess_lincomb removed',-1)
!       call starting_guess_lincomb(D(1),config%opt)
    else
       write(*,*) 'Optimize first density for H1 operator'
       write(config%lupri,*) 'Optimize first density for H1 operator'
       write(config%lupri,*)
       call starting_guess_h1(config,H1,D(1),S)
    !ELSE
    !   ! Option to start from the fitted density. Works only in Coulomb,
    !   ! i.e. the Hartree approximation. Could prove useful for projecting
    !   ! from a small to a larger (normal) basis via the fitted density. 
    !   ! Currently turned off. \SR
    !   CALL SET_ITER(0)
    !   !Stinne, Simen jan. 2007: It is not necessary to have F as an argument here.
    !   call starting_guess_fit_density(S,D)
    ENDIF

  end subroutine starting_guess_entry

! commentet out by Thomas Kjaergaard - no test case - not tested ....  
!!$!> \brief Asymmetrizes the starting guess when by mixing HOMO and LUMO
!!$!> \author C. Nygaard
!!$!> \date March 2010
!!$!> \param Cmo The MO orbital coefficients 
!!$  subroutine asymmetrize_starting_guess (Cmo, decomp)
!!$
!!$  implicit none
!!$
!!$  !I/O:
!!$  type(decompItem), intent(in) :: decomp
!!$  type(matrix), intent(inout)  :: Cmo
!!$  !Other
!!$  type(matrix)                 :: homo, lumo
!!$  integer                      :: i, ndim
!!$  real(realk)                  :: elms(2)
!!$
!!$  ndim = Cmo%nrow
!!$
!!$  call mat_init (homo, ndim, 1)
!!$  call mat_init (lumo, ndim, 1)
!!$  
!!$  write (decomp%lupri, *) "Assymetrized starting guess for UHF chosen"
!!$  write (decomp%lupri, *) "mixing HOMO and LUMO"
!!$
!!$  call mat_section (Cmo, 1, ndim, decomp%nocc, decomp%nocc, homo)
!!$  call mat_section (Cmo, 1, ndim, decomp%nocc+1, decomp%nocc+1, lumo)
!!$
!!$  call mat_mix_homolumo (homo, lumo)
!!$
!!$  do i=1,ndim
!!$    call mat_get_ab_elms (homo, i, 1, elms)
!!$    call mat_create_ab_elms (i, decomp%nocc, elms, Cmo)
!!$    call mat_get_ab_elms (lumo, i, 1, elms)
!!$    call mat_create_ab_elms (i, decomp%nocc+1, elms, Cmo)
!!$  enddo
!!$
!!$  call mat_free (homo)
!!$  call mat_free (lumo)
!!$
!!$  end subroutine asymmetrize_starting_guess


!> \brief Obtain initial guess by diagonalizing one-el. Hamiltonian
!> \author L. Thogersen
!> \date 2003
!> \param lupri Logical unit number for output
!> \param H1 The one-electron part of the Fock matrix
!> \param D The AO density matrix
  subroutine starting_guess_h1(config,H1,D,S)
    implicit none
    type(configItem),intent(in) :: config
    type(matrix), intent(in)    :: H1,S
    type(matrix),intent(inout)  :: D
    type(Matrix)                :: Cmo
    real(realk), pointer    :: eival(:)
    integer :: cycles

    call mem_alloc(eival,S%nrow*2) ! allow for unrestricted.
    call mat_init(Cmo,S%nrow,S%nrow)
    call mat_diag_f(H1,S,eival,Cmo)

!    if (config%decomp%cfg_unres .and. config%opt%cfg_asym) then
!      call asymmetrize_starting_guess (Cmo, config%decomp)
!    endif

    call mat_density_from_orbs(Cmo,D,config%decomp%nocc,config%decomp%nocca,config%decomp%noccb)
    call mem_dealloc(eival)
    call mat_free(Cmo)
  end subroutine starting_guess_h1

!!$  !> \brief Obtain initial guess from the fitted density.
!!$  !> \author S. Reine
!!$  !> \date 2005
!!$  subroutine starting_guess_fit_density(decomp,D)
!!$    implicit none
!!$    !> Contains matrices from OAO decomposition of overlap matrix
!!$    type(decompItem),intent(in) :: decomp
!!$    !> Initial density matrix (output)
!!$    type(matrix)             :: D
!!$    type(matrix)             :: F
!!$    type(matrix)             :: cmo
!!$    integer                  :: ndim
!!$    real(realk)              :: Etotal
!!$    logical :: getd, getdv
!!$    real(realk), pointer :: orb(:), eival(:)
!!$
!!$    ndim=D%nrow
!!$    call mat_init(F,ndim,ndim)
!!$    call lsquit('CALL di_get_fock(D,F,Etotal) replaced with a quit statement in starting_guess_fit_density',-1)
!!$    call mat_init(cmo,ndim,ndim)
!!$    call mem_alloc(eival,ndim*2) ! allow for unrestricted
!!$    call mat_diag_f(F,decomp%S,eival,cmo)
!!$    call mem_dealloc(eival)
!!$    call mat_density_from_orbs(cmo,D,decomp%nocc,decomp%nocca,decomp%noccb)
!!$    call mat_free(cmo)
!!$    call mat_free(F)
!!$  end subroutine starting_guess_fit_density

!!$  !> \brief Obtain initial guess from linear combination of saved densities.
!!$  !> \author S. Host
!!$  !> \date February 2010 
!!$  subroutine starting_guess_lincomb(D,opt)
!!$  implicit none
!!$      !> Initial density (output)
!!$      type(matrix), intent(inout) :: D
!!$      !> Contains general settings for SCF optimization
!!$      type(optItem), intent(in)   :: opt
!!$      type(matrix)                :: D1, D2
!!$      integer                     :: idum, ldum, D1lun, D2lun
!!$      logical                     :: D1_exists, D2_exists,OnMaster
!!$      OnMaster = .TRUE.
!!$      INQUIRE(file='D1',EXIST=D1_exists) 
!!$      INQUIRE(file='D2',EXIST=D2_exists)
!!$      if (.not. D1_exists) then
!!$         write(opt%lupri,*) 'File D1 must be present with .HESONLY'
!!$         call lsquit('File D1 must be present with .HESONLY',opt%lupri)
!!$      else if (.not. D2_exists) then
!!$         write(opt%lupri,*) 'File D2 must be present with .HESONLY'
!!$         call lsquit('File D2 must be present with .HESONLY',opt%lupri)
!!$      endif
!!$      call mat_init(D1,D%nrow,D%ncol)
!!$      call mat_init(D2,D%nrow,D%ncol)
!!$      D1lun = -1 ; D2lun = -1
!!$      CALL LSOPEN(D1lun,'D1','OLD','UNFORMATTED')
!!$      CALL LSOPEN(D2lun,'D2','OLD','UNFORMATTED')
!!$      call mat_read_from_disk(D1lun,D1,OnMaster)
!!$      call mat_read_from_disk(D2lun,D2,OnMaster)
!!$      call LSCLOSE(D1lun,'KEEP')
!!$      call LSCLOSE(D2lun,'KEEP')
!!$
!!$      call mat_add(opt%cfg_weight_param,D1,1.0E0_realk-opt%cfg_weight_param,D2,D)
!!$
!!$      call mat_free(D1)
!!$      call mat_free(D2)
!!$   end subroutine starting_guess_lincomb

END MODULE initial_guess


