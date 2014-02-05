!> @file 
!> Contains general routine to calculate energy
!> and it's geometric derivative
module Energy_and_deriv
!
use precision
use matrix_module, only: matrix
use TYPEDEFTYPE, only: lsitem
use configurationType, only: configitem
use TYPEDEF, only: typedef_set_default_setting, typedef_init_setting, &
     & typedef_free_setting
use memory_handling, only: mem_alloc,mem_dealloc
use scfloop_module, only: scfloop
use matrix_operations, only: mat_init, mat_free, mat_diag_f
use basis_type, only: free_basissetinfo
use ks_settings, only: ks_init_incremental_fock, ks_free_incremental_fock
use decompMod, only: decomp_shutdown, decomp_init, decomposition
use initial_guess, only: get_initial_dens
use dec_typedef_module, only: DECinfo
use lsdalton_fock_module, only: lsint_fock_data
use matrix_operations_aux, only: mat_density_from_orbs
use integralinterfaceMod, only: II_get_molecular_gradient,&
     & II_get_nucpot,II_get_overlap,II_get_h1
use lsdalton_rsp_mod,only: get_excitation_energy, GET_EXCITED_STATE_GRADIENT
use dec_main_mod!, only: get_total_mp2energy_from_inputs, get_mp2gradient_and_energy_from_inputs
use optimlocMOD, only: optimloc
use screen_mod, only: screen_free, screen_init
private
public :: Get_Energy, Get_Gradient, get_num_grad
!
contains
! Get energy: calculates energy in a general way
!> \currently used for energies at modified goemetries
!> \ covers SCF and DEC energies + dispersion for DFT
!> \author \latexonly V. Rybkin  \endlatexonly
!> \date 2012-27-02
  Subroutine Get_Energy(E,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
!
!
    !
    Implicit none
    Type(lsitem),target :: ls   ! General information,used only to get E and gradient
    Type(Matrix), intent(inout) :: F(1),D(1)       ! Fock,density
    Type(Matrix), intent(inout),target :: S  ! overlap matrices
    Type(Matrix), intent(inout),target :: H1 ! One electron matrix
    Type(Matrix), intent(inout) :: C       ! Orbitals
    Type(Matrix) :: CMO       ! Orbitals for Fock matrix dynamics
    Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk) :: E(1), PotNuc   ! Electronic energy and nuclear potential
    Real(realk) :: Eerr ! For DEC: Estimated intrinsic energy error
    Real(realk), allocatable :: eival(:)
!    Real(realk), pointer :: ExcitE(:)
    Integer :: NAtoms,i,lupri,luerr,nbast
    Real(realk) :: DUMMY(1,1),ExcitE
    Logical :: do_decomp,integraltransformGC


       Eerr = 0E0_realk
       nbast = D(1)%nrow
       do_decomp =(config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
            & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
            & config%decomp%cfg_check_converged_solution .or. &
            & config%decomp%cfg_rsp_nexcit > 0 .or. config%integral%locallink) 
       integraltransformGC = ls%setting%integraltransformGC
       if (do_decomp) then
          call decomp_shutdown(config%decomp)
       endif

       !call mat_zero(F)
       !call mat_zero(H1) 
       !call mat_zero(S)
       !call mat_zero(D)
       !
       call typedef_free_setting(ls%setting)
       call screen_free()
       call typedef_init_setting(ls%setting)
       call screen_init()
       ! Empirical dispersion correction in case of dft
       !CALL II_DFTDISP(LS%SETTING,DUMMY,1,1,0,LUPRI,1)

       !
       !   Setting DFT grid equal to zero in order to recalculate
       !   it at new geometry
       !
       If (ls%input%do_dft) then
          ls%input%dalton%DFT%griddone = 0
       Endif
       call typedef_set_default_setting(ls%setting,ls%input)

       if ((config%opt%cfg_start_guess == 'TRILEVEL')&
            &.or.(config%opt%cfg_start_guess == 'ATOMS')&
            &.or.config%decomp%cfg_gcbasis) then
          !     Not working properly for geometry-optimization
          !     config%diag%cfg_restart = .FALSE.
          IF(ls%input%BASIS%GCtransAlloc)THEN
             IF(ls%input%BASIS%GCtrans%natomtypes .NE. 0)THEN
                call free_basissetinfo(ls%input%BASIS%GCtrans)
             ENDIF
          ENDIF
          IF(ls%input%BASIS%VALENCE%natomtypes .NE. 0)THEN
             call free_basissetinfo(ls%input%BASIS%VALENCE)
          ENDIF
          if(config%decomp%cfg_gcbasis) call trilevel_basis(config%opt,ls)
       endif
       ls%setting%integraltransformGC = integraltransformGC
       !
       ! New nuclear repulsion
       PotNuc = 0E0_realk
       CALL II_get_nucpot(lupri,luerr,ls%setting,PotNuc)
       config%opt%potnuc = POTNUC
       ls%input%potnuc = POTNUC
       ! New energy
       Write(*,*)'CALLING FOR NEW ENERGY!'
       Call II_get_overlap(lupri,luerr,ls%setting,S)
       Call II_get_h1(lupri,luerr,ls%setting,H1)
       lsint_fock_data%ls => ls
       lsint_fock_data%H1 => H1
       lsint_fock_data%lupri = lupri
       lsint_fock_data%luerr = luerr

       ! We get initial density unless it's propagated when dynamics is done
       If (.NOT. config%dynamics%Start_propagation) then
           ! Everything else (e.g. optimization)
           Call get_initial_dens(H1,S,D,ls,config)
       Else
           config%decomp%S => S ! For decomp. to work for extr.guess
           ! For Fock matrix dynamics we diagonalize extrapolated Fock matrix
           Write(*,*)'FMD is done'
           If (config%dynamics%FockMD) then
              Call mat_init(CMO,D(1)%nrow,D(1)%ncol)
              allocate(eival(D(1)%nrow))
              Call mat_diag_f(F(1),S,eival,CMO)
              deallocate(eival)
              Call mat_density_from_orbs(CMO,D(1),config%decomp%nocc,config%decomp%nocca,config%decomp%noccb)
              Call mat_free(CMO)
           Endif
       Endif

       if (config%opt%cfg_incremental) call ks_init_incremental_fock(nbast)
       if (do_decomp) then
          call decomp_init(nbast,config%decomp)
          config%decomp%S => S
          call decomposition(config%decomp)
       else if (config%opt%cfg_start_guess == 'TRILEVEL') then
          call mat_free(config%decomp%lcv_CMO)
       endif

       if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then !FIXME: put this somewhere else!
          call mat_init(config%av%Fprev,nbast,nbast)
          call mat_init(config%av%Dprev,nbast,nbast)
       endif

       Call scfloop(H1,F,D,S,E,ls,config)
       if (config%opt%cfg_incremental) call ks_free_incremental_fock()
       ! Add a Grimme correction to the energy
       ! E = E + ls%setting%EDISP

       !
       ! We check whether it is a dec-calculation with local orbitals
       !lcm basis
       if (config%decomp%cfg_lcm) then
          ! get orbitals
          allocate(eival(nbast))
          call mat_diag_f(F(1),S,eival,C)
          deallocate(eival)

         ! localize orbitals
          call leastchange_lcm(config%decomp,C,config%decomp%nocc,ls)
       endif

       if (config%decomp%cfg_mlo) &
      & call optimloc(C,config%decomp%nocc,config%decomp%cfg_mlo_m,ls,config%davidOrbLoc)

       If (config%doDEC.AND.(.NOT.config%noDecEnergy)) then
          ! Get dec energy
          call get_total_mp2energy_from_inputs(ls,F(1),D(1),S,C,E(1),Eerr)
       elseif(config%doESGopt)then
          call get_excitation_energy(ls,config,F(1),D(1),S,ExcitE,&
               & config%decomp%cfg_rsp_nexcit)       
          Write(lupri,'(A,ES20.9)')'Ground state SCF Energy:',E(1)
          Write(lupri,'(A,ES20.9)')'Excitation Energy      :',ExcitE
          E(1) = E(1) + ExcitE
          Write(lupri,*)'==============================================='
          Write(lupri,'(A,ES20.9)')'Exicted state Energy   :',E(1)
          Write(lupri,*)'==============================================='
       Endif
       !
    !
  End subroutine Get_Energy
! Get energy: calculates molecular gradient in a general way
!> \currently used for energies at modified goemetries
!> \ covers SCF and DEC energies + dispersion for DFT
!> \author \latexonly V. Rybkin  \endlatexonly
!> \date 2012-27-02
  Subroutine Get_Gradient(E,Eerr,lupri,NAtoms,S,F,D,ls,config,C,Gradient)
    !
    ! Calls II_get_molecular_gradient  
    !
    Implicit none
    Integer :: lupri,NAtoms,i,j
    Real(realk) :: E     ! MP2 energy, auxiliary
    Type(Matrix), intent(inout) :: S  ! overlap matrix
    Type(Matrix), intent(inout) :: F,D   ! Fock and density matrix
    Type(Matrix), intent(inout) :: C       ! Orbitals
    Type(lsitem) :: ls
    Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk) :: Gradient(3,NAtoms)
    real(realk) :: Eerr  ! For DEC: Estimated intrinsic energy error

    Gradient = 0E0_realk
    Eerr     = 0E0_realk
    ! Calculate gradient

    ! Check whether it is a dec calculation
    If (DECinfo%doDEC) then
       ! Gradient from DEC (currently only MP2)
       Call get_mp2gradient_and_energy_from_inputs(ls,F,D,S,C,Natoms,gradient,E,Eerr)
    elseif(config%doESGopt)then
       call GET_EXCITED_STATE_GRADIENT(ls,config,F,D,S,Gradient,Natoms)
    else
       ! HF or DFT gradient
       Call II_get_molecular_gradient(Gradient,lupri,F,D,ls%setting,ls%input%do_dft,.TRUE.)
    Endif
    
  End subroutine Get_Gradient

!
End module Energy_and_deriv
