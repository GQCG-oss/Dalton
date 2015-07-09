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
     & typedef_free_setting, copy_molecule
use memory_handling, only: mem_alloc,mem_dealloc
use scfloop_module, only: scfloop
use matrix_operations, only: mat_init, mat_free, mat_diag_f, mat_assign
use basis_type, only: free_basissetinfo
use basis_typetype, only: VALBasParam,GCTBasParam
use ks_settings, only: ks_init_incremental_fock, ks_free_incremental_fock
use decompMod, only: decomp_shutdown, decomp_init, decomposition
use initial_guess, only: get_initial_dens
use dec_typedef_module, only: DECinfo
use lsdalton_fock_module, only: lsint_fock_data
use matrix_operations_aux, only: mat_density_from_orbs
use matrix_util, only: save_fock_matrix_to_file
use integralinterfaceMod, only: II_get_molecular_gradient,&
     & II_get_nucpot,II_get_overlap,II_get_h1,II_precalc_ScreenMat,&
     & II_get_fock_mat
use lsdalton_rsp_mod,only: get_excitation_energy, GET_EXCITED_STATE_GRADIENT
use dec_main_mod
use ls_util, only: ls_print_gradient
use molecule_typetype, only: moleculeinfo
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



       Eerr   = 0E0_realk
       ExcitE = 0E0_realk    !Zeroing to initialize
       nbast  = D(1)%nrow
       do_decomp =.TRUE. !(config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
!           & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
!           & config%decomp%cfg_check_converged_solution .or. &
!           & config%decomp%cfg_rsp_nexcit > 0 .or. config%integral%locallink) 
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
       call screen_init()
       ls%lupri = lupri
       ls%luerr = luerr
       ls%optlevel = 3
       call typedef_init_setting(ls%setting)

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
          IF(ls%input%BASIS%WBASIS(GCTBasParam))THEN
             call free_basissetinfo(ls%input%BASIS%BINFO(GCTBasParam))
          ENDIF
          ls%input%BASIS%WBASIS(GCTBasParam) = .FALSE.
          IF(ls%input%BASIS%WBASIS(VALBasParam))THEN
             call free_basissetinfo(ls%input%BASIS%BINFO(VALBasParam))
          ENDIF
          ls%input%BASIS%WBASIS(VALBasParam) = .FALSE.
          if(config%decomp%cfg_gcbasis) call trilevel_basis(config%opt,ls)
       endif
       ls%setting%integraltransformGC = integraltransformGC
       !
       ! Precalculate screening matrices
       call II_precalc_ScreenMat(LUPRI,LUERR,ls%SETTING)
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
          ! Need to free memory allocated previously in LSDALTON driver
          call mat_free(config%av%Fprev)
          call mat_free(config%av%Dprev)
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
          call get_total_CCenergy_from_inputs(ls,F(1),D(1),C,E(1),Eerr)
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
    Integer :: lupri,NAtoms,i,j,nbast
    Real(realk) :: E     ! MP2 energy, auxiliary
    Type(Matrix), intent(inout) :: S  ! overlap matrix
    Type(Matrix), intent(inout) :: F,D   ! Fock and density matrix
    Type(Matrix), intent(inout) :: C       ! Orbitals
    Type(lsitem) :: ls
    Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk) :: Gradient(3,NAtoms)
    real(realk) :: Eerr  ! For DEC: Estimated intrinsic energy error
    real(realk) :: step  ! For numerical gradient
    Gradient = 0E0_realk
    Eerr     = 0E0_realk
    ! Calculate gradient
#ifdef MOD_UNRELEASED
    IF (.NOT.config%response%tasks%doNumGrad) THEN !Analytical gradient
#endif
      ! Check whether it is a dec calculation
      If (DECinfo%doDEC) then
         ! Gradient from DEC (currently only MP2)
         Call get_mp2gradient_and_energy_from_inputs(ls,F,D,C,Natoms,gradient,E,Eerr)
      elseif(config%doESGopt)then
         call GET_EXCITED_STATE_GRADIENT(ls,config,F,D,S,Gradient,Natoms)
      else
         ! HF or DFT gradient
         Call II_get_molecular_gradient(Gradient,lupri,F,D,ls%setting,ls%input%do_dft,.TRUE.)
      Endif
#ifdef MOD_UNRELEASED
    ELSE !Numerical gradient
      step = 1E-5_realk
      call get_num_grad(step,lupri,ls%luerr,ls,S,F,D,C,config,Gradient)
    ENDIF
#endif
    !
  End subroutine Get_Gradient


! Calculates the numerical geometrical gradient of the energy
!> \author: P. Merlot (almost copy of the work by S. Reine and K. Dankel)
!> \date 2013-05-16
subroutine get_num_grad(h,lupri,luerr,ls,S,F,D,C,config,numerical_gradient)
implicit none
!
real(realk), intent(in) :: h
integer, intent(in)     :: lupri,luerr
type(lsitem)            :: ls
Type(Matrix)            :: S,F,D,C ! overlap matrix, Fock and density matrix, Orbitals
type(configItem)        :: config
real(realk)             :: numerical_gradient(3,ls%INPUT%MOLECULE%nAtoms)
!
Type(Matrix)            :: H1
real(realk)             :: E(1),Emin,Eplus, Eerr
integer                 :: i, j, nbast, nAtoms
Type(Matrix)            :: Fmat(1),Dmat(1)

nbast=D%nrow
CALL mat_init(H1,nbast,nbast)
nAtoms = ls%INPUT%MOLECULE%nAtoms

call mat_init(Dmat(1),nbast,nbast)
call mat_assign(Dmat(1),D)
call mat_init(Fmat(1),nbast,nbast)
call mat_assign(Fmat(1),F)


do i=1,ls%INPUT%MOLECULE%nAtoms
   do j=1, 3
      write (*,*) "atom index:",i,"  coord:",j
      write (lupri,*) "atom index:",i,"  coord:",j
      ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)-h 
      CALL get_energy(E,Eerr,config,H1,Fmat,Dmat,S,ls,C,nAtoms,lupri,luerr)
      Emin=E(1)
      
      ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)+(2*h)
      CALL get_energy(E,Eerr,config,H1,Fmat,Dmat,S,ls,C,nAtoms,lupri,luerr)
      Eplus=E(1)

      ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)-h
      
      numerical_gradient(j,i)=(Eplus-Emin)/(2*h)
    enddo
enddo
CALL LS_PRINT_GRADIENT(lupri,ls%setting%molecule(1)%p,numerical_gradient,nAtoms,'NUM_GRAD')

CALL mat_free(H1)
CALL mat_free(Dmat(1))
CALL mat_free(Fmat(1))
end subroutine get_num_grad

End module Energy_and_deriv
