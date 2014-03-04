!> @file 
!> Contains response driver.

MODULE lsdalton_rsp_mod
!> \brief Driver for stand-alone f90 linear scaling SCF.
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2008-10-26
!>
!> PROPERTIES SECTION
!>
use precision
use configurationType, only: configitem
use TYPEDEFTYPE,   only: LSSETTING,lsitem
use matrix_module, only: matrix
#ifdef VAR_RSP
use response_wrapper_module,only: get_dipole_moment, &
     & calculate_and_store_transition_density_matrices, &
     & oparesponse_driver, tparesponse_driver,&
     & esgresponse_driver, esdresponse_driver, &
     & mcdresponse_driver, free_transition_density_matrices, &
     & alpharesponse_driver, betaresponse_driver, &
     & gammaresponse_driver, dtparesponse_driver, &
     & nmrshieldresponse_driver, get_excitation_energies
use lstiming,      only: lstimer
use lsdalton_rsp_contribs,  only: rsp_oneave
use rspsolver,     only: rsp_molcfg, init_rsp_molcfg
use lsdalton_rsp_equations, only: rsp_eq_sol_empty
use rsp_util,      only: util_save_MOinfo,util_free_MOstuff
!use molecule_type, only: MOLECULE_PT,MOLECULEINFO
!use matrix_defop,  only: operator(*)
use memory_handling, only: mem_alloc,mem_dealloc
use integralinterfaceMod, only: ii_get_molecular_gradient
#endif

private
public :: lsdalton_response, get_excitation_energy, &
     & GET_EXCITED_STATE_GRADIENT, LS_rsp_eq_sol_empty
contains
  SUBROUTINE lsdalton_response(ls,config,F,D,S)
    implicit none
    TYPE(lsitem),target     :: ls
    type(configItem),target :: config
    TYPE(Matrix)            :: F,D,S
#ifdef VAR_RSP
    integer                 :: lu_pri,luerr,nbast
    logical                 :: dodft
    real(realk)             :: Tstart,Tend,t1,t2 !,ten,tstr,E,gradnrm
    real(realk)             :: DipoleMoment(3)
    type(rsp_molcfg)        :: molcfg
    ! Molecular gradient
    real(realk), pointer   :: Grad(:,:)

    lu_pri = config%lupri
    if (config%opt%calctype == config%opt%dftcalc) then
       dodft = .true.
    else
       dodft = .false.
    endif
    call CPU_TIME(tstart)
    call LSTIMER('START',t1,t2,LU_PRI)

    !create config struct to be passed to rsp_contribs / rsp_equations
    !    molcfg = rsp_molcfg(0E0_realk*S,ls%setting%MOLECULE(1)%p%Natoms, &
    !         & config%decomp%lupri,config%decomp%luerr, &
    !         & ls%setting,config%decomp,config%response%rspsolverinput)
    call init_rsp_molcfg(molcfg,S,ls%setting%MOLECULE(1)%p%Natoms, &
         & config%decomp%lupri,config%decomp%luerr, &
         & ls%setting,config%decomp,config%response%rspsolverinput)

    ! Kasper K, we ALWAYS calculate the permanent electric dipole for LSDALTON -
    ! also if no response properties have been requested.
    call Get_dipole_moment(molcfg,F,D,S,.true.,DipoleMoment)


    if(config%response%tasks%doResponse)Then
       !ajt This call generates Cmo/orbe in rsp_util, which rsp_solver needs
       if(config%response%RSPSOLVERinput%rsp_mo_precond) then
          call util_save_MOinfo(F,S,config%decomp%nocc) 
       endif

       ! Determine excited states and one-photon absorption
       if(config%response%tasks%doOPA)then
          ! Determine transition density matrices and corresponding excitation energies.
          ! The transition density matrices and corresponding excitation energies
          ! are stored in solved_eqs(:)%D and solved_eqs(:)%freq(1) inside rsp_equations.
          ! IMPORTANT!!! This must be done BEFORE any response calculations involving excited states,
          ! i.e. residues of response functions etc.
          ! Therefore, this call should be the first thing in lsdalton_response!
         call calculate_and_store_transition_density_matrices(molcfg,F,D,S)
          ! Calculate one-photon absorption strengths for the excited states determined above
          call OPAresponse_driver(molcfg,F,D,S)
       endif

       ! KK: Place all response properties involving excited states here:
       ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

       ! TPA
       if(config%response%tasks%doTPA)then
          call TPAresponse_driver(molcfg,F,D,S,config%response%tpainput)
       endif

       ! Excited state gradient
       if(config%response%tasks%doESG)then
          if(.NOT.config%doESGopt)then
             call ESGresponse_driver(molcfg,F,D,S,config%response%ESGinput)
          endif
          if(config%response%esginput%specific_states_in_input) &
               & call mem_dealloc(config%response%esginput%ExStates)
       endif

       ! Excited state dipole moment
       if(config%response%tasks%doESD)then
          call ESDresponse_driver(molcfg,F,D,S,config%response%ESDinput,DipoleMoment)
       endif

       ! Magnetic Circular Dichroism MCD
       if(config%response%tasks%doMCD)then
          call MCDresponse_driver(molcfg%lupri,molcfg%setting,molcfg%decomp,molcfg%solver,F,D,S,config%response%MCDinput)
       endif

       ! MCD, three-photon absorption etc.
       ! etc.

       ! Now all response properties involving excited states have been determined
       ! and we may delete the transition moment density matrices
       ! stored in solved_eqs(:)%D inside rsp_equations.
       ! IMPORTANT!!! This call has to be placed AFTER all response properties
       ! involving excited states!!!
       if(config%response%tasks%doOPA)then
          call free_transition_density_matrices(molcfg)
       endif

       ! Molecular gradient
       ! 29/10/11 Vladimir Rybkin: need Grad as input argument to 
       ! ii_get_molecular gradient elsewhere, therefore it's declared
       ! and allocated/deallocated here as well
       if(config%response%tasks%doGrad)then
          call mem_alloc(Grad,3,config%Molecule%NAtoms)
          call ii_get_molecular_gradient(Grad,lu_pri,F,D, &
               & ls%setting,dodft,.true.)
          call mem_dealloc(Grad)
       endif
    
       ! Polarizability
       if(config%response%tasks%doALPHA)then
          call ALPHAresponse_driver(molcfg,F,D,S,config%response%alphainput)
       endif

       ! 1st hyperpolarizability
       if(config%response%tasks%doBETA)then
          call BETAresponse_driver(molcfg,F,D,S,config%response%betainput,DipoleMoment)
       endif

       ! 2nd hyperpolarizability
       if(config%response%tasks%doGAMMA)then
          call GAMMAresponse_driver(molcfg,F,D,S,config%response%gammainput)
       endif

       ! Damped two-photon absorption
       if(config%response%tasks%doDTPA)then
          call DTPAresponse_driver(molcfg,F,D,S,config%response%dtpainput)
       endif

       if(config%response%tasks%doNMRshield) then
          call NMRshieldresponse_driver(molcfg,F,D,S)
       endif
       ! Free any allocated Cmo/orbe in rsp_util
       if(config%response%RSPSOLVERinput%rsp_mo_precond)then
          call util_free_MOstuff()
       endif
    endif

    call LSTIMER('LSDALTON RSP',t1,t2,LU_PRI,.TRUE.)
    call CPU_TIME(tend)
    WRITE(lu_pri,*) "*****************************************************"
    Write(lu_pri,*) "**     CPU-TIME USED IN LSDALTON RESPONSE: ",tend-tstart,"   **"
    WRITE(lu_pri,*) "*****************************************************"

    ! Clear saved solutions of response equations, stored in
    call rsp_eq_sol_empty()  !rsp_eq_sol in module rsp_equations

#endif
  END SUBROUTINE LSDALTON_RESPONSE

!    nexci = molcfg%decomp%cfg_rsp_nexcit
  SUBROUTINE get_excitation_energy(ls,config,F,D,S,ExcitE,nexcit)
    implicit none
    TYPE(lsitem),target     :: ls
    type(configItem),target :: config
    TYPE(Matrix)            :: F,D,S
    real(realk)             :: ExcitE
    integer                 :: nexcit
#ifdef VAR_RSP
    integer                 :: lupri,luerr,nbast
    type(rsp_molcfg)        :: molcfg
    real(realk),pointer     :: Excit(:)
    !integer
    integer :: i
    lupri = config%lupri
    call init_rsp_molcfg(molcfg,S,ls%setting%MOLECULE(1)%p%Natoms, &
         & config%decomp%lupri,config%decomp%luerr, &
         & ls%setting,config%decomp,config%response%rspsolverinput)
    call mem_alloc(Excit,molcfg%decomp%cfg_rsp_nexcit)
    if(config%response%tasks%doResponse)Then
       !ajt This call generates Cmo/orbe in rsp_util, which rsp_solver needs
       if(config%response%RSPSOLVERinput%rsp_mo_precond) then
          call util_save_MOinfo(F,S,config%decomp%nocc) 
       endif
       ! Determine excited states
       if(config%response%tasks%doOPA)then
          ! Determine transition density matrices and corresponding excitation energies.
          ! The transition density matrices and corresponding excitation energies
          ! are stored in solved_eqs(:)%D and solved_eqs(:)%freq(1) inside rsp_equations.
          ! IMPORTANT!!! This must be done BEFORE any response calculations involving excited states,
          ! i.e. residues of response functions etc.
          ! Therefore, this call should be the first thing in lsdalton_response!
         call calculate_and_store_transition_density_matrices(molcfg,F,D,S)
          ! Calculate one-photon absorption strengths for the excited states determined above         
         call get_excitation_energies(Excit,nexcit)
         WRITE(lupri,'(A,I5,A)')'The ',nexcit,' excited state vertical excitation energies'
         do I=1,nexcit
            WRITE(lupri,'(F16.8,A)')Excit(I),'a.u.'
         enddo
         WRITE(lupri,'(A,I4)')'We choose state number',config%response%ESGinput%ExStates(1)
         ExcitE = Excit(config%response%ESGinput%ExStates(1))
       endif
       ! Free any allocated Cmo/orbe in rsp_util
       if(config%response%RSPSOLVERinput%rsp_mo_precond)then
          call util_free_MOstuff()
       endif
    endif
    call mem_dealloc(Excit)
!    ! Clear saved solutions of response equations, stored in
!    call rsp_eq_sol_empty()  !rsp_eq_sol in module rsp_equations
!    we keep them for calculation of EXCITED_STATE_GRADIENT
#endif
  END SUBROUTINE GET_EXCITATION_ENERGY

  SUBROUTINE GET_EXCITED_STATE_GRADIENT(ls,config,F,D,S,Grad,Natoms)
    implicit none
    TYPE(lsitem),target     :: ls
    type(configItem),target :: config
    TYPE(Matrix)            :: F,D,S
    integer                 :: Natoms
    real(realk)             :: Grad(3,Natoms)
#ifdef VAR_RSP
    integer                 :: lupri,luerr,nbast
    logical                 :: dodft
    real(realk)             :: Tstart,Tend,t1,t2 !,ten,tstr,E,gradnrm
    real(realk)             :: DipoleMoment(3)
    type(rsp_molcfg)        :: molcfg
    ! Molecular gradient

    lupri = config%lupri
    call init_rsp_molcfg(molcfg,S,ls%setting%MOLECULE(1)%p%Natoms, &
         & config%decomp%lupri,config%decomp%luerr, &
         & ls%setting,config%decomp,config%response%rspsolverinput)
    ! The excited states should already be stored in rsp_eq_sol
    !ajt This call generates Cmo/orbe in rsp_util, which rsp_solver needs
    if(config%response%RSPSOLVERinput%rsp_mo_precond) then
       call util_save_MOinfo(F,S,config%decomp%nocc) 
    endif
    ! call Excited state gradient
    if(config%response%tasks%doESG)then
       call ESGresponse_driver(molcfg,F,D,S,config%response%ESGinput,GRAD)
    endif
    ! Free any allocated Cmo/orbe in rsp_util
    if(config%response%RSPSOLVERinput%rsp_mo_precond)then
       call util_free_MOstuff()
    endif
!    ! Clear saved solutions of response equations, stored in
    call rsp_eq_sol_empty()  !rsp_eq_sol in module rsp_equations
#endif
  END SUBROUTINE GET_EXCITED_STATE_GRADIENT

  SUBROUTINE LS_rsp_eq_sol_empty
#ifdef VAR_RSP
    ! Clear saved solutions of response equations, stored in
    call rsp_eq_sol_empty()  !rsp_eq_sol in module rsp_equations
#endif
  END SUBROUTINE LS_RSP_EQ_SOL_EMPTY

end MODULE lsdalton_rsp_mod
