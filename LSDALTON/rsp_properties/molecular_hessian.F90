!> @file
!> Contains the Hessian module, with routines specific to the molecular Hessian.
!> \brief Contains Hessian specific routines.
MODULE molecular_hessian_mod
#ifdef MOD_UNRELEASED
  use precision ! realk
  use matrix_module, only: matrix
  use matrix_Operations, only: mat_mul, mat_free, mat_zero, mat_tr, &
       & mat_sqnorm2, mat_init
  use typedeftype, only: LSsetting,geoHessianConfig
  use memory_handling, only: mem_alloc,mem_dealloc
  use matrix_util, only: VerifyMatrices
  use lstiming, only: lstimer
  use integralinterfaceMod, only: ii_get_geoderivexchange, ii_get_overlap, &
       & ii_get_geoderivoverlap
#endif
#ifdef BUILD_GEN1INT_LSDALTON
  use gen1int_host
#endif

  private
  public :: dummy_subroutine_hessian,&
       & geohessian_set_default_config,&
       & get_molecular_hessian,&
       & get_first_geoderiv_overlap,&
       & get_first_geoderiv_refDmat

CONTAINS

  SUBROUTINE dummy_subroutine_hessian()
  END SUBROUTINE dummy_subroutine_hessian

#ifdef MOD_UNRELEASED
  !> \brief Set default settings for the geometrical Hessian calculation (default: run nothing)
  !> \author Patrick Merlot
  !> \date October,26 2012
  SUBROUTINE geohessian_set_default_config(geoHessianConf)
    !
    IMPLICIT NONE
    Type(geoHessianConfig),INTENT(INOUT)  :: geoHessianConf
    !
    geoHessianConf%do_geoHessian = .false.
    geoHessianConf%testContrib   = .false.
    geoHessianConf%DebugGen1Int  = .false.
    geoHessianConf%IntPrint      = 1
  END SUBROUTINE geohessian_set_default_config

  !> \brief Calculates the Hessian tensor
  !> \author \latexonly P. Merlot  \endlatexonly
  !> \date 2012-09-10
  !> \param Hessian The Hessian tensor
  !> \param F The Fock/Kohn-Sham matrix
  !> \param D The "reference" Density matrix D0
  !> \param setting Integral evalualtion settings
  !> \param lupri Default print unit
  !> \param luerr Default error print unit
  SUBROUTINE get_molecular_hessian(Hessian,Natoms,F,D,setting,config_hessian,lupri,luerr)
    !
    IMPLICIT NONE
    Real(realk),INTENT(INOUT)         :: Hessian(3*Natoms,3*Natoms)
    Type(LSSETTING),INTENT(INOUT)     :: setting
    Type(geoHessianConfig),INTENT(IN) :: config_hessian
    Integer,INTENT(IN)                :: Natoms,lupri,luerr
    Type(matrix),INTENT(IN)           :: F,D
    Type(matrix)           :: S
    Type(matrix),pointer   :: Sa(:),Ka(:),Da(:) !Sx,Sy,Sz for each atom
    Type(matrix),pointer   :: ha(:)
    Integer                :: i,nbast,iprint,ndmat
    Real(realk)            :: ts,te,sum
    !
    call lstimer('START ',ts,te,lupri)
    ndmat = 1
    nbast = D%nrow
    Hessian = 0E0_realk
    iprint = config_hessian%intprint

    ! calculate the overlap matrix S
    call mat_init(S,nbast,nbast)
    call II_get_overlap(lupri,luerr,setting,S)

    ! calculate the 1st geo. deriv. of the overlap matrix S^a
    call mem_alloc(Sa,3*Natoms)
     DO i=1,3*Natoms
        call mat_init(Sa(i),nbast,nbast)
     ENDDO
    call get_first_geoderiv_overlap(Sa,Natoms,setting,lupri,luerr)

    ! calculate the 1st geo. deriv. of the reference density matrix D_0^a
    call mem_alloc(Da,3*Natoms)
    DO i=1,3*Natoms
       call mat_init(Da(i),nbast,nbast)
    ENDDO
    call get_first_geoderiv_refDmat(Da,D,Sa,Natoms,setting,lupri,luerr,iprint)


    ! calculate the 1st geo. deriv. of the one-electron Hamiltonian
    call mem_alloc(ha,3*Natoms)
    DO i=1,3*Natoms
       call mat_init(ha(i),nbast,nbast)
    ENDDO
    ! call gen1int_host_get_first_geoderiv_h1(setting,ha,Natoms,lupri)

    ! calculate the 2nd geo. deriv. of the overlap matrix Sab
    ! call gen1int_host_get_second_geoderiv_overlap(ls_setting, Sxy, natoms,io_viewer)

    ! calculate the 2nd geo. deriv. of the one-electron Hamiltonian
    ! call gen1int_host_get_second_geoderiv_h1_expval(ls_setting, D,h1xy, natoms,io_viewer)

    ! calculate the 2nd geo. deriv. of the nuc-nuc repulsion

    ! generate the Right-Hand Side of the 1st order geo. response equation
    ! Solve the 1st order geo. resp. eq., and get the response vectors
    ! Calcualte 2nd deriv. of the density matrix and of D_{2n+1}

    ! CAN CHECK HESSIAN FOR ONE-ELECTRON SYSTEMS

    ! calculate the 1st geo. deriv. of the Coulomb contrib.
    ! calculate the 1st geo. deriv. of the Exchange contrib.
    call mem_alloc(Ka,3*Natoms)
    DO i=1,3*Natoms
       call mat_init(Ka(i),nbast,nbast)
    ENDDO
    call II_get_geoderivExchange(Ka,D,Natoms,setting,lupri,luerr)
    IF (IPRINT .GE. 3) THEN
       sum = 0.0E0_realk
       DO i=1,3*Natoms
          sum = sum + mat_sqnorm2(Ka(i))
       ENDDO
       WRITE(LUPRI,*) '   - Cumul. norm of Ka: ', sum
       WRITE(*,*)     '   - Cumul. norm of Ka: ', sum
    ENDIF

    ! calculate the 2nd geo. deriv. of the Coulomb contrib.
    ! calculate the 2nd geo. deriv. of the Exchange contrib.

    ! ! calculate the RHS of the response equations

    ! ! solve the response equations
    ! !> LINEQ == TRUE to solve with a RHS ( E(2)-w(I)S(2))X(I) - GD = 0 )
    ! DO i=1,3*Natoms
    !    ! Initialize solver parameters.
    !    call init_prop_molcfg(molcfg,S,Natoms,lupri,luerr,setting,config%decomp,config%solver)
    !    !!!   rsp_init(ntrial, nrhs, nsol, nomega, nstart)
    !    !call rsp_init(rsp_number_of_current_trial,rsp_number_of_rhs,&
    !    !        &rsp_number_of_sols,rsp_number_of_omegas,rsp_number_of_startvecs)
    !    call rsp_init(nexci_max,1,nexci_max,nexci_max,nstart)
    !    ! Calling solver. 
    !    LINEQ = .TRUE.
    !    nb_eq = 1
    !    ExEnergies = 0
    !    call rsp_solver(molcfg,D,S,F,LINEQ,nb_eq,RHS,ExEnergies,Xa)
    !    ! At this point Dx contain the excitation vectors, not the transition densities.
                     
    !     ENDDO


    ! free memory
     DO i=1,3*Natoms
        call mat_free(Sa(i))
        call mat_free(Ka(i))
        call mat_free(ha(i))
        call mat_free(Da(i))
     ENDDO

     call mem_dealloc(ha)
     call mem_dealloc(Ka)
     call mem_dealloc(Da)
     call mem_dealloc(Sa)
     call mat_free(S)
     call lstimer('geoHessian',ts,te,lupri)
  END SUBROUTINE get_molecular_hessian




!   !> \brief Calculates the Right Hand Side for 1st order HF density-matrix 
!   !> \brief response equation for geometrical perturbation
!   !> \latexonly
!   !>        $RHS = -\textit{P_A} \[ \textbf{h}^a \textbf{DS} + \textbf{G(D)DS}^a + \textbf{FDS}^a 
!   !>                + \textbf{K(D^a_{2n+1})}  \] $
!   !>  ref. http://dx.doi.org/10.1063.1.1415082
!   !> \endlatexonly
!   !> \author \latexonly P. Merlot  \endlatexonly
!   !> \date 2012-09-11
!   !> \param RHS The Right Hand Side of 1st order response equation
!   !> \param S The overlap matrix
!   !> \param S_a The 1st order geo. deriv. of the S
!   !> \param F The Fock/Kohn-Sham matrix
!   !> \param D The reference density matrix
!   !> \param Da The 1st geo. deriv. of the reference matrix D0
!   !> \param Jmat_a The 1st geo. deriv. of Jmat
!   !> \param Kmat_a The 1st geo. deriv. of Kmat
!   !> \param setting Integral evalualtion settings
!   !> \param lupri Default print unit
!   !> \param luerr Default error print unit
!   SUBROUTINE get_geom_first_order_RHS_HF(RHS,Natoms,S,S_a,F,D,Da,Jmat_a,Kmat_a,setting,lupri,luerr)
!     !
!     IMPLICIT NONE
!     Integer,INTENT(IN)            :: Natoms,lupri,luerr
!     TYPE(LSSETTING),intent(INOUT) :: setting
!     Type(matrix),INTENT(INOUT)    :: RHS(3*Natoms) ! derivative along x,y and z for each atom
!     Type(matrix),INTENT(IN)       :: S ! Overlap matrix
!     Type(matrix),INTENT(IN)       :: S_a(3*Natoms) ! 1st geo. deriv. of S
!     Type(matrix),intent(IN)       :: D ! Density matrix
!     Type(matrix),intent(IN)       :: F ! Fock matrix
!     Type(matrix),INTENT(IN)       :: Da(3*Natoms) ! 1st geo. deriv. of D
!     Type(matrix),INTENT(IN)       :: Jmat_a(3*Natoms) ! 1st geo. deriv. of Jmat
!     Type(matrix),INTENT(IN)       :: Kmat_a(3*Natoms) ! 1st geo. deriv. of Kmat
!     !
!     Integer          :: nbast
!     !
!     nbast = S%nrow

!   END  SUBROUTINE get_geom_first_order_RHS_HF


  !> \brief Calculates the first geometric derivative of the overlap matrix
  !> \author \latexonly P. Merlot  \endlatexonly
  !> \date 2012-09-11
  !> \param Sa The 3*Natoms first geometric derivative components of the overlap matrix
  !> \param setting Integral evalualtion settings
  !> \param lupri Default print unit
  !> \param luerr Unit for error printing
  SUBROUTINE get_first_geoderiv_overlap(Sa,Natoms,setting,lupri,luerr)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)            :: Natoms,lupri,luerr
    Type(LSSETTING),intent(INOUT) :: setting
    Type(matrix),INTENT(INOUT)    :: Sa(3*Natoms) ! derivative along x,y and z for each atom
    Real(realk)                   :: ts,te
    call lstimer('START ',ts,te,lupri)
    call II_get_geoderivOverlap(Sa,Natoms,setting,lupri,luerr)
    call lstimer('Sa_build',ts,te,lupri)
  END SUBROUTINE get_first_geoderiv_overlap


  !> \brief Calculates the first geometric derivative of the reference density matrix
  !> \author \latexonly P. Merlot  \endlatexonly
  !> \date 2012-09-11
  !> \param D The reference density matrix
  !> \param Da The 3*Natoms geo. deriv. components of the ref. density matrix: Da = -D Sa D
  !> \param Sa The 3*Natoms geometric derivative components of the overlap matrix
  !> \param setting Integral evalualtion settings
  !> \param lupri Default print unit
  !> \param luerr Unit for error printing
  !> \param iprint the printlevel, determining how much output should be generated
  SUBROUTINE get_first_geoderiv_refDmat(Da,D,Sa,Natoms,setting,lupri,luerr,iprint)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)            :: Natoms,lupri,luerr,iprint
    TYPE(LSSETTING),intent(INOUT) :: setting
    Type(matrix),INTENT(IN)       :: D
    Type(matrix),INTENT(IN)       :: Sa(3*Natoms)
    Type(matrix),INTENT(INOUT)    :: Da(3*Natoms) ! derivative along x,y and z for each atom
    !
    Type(matrix)                  :: temp
    Integer                       :: i,nbast
    Real(realk)                   :: sum
    Real(realk)                   :: ts,te
    !
    call lstimer('START ',ts,te,lupri)
    nbast = D%nrow
    call mat_init(temp,nbast,nbast)
    DO i=1,3*Natoms
       call mat_mul(D,Sa(i),'N','N',1E0_realk,0E0_realk,temp)
       call mat_mul(temp,D,'N','N',-1E0_realk,0E0_realk,Da(i)) ! Da = -D Sa D
    ENDDO
    IF (IPRINT .GE. 3) THEN
       sum = 0.0E0_realk
       DO i=1,3*Natoms
          sum = sum + mat_sqnorm2(Da(i))
       ENDDO
       WRITE(LUPRI,*) '   - Cumul. norm of Da: ', sum
       WRITE(*,*)     '   - Cumul. norm of Da: ', sum
    ENDIF
    call mat_free(temp)
    call lstimer('Da_build',ts,te,lupri)
  END SUBROUTINE get_first_geoderiv_refDmat


#endif
END MODULE molecular_hessian_mod
