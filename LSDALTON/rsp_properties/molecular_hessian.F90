!> @file
!> Contains the Hessian module, with routines specific to the molecular Hessian.
!> \brief Contains Hessian specific routines.
MODULE molecular_hessian_mod
    use precision ! realk
    use matrix_module,          only: matrix, matrixp
    use matrix_Operations,      only: mat_mul, mat_daxpy, &
                                    & mat_free, mat_scal, &
                                    & mat_zero, mat_tr, &
                                    & mat_sqnorm2, mat_init, &  
                                    & mat_set_from_full, &
                                    & mat_assign
    use typedeftype,            only: LSsetting, geoHessianConfig
    use configurationType,      only: ConfigItem
    use memory_handling,        only: mem_alloc, mem_dealloc
    use matrix_util,            only: VerifyMatrices
    use lstiming,               only: lstimer
    use integralinterfaceMod,   only: II_get_overlap, &
                                    & II_get_geoderivoverlap, &
                                    & II_get_geoderivKinetic, &
                                    & II_get_geoderivnucel, &
                                    & II_get_geoderivCoulomb, &
                                    & II_get_geoderivexchange, &
                                    & II_get_coulomb_mat, &
                                    & II_get_exchange_mat
    use RSPsolver,              only: rsp_molcfg,&
                                    & init_rsp_molcfg,&
                                    & rsp_init,&
                                    & rsp_solver
    use RSP_util,               only: util_save_MOinfo,&
                                    & util_free_MOstuff
#ifdef BUILD_GEN1INT_LSDALTON
  use gen1int_host
#endif

  private   ::  get_first_order_rsp_vectors

  public    ::  dummy_subroutine_hessian,&
                & geohessian_set_default_config,&
                & get_molecular_hessian, &
                & get_first_geoderiv_overlap,&
                & get_first_geoderiv_refDmat,&
                & get_first_geoderiv_H1_mat, &
                & get_first_geoderiv_Coulomb_mat, &
                & get_first_geoderiv_exchange_mat, &
                & get_first_geoderiv_twoElectron_mat, &
                & get_first_geoderiv_Fock_mat,&
                & get_geom_first_order_RHS_HF

CONTAINS

  !> \brief Set default settings for the geometrical Hessian calculation 
  !>        (default: do not compute the geometrical molecular Hessian)
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
  !> \param config Info/Settings/Data for entire calculation
  !>               (defaults or read from input file)
  !> \param lupri Default print unit
  !> \param luerr Default error print unit
  SUBROUTINE get_molecular_hessian(Hessian,Natoms,F,D,setting,config,lupri,luerr)
    !
    IMPLICIT NONE
    Real(realk),INTENT(INOUT)         :: Hessian(3*Natoms,3*Natoms)
    Type(LSSETTING),INTENT(INOUT)     :: setting
    Type(ConfigItem),INTENT(IN)       :: config
    Integer,INTENT(IN)                :: Natoms,lupri,luerr
    Type(matrix),INTENT(IN)           :: F
    Type(matrix),INTENT(IN)           :: D
    Type(matrix)           :: S
    Type(matrix),pointer   :: Sa(:),Da(:),ha(:) !Sx,Sy,Sz for each atom
    Type(matrix),pointer   :: Ga(:),GDa(:),Fa(:),RHS_HF(:)
    Type(matrix),pointer   :: Xa(:)
    Integer                :: i,nbast,iprint,ndmat
    Real(realk)            :: ts,te
    !
    call lstimer('START ',ts,te,lupri)
    ndmat = 1
    nbast = D%nrow
    Hessian = 0E0_realk
    iprint = config%geoHessian%intprint

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

    ! calculate the 1st geo. deriv. of the one-electron Hamiltonian h^a
    call mem_alloc(ha,3*Natoms)
     DO i=1,3*Natoms
        call mat_init(ha(i),nbast,nbast)
    ENDDO
    call get_first_geoderiv_H1_mat(ha,Natoms,setting,lupri,luerr,iprint)

    ! calculate 1st geo. deriv. of the Coulomb and Exchange contrib.: G^a(D)=J^a(D)-K^a(D)
    call mem_alloc(Ga,3*Natoms)
    DO i=1,3*Natoms
        call mat_init(Ga(i),nbast,nbast)
    ENDDO
    call get_first_geoderiv_twoElectron_mat(Ga,D,Natoms,setting,lupri,luerr)

    ! calculate G(D^a)=J(D^a)-K(D^a)
    call mem_alloc(GDa,3*Natoms)
    DO i=1,3*Natoms
        call mat_init(GDa(i),nbast,nbast)
    ENDDO
    call get_twoElectron_mat_of_first_geoderiv_refDmat(GDa,Da,Natoms,setting,lupri,luerr,iprint)

    ! calculate the 1st geo. deriv. of the Fock matrix F^a
    call mem_alloc(Fa,3*Natoms)
    DO i=1,3*Natoms
       call mat_init(Fa(i),nbast,nbast)
    ENDDO
    call get_first_geoderiv_Fock_mat(Fa,D,ha,Ga,GDa,Natoms,setting,lupri,luerr)

    ! generate the Right-Hand Side of the 1st order geo. response equation
    call mem_alloc(RHS_HF,3*Natoms)
    DO i=1,3*Natoms
       call mat_init(RHS_HF(i),nbast,nbast)
    ENDDO
    call get_geom_first_order_RHS_HF(RHS_HF,Natoms,S,D,F,Sa,Da,Fa,&
                                    & setting,lupri,luerr)

    ! Solve the 1st order geo. resp. eq., and get the response vectors
    !> LINEQ == TRUE to solve with a RHS ( E(2)-w(I)S(2))X(I) - GD = 0 )
    call mem_alloc(Xa,3*Natoms)
     DO i=1,3*Natoms
        call mat_init(Xa(i),nbast,nbast)
        call mat_zero(Xa(i))
     ENDDO
    call get_first_order_rsp_vectors(Xa,S,D,F,RHS_HF,Natoms,setting,config,lupri,luerr,iprint)

    ! Calcualte 2nd deriv. of the density matrix and of D_{2n+1} 

    ! calculate the 2nd geo. deriv. of the overlap matrix Sab
    ! call gen1int_host_get_second_geoderiv_overlap(ls_setting, Sxy, natoms,io_viewer)

    ! calculate the 2nd geo. deriv. of the one-electron Hamiltonian
    ! call gen1int_host_get_second_geoderiv_h1_expval(ls_setting, D,h1xy, natoms,io_viewer)

    ! calculate the 2nd geo. deriv. of the nuc-nuc repulsion

    ! calculate the 2nd geo. deriv. of the Coulomb contrib.
    ! calculate the 2nd geo. deriv. of the Exchange contrib.


    ! free memory
     DO i=1,3*Natoms
        call mat_free(Sa(i))
        call mat_free(Da(i))
        call mat_free(ha(i))
        call mat_free(Ga(i))
        call mat_free(GDa(i))
        call mat_free(Fa(i))
        call mat_free(Xa(i))
        call mat_free(RHS_HF(i))
     ENDDO
     call mem_dealloc(Sa)
     call mem_dealloc(Da)
     call mem_dealloc(ha)
     call mem_dealloc(Ga)
     call mem_dealloc(GDa)
     call mem_dealloc(Fa)
     call mem_dealloc(Xa)
     call mem_dealloc(RHS_HF)
     call mat_free(S)
     call lstimer('geoHess_build',ts,te,lupri)
  END SUBROUTINE get_molecular_hessian




   !> \brief Calculates the Right Hand Side for 1st order HF density-matrix 
   !> \brief response equation for geometrical perturbation
   !> \latexonly
   !>        $RHS = -\textit{P_A} \[ \textbf{h}^a \textbf{DS} + \textbf{G(D)DS}^a + \textbf{FDS}^a 
   !>                + \textbf{K(D^a_{2n+1})}  \] $
   !>  ref. http://dx.doi.org/10.1063.1.1415082
   !> \endlatexonly
   !> \author P. Merlot
   !> \date 2013-07-11
   !> \param RHS The Right Hand Side of 1st order response equation
   !> \param S The overlap matrix
   !> \param D The reference density matrix D0
   !> \param F The Fock/Kohn-Sham matrix
   !> \param Sa The 1st order geo. deriv. of the S
   !> \param Da The 1st geo. deriv. of the reference matrix D0
   !> \param Fa The 1st order geo. deriv. of the Fock/Kohn-Sham matrix
   !> \param setting Integral evalualtion settings
   !> \param lupri Default print unit
   !> \param luerr Default error print unit
   SUBROUTINE get_geom_first_order_RHS_HF(RHS,Natoms,S,D,F,Sa,Da,Fa,setting,lupri,luerr)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)            :: Natoms,lupri,luerr
    Type(LSSETTING),intent(INOUT) :: setting
    Type(matrix),INTENT(INOUT)    :: RHS(3*Natoms) ! derivative along x,y and z for each atom
    Type(matrix),INTENT(IN)       :: S ! Overlap matrix
    Type(matrix),intent(IN)       :: D ! Density matrix
    Type(matrix),intent(IN)       :: F ! Fock matrix
    Type(matrix),INTENT(IN)       :: Sa(3*Natoms) ! 1st geo. deriv. of S
    Type(matrix),INTENT(IN)       :: Da(3*Natoms) ! 1st geo. deriv. of D
    Type(matrix),INTENT(IN)       :: Fa(3*Natoms) ! 1st geo. deriv. of F
    !
    Integer          :: nbast,i
    Type(matrix)     :: temp
    Real(realk)      :: ts,te
    !
    nbast = S%nrow
    ! RHS = 0.5 { SDFa - FaDS + SaDF + SDaF - FDaS  - FDSa }
    call lstimer('START ',ts,te,lupri)
    nbast = D%nrow
    call mat_init(temp,nbast,nbast)
    call mat_mul(S,D,'N','N',1E0_realk,0E0_realk,temp) ! temp = SD
    DO i=1,3*Natoms
        call mat_mul(temp,Fa(i),'N','N', 1E0_realk,0E0_realk,RHS(i)) ! RHS = SDFa
        call mat_mul(Fa(i),temp,'N','T',-1E0_realk,1E0_realk,RHS(i)) ! RHS = SDFa - FaDS
    ENDDO
    call mat_mul(D,F,'N','N',1E0_realk,0E0_realk,temp) ! temp = DF
    DO i=1,3*Natoms
        call mat_mul(Sa(i),temp,'N','N', 1E0_realk,1E0_realk,RHS(i)) ! RHS += SaDF
        call mat_mul(temp,Fa(i),'T','N',-1E0_realk,1E0_realk,RHS(i)) ! RHS -= FDSa 
    ENDDO
    DO i=1,3*Natoms
        call mat_mul(S,Da(i),'N','N', 1E0_realk,0E0_realk,temp)   ! temp = SDa(i)
        call mat_mul(temp,F ,'N','N', 1E0_realk,1E0_realk,RHS(i)) ! RHS += FDSa 
        call mat_mul(F,Da(i),'N','N',-1E0_realk,0E0_realk,temp)   ! temp = -FDa(i)
        call mat_mul(temp,S ,'N','N', 1E0_realk,1E0_realk,RHS(i)) ! RHS -= FDaS 
    ENDDO

    DO i=1,3*Natoms
        ! RHS = 0.5 { ... }
        !call mat_scal(0.5E0_realk,RHS(i))
        call mat_scal(0.125E0_realk,RHS(i))
    ENDDO

    call mat_free(temp)
    call lstimer('RHS_HF_build',ts,te,lupri)
    !     
   END  SUBROUTINE get_geom_first_order_RHS_HF


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

  !> \brief Calculates the first geometric derivative of the one-electron Hamiltonian
  !> Nuclear-electron attraction + kinetic energy gradient
  !> \author P. Merlot
  !> \date 2012-09-11
  !> \param ha The 3*Natoms first geometric derivative components of the one-electron Hamiltonian
  !> \param Natoms Nb. of atoms in the molecule
  !> \param Dmat The density matrix
  !> \param ndmat The number of density matrices
  !> \param setting Integral evalualtion settings
  !> \param lupri Default print unit
  !> \param luerr Unit for error printing
  !> \param iprint the printlevel, determining how much output should be generated
  SUBROUTINE get_first_geoderiv_H1_mat(ha,Natoms,setting,lupri,luerr,iprint)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)              :: Natoms,lupri,luerr,iprint
    Type(LSSETTING),intent(INOUT)   :: setting
    Type(matrix),INTENT(INOUT)      :: ha(3*Natoms) ! derivative along x,y and z for each atom
    !
    Real(realk)                     :: ts,te,sum
    Integer                         :: i, nbast
    Type(matrix), pointer           :: Va(:) 
    !
    call lstimer('START ',ts,te,lupri)
    nbast = ha(1)%nrow
    call mem_alloc(Va,3*Natoms)
     DO i=1,3*Natoms
        call mat_init(Va(i),nbast,nbast)
     ENDDO
    call II_get_geoderivKinetic(ha,Natoms,setting,lupri,luerr)
    call II_get_geoderivnucel(Va,Natoms,setting,lupri,luerr)
    DO i=1,3*Natoms
        call mat_daxpy(1.0E0_realk,Va(i),ha(i))
        call mat_free(Va(i))
     ENDDO
    IF (iprint .GE. 3) THEN
       sum = 0.0E0_realk
       DO i=1,3*Natoms
          sum = sum + mat_sqnorm2(ha(i))
       ENDDO
       WRITE(LUPRI,*) '   - Cumul. norm of H1a: ', sum
       WRITE(*,*)     '   - Cumul. norm of H1a: ', sum
    ENDIF
    call mem_dealloc(Va)
    call lstimer('ha_build',ts,te,lupri)
  END SUBROUTINE get_first_geoderiv_H1_mat

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


  !> \brief Calculates the first geometric derivative of the Coulomb matrix
  !> \author P. Merlot
  !> \date 2013-07-10
  !> \param D The reference density matrix
  !> \param Ja The 3*Natoms geo. deriv. components of the Coulomb matrix: Ja
  !> \param setting Integral evalualtion settings
  !> \param lupri Default print unit
  !> \param luerr Unit for error printing
  !> \param iprint the printlevel, determining how much output should be generated
  SUBROUTINE get_first_geoderiv_Coulomb_mat(Ja,D,Natoms,setting,lupri,luerr,iprint)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)            :: Natoms,lupri,luerr,iprint
    TYPE(LSSETTING),intent(INOUT) :: setting
    Type(matrix),INTENT(IN)       :: D
    Type(matrix),INTENT(INOUT)    :: Ja(3*Natoms) ! derivative along x,y and z for each atom
    !
    Integer                 :: i,nbast
    Real(realk)             :: ts,te,sum
    !
    call lstimer('START ',ts,te,lupri)
    nbast = D%nrow
    call II_get_geoderivCoulomb(Ja,D,Natoms,setting,lupri,luerr)
    !call mat_scal(1E0_realk,Ja)
    IF (iprint .GE. 3) THEN
       sum = 0.0E0_realk
       DO i=1,3*Natoms
          sum = sum + mat_sqnorm2(Ja(i))
       ENDDO
       WRITE(LUPRI,*) '   - Cumul. norm of Ja: ', sum
       WRITE(*,*)     '   - Cumul. norm of Ja: ', sum
    ENDIF
    call lstimer('Ja_build',ts,te,lupri)
  END SUBROUTINE get_first_geoderiv_Coulomb_mat


  !> \brief Calculates the first geometric derivative of the exchange matrix
  !> \author P. Merlot
  !> \date 2013-07-10
  !> \param D The reference density matrix
  !> \param Ja The 3*Natoms geo. deriv. components of the exchange matrix: Ka
  !> \param setting Integral evalualtion settings
  !> \param lupri Default print unit
  !> \param luerr Unit for error printing
  !> \param iprint the printlevel, determining how much output should be generated
  SUBROUTINE get_first_geoderiv_exchange_mat(Ka,D,Natoms,setting,lupri,luerr,iprint)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)            :: Natoms,lupri,luerr,iprint
    TYPE(LSSETTING),intent(INOUT) :: setting
    Type(matrix),INTENT(IN)       :: D
    Type(matrix),INTENT(INOUT)    :: Ka(3*Natoms) ! derivative along x,y and z for each atom
    !
    Integer                 :: i,nbast
    Real(realk)             :: ts,te,sum
    !
    call lstimer('START ',ts,te,lupri)
    nbast = D%nrow
    call II_get_geoderivExchange(Ka,D,Natoms,setting,lupri,luerr)
    !call mat_scal(1E0_realk,Ka)
    IF (iprint .GE. 3) THEN
       sum = 0.0E0_realk
       DO i=1,3*Natoms
          sum = sum + mat_sqnorm2(Ka(i))
       ENDDO
       WRITE(LUPRI,*) '   - Cumul. norm of Ka: ', sum
       WRITE(*,*)     '   - Cumul. norm of Ka: ', sum
    ENDIF
    call lstimer('Ka_build',ts,te,lupri)
  END SUBROUTINE get_first_geoderiv_exchange_mat

  !> \brief Calculates the first geometric derivative of the 2-electron matrix
  !> \author P. Merlot
  !> \date 2013-07-10
  !> \param Ga The 3*Natoms first geometric derivative components: G^a(D)
  !> \param D The reference density matrix
  !> \param Natoms Nb. of atoms in the molecule
  !> \param Dmat The density matrix
  !> \param ndmat The number of density matrices
  !> \param setting Integral evalualtion settings
  !> \param lupri Default print unit
  !> \param luerr Unit for error printing
  SUBROUTINE get_first_geoderiv_twoElectron_mat(Ga,D,Natoms,setting,lupri,luerr)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)              :: Natoms,lupri,luerr
    Type(matrix),INTENT(IN)         :: D
    Type(LSSETTING),intent(INOUT)   :: setting
    Type(matrix),INTENT(INOUT)      :: Ga(3*Natoms) ! derivative along x,y and z for each atom
    !
    Type(matrix),pointer            :: Ja(:),Ka(:)
    Real(realk)                     :: ts,te 
    Integer                         :: i, nbast
    !
    call lstimer('START ',ts,te,lupri)
    nbast = Ga(1)%nrow
    call mem_alloc(Ja,3*Natoms)
    call mem_alloc(Ka,3*Natoms)
    DO i=1,3*Natoms
        call mat_zero(Ga(i))
        call mat_init(Ja(i),nbast,nbast)
        call mat_zero(Ja(i))
        call mat_init(Ka(i),nbast,nbast)
        call mat_zero(Ka(i))
    ENDDO

    ! calculate 1st geo. deriv. of the Coulomb and Exchange contrib.: J^a, K^a
    call get_first_geoderiv_Coulomb_mat( Ja,D,Natoms,setting,lupri,luerr,0)
    call get_first_geoderiv_exchange_mat(Ka,D,Natoms,setting,lupri,luerr,0)

    DO i=1,3*Natoms
        call mat_daxpy( 1.0E0_realk,Ja(i),Ga(i))
        call mat_daxpy(-1.0E0_realk,Ka(i),Ga(i))
        call mat_free(Ja(i))
        call mat_free(Ka(i))
    ENDDO
    call mem_dealloc(Ja)
    call mem_dealloc(Ka)
    call lstimer('Ga_build',ts,te,lupri)
  END SUBROUTINE get_first_geoderiv_twoElectron_mat

  !> \brief Calculates the first geometric derivative of the 2-electron matrix
  !> \author P. Merlot
  !> \date 2013-07-10
  !> \param GDa defined as G(D^a) = J(D^a) - K(D^a)
  !> \param Da The 3*Natoms geo. deriv. components of the ref. density matrix
  !> \param Natoms Nb. of atoms in the molecule
  !> \param Dmat The density matrix
  !> \param ndmat The number of density matrices
  !> \param setting Integral evalualtion settings
  !> \param lupri Default print unit
  !> \param luerr Unit for error printing
  !> \param iprint the printlevel, determining how much output should be generated
  SUBROUTINE get_twoElectron_mat_of_first_geoderiv_refDmat(GDa,Da,Natoms,setting,lupri,luerr,iprint)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)              :: Natoms,lupri,luerr,iprint
    Type(matrix),INTENT(IN)         :: Da(3*Natoms)
    Type(LSSETTING),intent(INOUT)   :: setting
    Type(matrix),INTENT(INOUT)      :: GDa(3*Natoms) ! derivative along x,y and z for each atom
    !
    Type(matrix), pointer           :: tempK(:)
    Real(realk)                     :: ts,te,sum
    Integer                         :: i, nbast, ndmat
    LOGICAL                         :: Dsym
    !
    call lstimer('START ',ts,te,lupri)
    nbast = Da(1)%nrow
    ndmat = 3*Natoms
    call mem_alloc(tempK,3*Natoms)
    DO i=1,3*Natoms
        call mat_zero(GDa(i))
        call mat_init(tempK(i),nbast,nbast)
        call mat_zero(tempK(i))
    ENDDO
    call II_get_coulomb_mat( lupri,luerr,setting,Da,GDa,ndmat)
    IF (iprint .GE. 3) THEN
       sum = 0.0E0_realk
       DO i=1,3*Natoms
          sum = sum + mat_sqnorm2(GDa(i))
       ENDDO
       WRITE(LUPRI,*) '   - Cumul. norm of J(Da): ', sum
       WRITE(*,*)     '   - Cumul. norm of J(Da): ', sum
    ENDIF
    Dsym = .FALSE.
    call II_get_exchange_mat(lupri,luerr,setting,Da,ndmat,Dsym,tempK)
    IF (iprint .GE. 3) THEN
       sum = 0.0E0_realk
       DO i=1,3*Natoms
          sum = sum + mat_sqnorm2(tempK(i))
       ENDDO
       WRITE(LUPRI,*) '   - Cumul. norm of K(Da): ', sum
       WRITE(*,*)     '   - Cumul. norm of K(Da): ', sum
    ENDIF
    DO i=1,3*Natoms
        call mat_daxpy(-1.0E0_realk,tempK(i),GDa(i))
        call mat_free(tempK(i))
    ENDDO
    call mem_dealloc(tempK)
    call lstimer('GDa_build',ts,te,lupri)
  END SUBROUTINE get_twoElectron_mat_of_first_geoderiv_refDmat

  !> \brief Calculates the first geometric derivative of the Fock matrix
  !> \author P. Merlot
  !> \date 2013-07-10
  !> \param Fa The 3*Natoms first geometric derivative components of F
  !> \param D The reference density matrix
  !> \param ha The 3*Natoms first geometric derivative components of H1
  !> \param Ga The 3*Natoms first geometric derivative components: G^a(D)
  !> \param Natoms Nb. of atoms in the molecule
  !> \param setting Integral evalualtion settings
  !> \param lupri Default print unit
  !> \param luerr Unit for error printing
  SUBROUTINE get_first_geoderiv_Fock_mat(Fa,D,ha,Ga,GDa,Natoms,setting,lupri,luerr)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)              :: Natoms,lupri,luerr
    Type(matrix),INTENT(IN)         :: D
    Type(matrix),INTENT(IN)         :: ha(3*Natoms),Ga(3*Natoms),GDa(3*Natoms)
    Type(LSSETTING),INTENT(INOUT)   :: setting
    Type(matrix),INTENT(INOUT)      :: Fa(3*Natoms) ! derivative along x,y and z for each atom
!
    Real(realk)                     :: ts,te 
    Integer                         :: i, nbast
    !
    call lstimer('START ',ts,te,lupri)
    DO i=1,3*Natoms
        call mat_daxpy(1.0E0_realk,ha(i), Fa(i))
        call mat_daxpy(1.0E0_realk,Ga(i), Fa(i))
        call mat_daxpy(1.0E0_realk,GDa(i),Fa(i))
    ENDDO
    call lstimer('Fa_build',ts,te,lupri)
  END SUBROUTINE get_first_geoderiv_Fock_mat


  !> \brief Solve the 1st-order resp. eq. to get the response vectors Xa
  !> \author P. Merlot
  !> \date 2013-07-11
  !> \param S The overlap matrix
  !> \param D The reference density matrix
  !> \param F The Fock/Kohn-Sham matrix
  !> \param RHS The Right Hand Side of 1st order response equation
  !> \param Natoms Nb. of atoms in the molecule
  !> \param setting Integral evalualtion settings
  !> \param config Info/Settings/Data for entire calculation
  !>               (defaults or read from input file)
  !> \param lupri Default print unit
  !> \param luerr Unit for error printing
  !> \param iprint the printlevel, determining how much output should be generated
  SUBROUTINE get_first_order_rsp_vectors(Xa,S,D,F,RHS,Natoms,setting,config,lupri,luerr,iprint)
    !
    IMPLICIT NONE
    Integer,INTENT(IN)                  :: Natoms,lupri,luerr,iprint
    Type(ConfigItem),INTENT(IN)         :: config
    Type(matrix),INTENT(IN)             :: S,D,F
    Type(matrix),INTENT(INOUT)          :: RHS(3*Natoms)
    Type(LSSETTING),INTENT(INOUT)       :: setting
    Type(matrix),INTENT(INOUT)          :: Xa(3*Natoms) ! derivative along x,y and z for each atom
    !
    Real(realk)                     :: ts,te, sum
    Integer                         :: i, nbast, nb_eq
    type(rsp_molcfg)                :: molcfg
    logical                         :: LINEQ_x
    Real(realk)                     :: laser_freq(1)
    Type(matrix)                    :: oneRHS(1),oneXa(1)
    !
    call lstimer('START ',ts,te,lupri)
    nbast = D%nrow
    call mat_init(oneRHS(1),nbast,nbast)
    call mat_init(oneXa(1),nbast,nbast)

    ! Initialize solver parameters.
    call init_rsp_molcfg(molcfg, S, Natoms,&
                        & lupri, luerr, setting,&
                        & config%decomp,config%response%rspsolverinput)

    !> ntrial: (Max.) number of trial vectors in a given iteration
    !> nrhs:    Number of right-hand sides. Only relevant for linear equations (always 1 for eigenvalue problem)
    !> nsol:    Number of solution (output) vectors
    !> nomega:  If LINEQ, number of laser freq.s (input). Otherwise number of excitation energies (output) 
    !> Number of start vectors. Only relevant for eigenvalue problem
    !!!  rsp_init(ntrial, nrhs, nsol, nomega, nstart)
    call rsp_init(1,      1,    1,    1,      0)

    call util_save_MOinfo(F,S,config%decomp%nocc) !nocc: Number of occupied orbitals (if restricted)

    ! Calling the repsonse solver
    LINEQ_x = .TRUE. ! solving linear system, not eigenvalue problem
    nb_eq = 1        ! one response vector at a time
    laser_freq(1) = 0.0E0_realk
    DO i=1,3*Natoms
        write(*,*) "Solving response vector #",i, "out of ", 3*Natoms
        call mat_assign(oneRHS(1), RHS(i))
        WRITE(*,*)     'norm of RHS(i): ', mat_sqnorm2(RHS(i))
        call mat_zero(oneXa(1))
!        call rsp_solver(molcfg, D, S, F, &
!                        & LINEQ_x, nb_eq, oneRHS, laser_freq, oneXa)
        call rsp_solver(molcfg, D, S, F, &
                        & LINEQ_x, nb_eq, RHS(i:i), laser_freq(1:1), oneXa(1))
        call mat_assign(Xa(i), oneXa(1))
        WRITE(*,*)     'norm of Xa(i): ', mat_sqnorm2(oneXa(1))
    ENDDO
    call util_free_MOstuff()
    IF (iprint .GE. 3) THEN
       sum = 0.0E0_realk
       DO i=1,3*Natoms
          sum = sum + mat_sqnorm2(Xa(i))
       ENDDO
       WRITE(LUPRI,*) '   - Cumul. norm of Xa: ', sum
       WRITE(*,*)     '   - Cumul. norm of Xa: ', sum
    ENDIF
    call mat_free(oneRHS(1))
    call mat_free(oneXa(1))
    call lstimer('Xa_build',ts,te,lupri)
  END SUBROUTINE get_first_order_rsp_vectors

END MODULE molecular_hessian_mod

