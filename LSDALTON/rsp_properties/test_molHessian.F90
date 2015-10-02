!> @file
!> Contains test subroutines for the Hessian module.
!> \brief Contains Hessian specific test routines.
MODULE test_molecular_hessian_mod
  use precision ! realk
  use matrix_module,          only: matrix, matrixp
  use matrix_Operations,      only: mat_init, mat_mul, mat_free,&
                                  & mat_zero, mat_tr, &
                                  & mat_trab
  use typedeftype,            only: LSsetting,geoHessianConfig
!  use typedef
  use memory_handling,        only: mem_dealloc, mem_alloc
  use matrix_util,            only: VerifyMatrices
  use lstiming,               only: lstimer
  use molecular_hessian_mod,  only: get_first_geoderiv_overlap, &
                                  & get_first_geoderiv_H1_mat, &
                                  & get_first_geoderiv_refdmat, &
                                  & get_first_geoderiv_Coulomb_mat, &
                                  & get_first_geoderiv_exchange_mat
  use integralinterfaceMOD,   only: II_get_J_gradient, &
                                  & II_get_K_gradient, &
                                  & II_get_reorthoNormalization
#ifdef BUILD_GEN1INT_LSDALTON
  use gen1int_host
#endif
  private
!  private ::  test_first_geoderiv_overlap,&
!            & test_first_geoderiv_H1_mat,&
!            & test_first_geoderiv_refDmat,&
!            & test_first_geoderiv_Coulomb,&
!            & test_first_geoderiv_Exchange
  public  ::  dummy_subroutine_hessian_test,&
            & test_Hessian_contributions

CONTAINS

  SUBROUTINE dummy_subroutine_hessian_test()
  END SUBROUTINE dummy_subroutine_hessian_test

!> \brief Tests the different contributions used to build the Hessian tensor
!> \author \latexonly P. Merlot  \endlatexonly
!> \date 2012-09-13
!> \param F The Fock matrix
!> \param D The reference density matrix
!> \param ndmat The number of density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
!> \param iprint the printlevel, determining how much output should be generated
SUBROUTINE test_HESSIAN_contributions(F,D,Natoms,ndmat,setting,lupri,luerr)
  !
  IMPLICIT NONE
  Integer,INTENT(IN)            :: Natoms,lupri,luerr
  TYPE(LSSETTING),intent(INOUT) :: setting
  Type(matrix),INTENT(IN)       :: F,D
  Type(matrix),pointer :: Sa(:),Da(:),genSa(:) ! derivative along x,y and z for each atom
  Type(matrix),pointer :: Ka(:)
  Integer              :: i,nbast,ndmat
  Real(realk)          :: thresh1,thresh2,thresh3
  Real(realk)          :: ts,te
  Real(realk),pointer  :: reOrtho(:,:),exchangeGrad(:,:)
  Type(matrix),target  :: temp1,temp2
  type(matrixp)        :: DFDmat(ndmat)
  Type(matrixp)        :: Dmat(ndmat)
  Real(realk)          :: diffx,diffy,diffz
  !
  IF(ndmat.NE.1)  call lsquit('option not verified in test_Hessian_contributions()',-1)
  call lstimer('START ',ts,te,lupri)
  nbast = D%nrow
  thresh1 =1.0E-15 ! VerifyMatrices() threshold testing Sa, H1a
  thresh2=1.0E-15 ! Testing Da within reOrtho. grad. term
  thresh3=1.0E-10 ! Testing Ja,Ka

  ! Testing the geometric first derivative of the overlap matrix vs gen1int matrices
  call test_first_geoderiv_overlap(Natoms,nbast,thresh1,setting,lupri,luerr)

  ! Testing the geometric first derivative of the reference density matrix 
  ! using a gradient contribution
  ! \latexonly
  !   $Tr(\textbf{FD^a_0}) =Tr(\textbf{FD_0S^aD_0})=Tr(\textbf{D_0FD_0S^a})$ 
  ! \endlatexonly
  call test_first_geoderiv_refDmat(F,D,Natoms,ndmat,thresh2,setting,lupri,luerr)

  ! Testing the geometric first derivative of the one-electron Hamil. matrix
  ! vs gen1int matrices
  call test_first_geoderiv_H1_mat(Natoms,nbast,thresh1,setting,lupri,luerr)

  ! Testing the  geometric first derivative of the Coulomb term using a gradient contribution
  ! \latexonly $\frac{1}{2} Tr( \textbf{D J^a(D)} )$  \endlatexonly
  call test_first_geoderiv_Coulomb(D,Natoms,ndmat,thresh3,setting,lupri,luerr)

  ! Testing the  geometric first derivative of the Exchange term using a gradient contribution
  ! \latexonly $\frac{1}{2} Tr( \textbf{D K^a(D)} )$  \endlatexonly
  call test_first_geoderiv_Exchange(D,Natoms,ndmat,thresh3,setting,lupri,luerr)

  
  call lstimer('testHessian',ts,te,lupri)
END SUBROUTINE test_HESSIAN_contributions

!> \brief Test the first derivative of the overlap matrix vs gen1int routine
!> \author \latexonly P. Merlot  \endlatexonly
!> \date 2012-09-19
!> \param Natoms The nb. of atoms in the molecule
!> \param nbast The rank of the density matrix (nb. atomic orbitals)
!> \param Thresh The threshold used for comparing the resulting matrix elements
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE test_first_geoderiv_overlap(Natoms,nbast,thresh,setting,lupri,luerr)
  !
  IMPLICIT NONE
  Integer,INTENT(IN)            :: Natoms,lupri,luerr,nbast
  TYPE(LSSETTING),intent(INOUT) :: setting
  Real(realk),intent(IN)        :: THRESH
  !
  Type(matrix),pointer :: genSa(:) 
  Type(matrix),pointer :: Sa(:) ! derivative along x,y and z for each atom
  Integer              :: i
  Real(realk)          :: ts,te
  !
  call lstimer('START ',ts,te,lupri)
  !------------------------------------------------------
  ! Testing the geometric first derivative of the overlap
  ! matrix 
  call mem_alloc(Sa,3*Natoms)
  DO i=1,3*Natoms
     call mat_init(Sa(i),nbast,nbast)
  ENDDO
  call get_first_geoderiv_overlap(Sa,Natoms,setting,lupri,luerr)

#ifdef BUILD_GEN1INT_LSDALTON 
! PATRICK would you pleas ensure that your fucking code compiles!!!!

  ! test only possible if gen1int available
  call mem_alloc(genSa,3*nAtoms)
  DO i=1,3*Natoms
     call mat_init(genSa(i),nbast,nbast)
     call mat_zero(genSa(i))
  ENDDO
  call gen1int_host_get_first_geoderiv_overlap(setting,genSa,nAtoms,lupri)
  DO i=1,3*Natoms
     call VerifyMatrices(Sa(i),genSa(i),'genSa',THRESH,lupri)
     call mat_free(genSa(i))
  ENDDO
  call mem_dealloc(genSa)
  write (lupri,*) '   - geometric first derivative of the overlap matrix: OK'
  write (*,*)     '   - geometric first derivative of the overlap matrix: OK'
#endif
  DO i=1,3*Natoms
     call mat_free(Sa(i))
  ENDDO
  call mem_dealloc(Sa)
  call lstimer('test_Sa',ts,te,lupri)
END SUBROUTINE test_first_geoderiv_overlap



!> \brief Test the first derivative of the one-electron Hamiltonian
!>        vs gen1int routine
!> \author P. Merlot
!> \date 2013-07-10
!> \param Natoms The nb. of atoms in the molecule
!> \param nbast The rank of the density matrix (nb. atomic orbitals)
!> \param Thresh The threshold used for comparing the resulting matrix elements
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE test_first_geoderiv_H1_mat(Natoms,nbast,thresh,setting,lupri,luerr)
  !
  IMPLICIT NONE
  Integer,INTENT(IN)            :: Natoms,lupri,luerr,nbast
  TYPE(LSSETTING),intent(INOUT) :: setting
  Real(realk),intent(IN)        :: THRESH
  !
  Type(matrix),pointer :: genHa(:) 
  Type(matrix),pointer :: ha(:) ! derivative along x,y and z for each atom
  Integer              :: i
  Real(realk)          :: ts,te
  !
  call lstimer('START ',ts,te,lupri)
  !------------------------------------------------------
  ! Testing the geometric first derivative matrix
  call mem_alloc(ha,3*Natoms)
  DO i=1,3*Natoms
     call mat_init(ha(i),nbast,nbast)
  ENDDO
  call get_first_geoderiv_H1_mat(ha,Natoms,setting,lupri,luerr,0)
  

#ifdef BUILD_GEN1INT_LSDALTON 
  ! test only possible if gen1int available
  call mem_alloc(genHa,3*nAtoms)
  DO i=1,3*Natoms
     call mat_init(genHa(i),nbast,nbast)
     call mat_zero(genHa(i))
  ENDDO
  call gen1int_host_get_first_geoderiv_h1(setting,genHa,Natoms,lupri)
  DO i=1,3*Natoms
     call VerifyMatrices(Sa(i),genHa(i),'genHa',THRESH,lupri)
     call mat_free(genHa(i))
  ENDDO
  call mem_dealloc(genHa)
  write (lupri,*) '   - geometric first derivative of the one-electron Hamiltonian: OK'
  write (*,*)     '   - geometric first derivative of the one-electron Hamiltonian: OK'
#endif
  DO i=1,3*Natoms
     call mat_free(ha(i))
  ENDDO
  call mem_dealloc(ha)
  call lstimer('test_ha',ts,te,lupri)
END SUBROUTINE test_first_geoderiv_H1_mat


!> \brief Tests the geometric first derivative of the reference
!>        density matrix using the gradient reOrthonomalization term
!> \author \latexonly P. Merlot  \endlatexonly
!> \date 2012-09-19
!> \param F The Fock matrix
!> \param D The reference density matrix
!> \param Natoms The nb. of atoms in the molecule
!> \param ndmat The number of density matrices
!> \param Thresh The threshold used for comparing the resulting matrix elements
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE test_first_geoderiv_refDmat(F,D,Natoms,ndmat,thresh,setting,lupri,luerr)
  !
  IMPLICIT NONE
  Integer,INTENT(IN)            :: Natoms,lupri,luerr
  TYPE(LSSETTING),intent(INOUT) :: setting
  Type(matrix),INTENT(IN)       :: F,D
  Real(realk)                   :: THRESH
  Type(matrix),pointer :: Sa(:),Da(:) ! derivative along x,y and z for each atom
  Integer              :: i,nbast,ndmat
  Real(realk)          :: ts,te
  Real(realk),pointer  :: reOrtho(:,:),exchangeGrad(:,:)
  Type(matrix),target  :: temp1,temp2
  type(matrixp)        :: DFDmat(ndmat)
  Type(matrixp)        :: Dmat(ndmat)
  Real(realk)          :: diffx,diffy,diffz
  !
  call lstimer('START ',ts,te,lupri)
  nbast = D%nrow
  call mat_init(temp1,nbast,nbast)
  call mat_init(temp2,nbast,nbast)
  call mat_mul(D,F,'N','N',1E0_realk,0E0_realk,temp1)
  call mat_mul(temp1,D,'N','N',1E0_realk,0E0_realk,temp2)
  DFDmat(1)%p => temp2  ! DFD matrix calculated
  call mat_free(temp1)
  !
  call mem_alloc(reOrtho,3,Natoms)
  call II_get_reorthoNormalization(reortho,DFDmat,ndmat,setting,lupri,luerr) ! Tr(FDa)=Tr(FDSaD)=Tr(DFDSa)
  !
  call mem_alloc(Da,3*Natoms)
  DO i=1,3*Natoms
     call mat_init(Da(i),nbast,nbast)
  ENDDO
  call mem_alloc(Sa,3*Natoms)
  DO i=1,3*Natoms
     call mat_init(Sa(i),nbast,nbast)
  ENDDO
  call get_first_geoderiv_overlap(Sa,Natoms,setting,lupri,luerr)
  call get_first_geoderiv_refDmat(Da,D,Sa,Natoms,setting,lupri,luerr,0) ! Da=DSaD
  DO i=1,Natoms
     diffx = ABS(mat_trAB(F,Da(3*(i-1)+1))+reOrtho(1,i))
     diffy = ABS(mat_trAB(F,Da(3*(i-1)+2))+reOrtho(2,i))
     diffz = ABS(mat_trAB(F,Da(3*(i-1)+3))+reOrtho(3,i))
     IF(diffx.GT.thresh) THEN
        write(*,*) "   - wrong contribution to Tr(FDa)=Tr(DFDSa) DIFF(x",i,"): ",diffx
        call lsquit('reOrtho(x,i) grad term not well reproduced with Sa',lupri)
     ELSEIF(diffy.GT.thresh) THEN
        write(*,*) "   - wrong contribution to Tr(FDa)=Tr(DFDSa) DIFF(y",i,"): ",diffy
        call lsquit('reOrtho(y,i) grad term not well reproduced with Sa',lupri)
     ELSEIF(diffz.GT.thresh) THEN
        write(*,*) "   - wrong contribution to Tr(FDa)=Tr(DFDSa) DIFF(z",i,"): ",diffz
        call lsquit('reOrtho(z,i) grad term not well reproduced with Sa',lupri)
     ENDIF
  ENDDO
  write (lupri,*) '   - geometric first derivative of the ref. density matrix: OK'
  write (*,*)     '   - geometric first derivative of the ref. density matrix: OK'          
  call mem_dealloc(reOrtho)
  call mat_free(temp2)
  DO i=1,3*Natoms
     call mat_free(Sa(i))
     call mat_free(Da(i))
  ENDDO
  call mem_dealloc(Sa)
  call mem_dealloc(Da)
  call lstimer('test_Da',ts,te,lupri)
END SUBROUTINE test_first_geoderiv_refDmat

!> \brief Tests the geometric first derivative of the Coulomb term
!>        using a gradient contribution
!> \author \latexonly P. Merlot  \endlatexonly
!> \date 2013-07-10
!> \param D The reference density matrix
!> \param ndmat The number of density matrices
!> \param setting Integral evalualtion settings
!> \param Thresh The threshold used for comparing the resulting matrix elements
!> \param lupri Default print unit
!> \param luerr Unit for error printing
!> \param iprint the printlevel, determining how much output should be generated
SUBROUTINE test_first_geoderiv_Coulomb(D,Natoms,ndmat,thresh,setting,lupri,luerr)
  !
  IMPLICIT NONE
  Integer,INTENT(IN)             :: Natoms,lupri,luerr
  Real(realk),INTENT(IN)         :: thresh
  TYPE(LSSETTING),intent(INOUT)  :: setting
  Type(matrix),INTENT(IN),target :: D
  Type(matrix),pointer :: Ja(:)  ! derivative along x,y and z for each atom
  Integer              :: i,nbast,ndmat
  Real(realk)          :: ts,te
  Real(realk),pointer  :: coulombGrad(:,:)
  Type(matrixp)        :: Dmat(ndmat)
  Real(realk)          :: diffx,diffy,diffz
  !
  IF(ndmat.NE.1)  call lsquit('option not verified in test_first_geoderivCoulomb()',-1)
  call lstimer('START ',ts,te,lupri)
  nbast = D%nrow

  ! Testing the geometric first derivative of the Coulomb 
  ! term using a gradient contribution
  ! \latexonly $\frac{1}{2} Tr( \textbf{D J^a(D)} )$  \endlatexonly
  Dmat(1)%p => D
  call mem_alloc(coulombGrad,3,Natoms)
  coulombGrad = 0E0_realk
  call  II_get_J_gradient(coulombGrad,Dmat,Dmat,ndmat,ndmat,setting,lupri,luerr)
  !
  call mem_alloc(Ja,3*Natoms)
  DO i=1,3*Natoms
     call mat_init(Ja(i),nbast,nbast)
  ENDDO
  call get_first_geoderiv_Coulomb_mat(Ja,D,Natoms,setting,lupri,luerr,4)
DO i=1,Natoms
      
     diffx = ABS(-0.5E0_realk*mat_trAB(D,Ja(3*(i-1)+1))+coulombGrad(1,i))
     diffy = ABS(-0.5E0_realk*mat_trAB(D,Ja(3*(i-1)+2))+coulombGrad(2,i))
     diffz = ABS(-0.5E0_realk*mat_trAB(D,Ja(3*(i-1)+3))+coulombGrad(3,i))
     IF(diffx.GT.thresh) THEN
        write(*,*) "   - wrong contribution to Tr(DJa) DIFF(x",i,")= ",diffx
        call lsquit('coulombGrad(x,i) grad term not well reproduced with Ja',lupri)
     ELSEIF(diffy.GT.thresh) THEN
        write(*,*) "   - wrong contribution to Tr(DJa) DIFF(y",i,"): ",diffy
        call lsquit('coulombGrad(y,i) grad term not well reproduced with Ja',lupri)
     ELSEIF(diffz.GT.thresh) THEN
        write(*,*) "   - wrong contribution to Tr(DJa) DIFF(z",i,"): ",diffz
        call lsquit('coulombGrad(z,i) grad term not well reproduced with Ja',lupri)
     ENDIF
  ENDDO
  write (lupri,*) '   - geometric first derivative of the Coulomb term: OK'
  write (*,*)     '   - geometric first derivative of the Coulomb term: OK'

  ! TESTS FINISHED: FREE MEMORY    
  DO i=1,3*Natoms
     call mat_free(Ja(i))
  ENDDO
  call mem_dealloc(Ja)
  call mem_dealloc(coulombGrad)
  call lstimer('test_Ja',ts,te,lupri)
END SUBROUTINE test_first_geoderiv_Coulomb

!> \brief Tests the geometric first derivative of the Exchage term
!>        using a gradient contribution
!> \author \latexonly P. Merlot  \endlatexonly
!> \date 2012-09-19
!> \param D The reference density matrix
!> \param ndmat The number of density matrices
!> \param setting Integral evalualtion settings
!> \param Thresh The threshold used for comparing the resulting matrix elements
!> \param lupri Default print unit
!> \param luerr Unit for error printing
!> \param iprint the printlevel, determining how much output should be generated
SUBROUTINE test_first_geoderiv_Exchange(D,Natoms,ndmat,thresh,setting,lupri,luerr)
  !
  IMPLICIT NONE
  Integer,INTENT(IN)             :: Natoms,lupri,luerr
  Real(realk),INTENT(IN)         :: thresh
  TYPE(LSSETTING),intent(INOUT)  :: setting
  Type(matrix),INTENT(IN),target :: D
  Type(matrix),pointer :: Ka(:)  ! derivative along x,y and z for each atom
  Integer              :: i,nbast,ndmat
  Real(realk)          :: ts,te
  Real(realk),pointer  :: exchangeGrad(:,:)
  Type(matrixp)        :: Dmat(ndmat)
  Real(realk)          :: diffx,diffy,diffz
  !
  IF(ndmat.NE.1)  call lsquit('option not verified in test_first_geoderivExchange()',-1)
  call lstimer('START ',ts,te,lupri)
  nbast = D%nrow

  ! Testing the  geometric first derivative of the Exchange 
  ! term using a gradient contribution
  ! \latexonly $\frac{1}{2} Tr( \textbf{D K^a(D)} )$  \endlatexonly
  Dmat(1)%p => D
  call mem_alloc(exchangeGrad,3,Natoms)
  exchangeGrad = 0E0_realk
  call  II_get_K_gradient(exchangeGrad,Dmat,Dmat,ndmat,ndmat,setting,lupri,luerr)
  !
  call mem_alloc(Ka,3*Natoms)
  DO i=1,3*Natoms
     call mat_init(Ka(i),nbast,nbast)
  ENDDO
  call get_first_geoderiv_exchange_mat(Ka,D,Natoms,setting,lupri,luerr,4)
  DO i=1,Natoms
     diffx = ABS(0.5E0_realk*mat_trAB(D,Ka(3*(i-1)+1))+exchangeGrad(1,i))
     diffy = ABS(0.5E0_realk*mat_trAB(D,Ka(3*(i-1)+2))+exchangeGrad(2,i))
     diffz = ABS(0.5E0_realk*mat_trAB(D,Ka(3*(i-1)+3))+exchangeGrad(3,i))
     IF(diffx.GT.thresh) THEN
        write(*,*) "   - wrong contribution to Tr(DKa) DIFF(x",i,"): ",diffx
        call lsquit('exchangeGrad(x,i) grad term not well reproduced with Ka',lupri)
     ELSEIF(diffy.GT.thresh) THEN
        write(*,*) "   - wrong contribution to Tr(DKa) DIFF(y",i,"): ",diffy
        call lsquit('exchangeGrad(y,i) grad term not well reproduced with Ka',lupri)
     ELSEIF(diffz.GT.thresh) THEN
        write(*,*) "   - wrong contribution to Tr(DKa) DIFF(z",i,"): ",diffz
        call lsquit('exchangeGrad(z,i) grad term not well reproduced with Ka',lupri)
     ENDIF
  ENDDO
  write (lupri,*) '   - geometric first derivative of the exchange term: OK'
  write (*,*)     '   - geometric first derivative of the exchange term: OK'

  ! TESTS FINISHED: FREE MEMORY    
  DO i=1,3*Natoms
     call mat_free(Ka(i))
  ENDDO
  call mem_dealloc(Ka)
  call mem_dealloc(exchangeGrad)
  call lstimer('test_Ka',ts,te,lupri)
END SUBROUTINE test_first_geoderiv_Exchange

!------------ END MODULE ------------
END MODULE test_molecular_hessian_mod
