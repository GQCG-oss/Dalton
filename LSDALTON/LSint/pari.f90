
!> @file 
!> Contains pari-atomic resolution-of-the-identity (PARI) specific routines

!> \brief Contains PARI routines
MODULE pari_mod
use precision
use Typedef  
use ls_Integral_Interface
use AtomSparse
use mat3d_mod
use Integralparameters
use molecule_module
use memory_handling
use lstiming
CONTAINS

!> \brief Calculates the PARI (pair-atomic resolution-of-the-identity) fitting coefficients 
!> \latexonly
!>        $c_\alpha^{ab} = (\alpha|\beta)^{-1} (\beta|ab), \quad \alpha\in{A\cup B},
!>                         \quad a\in A,\quad b\in B$
!> \endlatexonly
!> \author S. Reine and P. Merlot
!> \date 2010-07-22
SUBROUTINE getPariCoefficients(LUPRI,LUERR,SETTING,calpha_ab,orbitalInfo,regCSfull,auxCSfull) 
    implicit none
    TYPE(LSSETTING),intent(inout)         :: SETTING
    INTEGER,intent(IN)                    :: LUPRI,LUERR
    TYPE(MOLECULARORBITALINFO),INTENT(IN) :: orbitalInfo
    TYPE(MAT3D),intent(inout)             :: calpha_ab(orbitalInfo%nAtoms)
    TYPE(LSTENSOR),pointer                :: regCSfull,auxCSfull
        !
    Integer               :: iAtomA,iAtomB,nRegA,nRegB,nAuxA,nAuxB
    Integer               :: iAlpha,iGamma,iRegA,iRegB
    Integer               :: startRegA,startRegB,startAuxA,startAuxB
    Integer               :: endRegA,endRegB,endAuxA,endAuxB
    Integer               :: nAux
    Integer               :: info
    Real(realk),pointer   :: alpha_ab(:,:,:)
    Real(realk),pointer   :: alphaBetaFull(:,:)
    Real(realk),pointer   :: alpha_ab_red(:,:)
    Real(realk),pointer   :: alphaBeta(:,:),copy_alpBeta(:,:)
    Real(realk),pointer   :: eigValphaBeta(:)
!   TYPE(MOLECULEINFO)    :: atoms(orbitalInfo%nAtoms)
    !
    TYPE(MOLECULEINFO),pointer :: molecule
    TYPE(MOLECULEINFO),pointer :: ATOMS_A(:)
    Integer           :: nAtoms,iRegAfull,iRegBfull,iBeta
    Integer,pointer   :: indexAB(:,:)
    Integer           :: nRegAB,iRegAB
    Real(realk)       :: maxAbs
    Real(realk)       :: threshold = 1E-20_realk ! 1E-10 NOK while {1E-20, 1E-30} OK
    Character(80)     :: C_filename
    Real(realk),pointer :: alpha_charge(:)
    Real(realk),pointer :: K_red(:)
    Real(realk),pointer :: lambda_ab(:,:)        ! multipole constraint
    ! Real(realk),pointer :: lambda_ab(:)        ! charge constraint
    Real(realk),pointer :: charge_ab(:,:)
    Real(realk),pointer :: dipole_ab(:,:,:)
    Real(realk),pointer :: alpha_dipole(:,:)
    Real(realk),pointer :: Amat(:,:)
    Real(realk),pointer :: multipoleDiff(:,:)
    Real(realk)         :: charge_diff_max,xdipole_diff_max,ydipole_diff_max,zdipole_diff_max
    Integer             :: iauxa,ireg,iaux,iauxb,ksi1,ksi2
    CHARACTER(LEN=3)    :: nline
    INTEGER             :: nb_eigVal
    Real(realk)         :: minEigv,maxEigv,minTemp,maxTemp,conditionNum,tempCondNum
    Real(realk),pointer :: eigval(:)

    minEigV = 9999.90E0_realk
    maxEigV = -9999.90E0_realk
    conditionNum = abs(maxEigV)/abs(minEigV)

    threshold = SETTING%SCHEME%PARI_THRESHOLD*SETTING%SCHEME%THRESHOLD
    
    !write(*,*) 'used PARI_THRES',threshold

    molecule => SETTING%MOLECULE(1)%p      ! --- prepare for re-initialization of SETTING
    nAtoms = molecule%nAtoms
    DO iAtomA=1,orbitalInfo%nAtoms
       call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)
       call init_MAT3D(calpha_ab(iAtomA),nAuxA,nRegA,orbitalInfo%nBastReg)
    ENDDO

    C_filename = 'CALPHA_AB'
    IF (io_file_exist(C_filename,setting%IO)) THEN
       call io_read_mat3d(calpha_ab,nAtoms,C_filename,setting%io,lupri,luerr)
    ELSE
       IF (setting%scheme%pari_dipole) THEN 
          write(lupri,'(1X,A)') 'Charge- and dipole-constrained PARI coefficients'
       ELSEIF (setting%scheme%pari_charge) THEN
          write(lupri,'(1X,A)') 'Charge-constrained PARI coefficients'
       ENDIF

       allocate(ATOMS_A(nAtoms))  ! --- each molecule made of only one atom

       call pari_set_atomic_fragments(molecule,ATOMS_A,nAtoms,lupri)
       
       ! --- For each "Pair of atoms A and B",
       DO iAtomA=1,orbitalInfo%nAtoms
          call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)
          DO iAtomB=iAtomA,orbitalInfo%nAtoms
             call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
             nAux      = nAuxA + nAuxB
             IF (iAtomA.EQ.iAtomB) nAux = nAuxA
             call mem_alloc(alpha_ab,nAux,nRegA,nRegB)
             call mem_alloc(alphaBeta,nAux,nAux)
             call ls_dzero(alpha_ab,nAux*nRegA*nRegB)
             call ls_dzero(alphaBeta,nAux*nAux)

             ! --- extract (alpha | beta) and (alpha | ab) with 
             ! --- alpha in {A U B}, a in {A}, b in {B} from the full matrices
             alpha_ab = 0E0_realk
             call pari_alphaab(  alpha_ab,setting,molecule,atoms_a,iAtomA,iAtomB,&
                  &                    nAuxA,nAuxB,nRegA,nRegB,regCSfull,auxCSfull,lupri,luerr)

             call pari_alphaBeta(alphaBeta,setting,molecule,atoms_a,iAtomA,iAtomB,&
                  &                    nAuxA,nAuxB,regCSfull,auxCSfull,lupri,luerr)
             
!              !write(*,*) '(alpha|beta) on AB'
!              iaux = 0
!              DO iauxa=1,nAux
!                 DO iauxb=1,nAux
!                    iaux=iaux+1
!                    nline='no'
!                    if(iaux .eq. nAux) THEN
!                       iaux=0
!                       nline='yes'
!                    ENDIF
!                    ! write(*,'(E12.4)',advance=nline)  alphaBeta(iauxa,iauxb)
!                 ENDDO
!              ENDDO

#if 0
             ! Checking eigenvalues of (alpha|beta) matrices
             call mem_alloc(copy_alpBeta,nAux,nAux)
             call ls_dzero(copy_alpBeta,nAux*nAux)
             call mem_alloc(eigValphaBeta,nAux)
             call ls_dzero(eigValphaBeta,nAux)
             copy_alpBeta = alphaBeta
             call II_get_eigv_square_mat(lupri,luerr,copy_alpBeta,eigValphaBeta,nAux) 
             !
             call check_min_max_Array_elem(minTemp,maxTemp,tempCondNum,eigValphaBeta,nAux,lupri,luerr)
             minEigV=min(minEigV,minTemp)
             maxEigV=max(maxEigV,maxTemp)
             conditionNum = max(conditionNum,tempCondNum)
             call mem_dealloc(copy_alpBeta)
             call mem_dealloc(eigValphaBeta)
#endif

             ! --- Solve only for (alpha|ab) > threshold, make first reduction of integrals according to
             ! --- size of |(alpha|ab)|
             nRegAB = 0
             DO iRegB=1,nRegB
                DO iRegA=1,nRegA
                   maxAbs = abs(alpha_ab(1,iRegA,iRegB))
                   DO iAlpha=2,nAux
                      maxAbs = max(maxAbs,abs(alpha_ab(iAlpha,iRegA,iRegB)))
                   ENDDO
                   IF (maxAbs.GT.threshold) THEN
                      nRegAB = nRegAB + 1
                   ENDIF
                ENDDO
             ENDDO
             ! write(*,'(A,I2,A,I2,A,I2,A,I2,A,I2)') 'iA=',iAtomA,', iB=',iAtomB,', nRegAB=',nregAB,', nRegA/nRegB=',nRegA,'/',nRegB
             IF (nRegAB.GT. 0) THEN
                call mem_alloc(indexAB,2,nRegAB)

                IF (setting%scheme%pari_dipole) THEN   
                   ! We need 4 additional dimensions when we constrain the dipole (charge + dipole-{x,y,z})
                   call mem_alloc(lambda_ab,4,nRegAB)
                   call mem_alloc(dipole_ab,nRegA,nRegB,4)
                   !  Calculate regular charges and dipoles (dipole_ab = {<ab>,<xab>,<yab>,<zab>})
                   call pari_dipole_ab(dipole_ab,setting,atoms_a(iAtomA),atoms_a(iAtomB),&
                        &                          nRegA,nRegB,lupri,luerr)

!                    write(*,*) 'dipole_ab(:,1)'
!                    iaux = 0
!                    DO iregA=1,nregA
!                       DO iregB=1,nregB
!                          DO iauxb=1,4
!                             iaux=iaux+1
!                             nline='no'
!                             if(iaux .eq. 4) THEN
!                                iaux=0
!                                nline='yes'
!                             ENDIF
!                             write(*,'(E12.4)',advance=nline)  dipole_ab(irega,iregab,iaux)
!                          ENDDO
!                       ENDDO
!                    ENDDO



                   nRegAB = nRegAB + 4
                   ! We need one additional dimension when we constrain the charge
                ELSEIF (setting%scheme%pari_charge) THEN
                   call mem_alloc(lambda_ab,1,nRegAB)
                   call mem_alloc(charge_ab,nRegA,nRegB)
                   !           Calculate regular charges (charge_ab = <ab>)
                   call pari_charge_ab(charge_ab,setting,atoms_a(iAtomA),atoms_a(iAtomB),&
                        &                          nRegA,nRegB,lupri,luerr)
                   nRegAB = nRegAB + 1
                ENDIF ! if pari_charge or pari_dipole

                call mem_alloc(alpha_ab_red,nAux,nRegAB)
                alpha_ab_red = 0E0_realk
                iRegAB = 0
                DO iRegB=1,nRegB
                   DO iRegA=1,nRegA
                      maxAbs = abs(alpha_ab(1,iRegA,iRegB))
                      DO iAlpha=2,nAux
                         maxAbs = max(maxAbs,abs(alpha_ab(iAlpha,iRegA,iRegB)))
                      ENDDO
                      IF (maxAbs.GT.threshold) THEN
                         iRegAB = iRegAB + 1
                         indexAB(1,iRegAB) = iRegA
                         indexAB(2,iRegAB) = iRegB
                         alpha_ab_red(:,iRegAB) = alpha_ab(:,iRegA,iRegB)
                         !               Put charge_ab, dipole_ab, ... into lambda_ab
                         IF (setting%scheme%pari_dipole) THEN 
                            lambda_ab(1,iRegAB) = dipole_ab(iRegA,iRegB,1) ! charge
                            lambda_ab(2,iRegAB) = dipole_ab(iRegA,iRegB,2) ! x-dipole
                            lambda_ab(3,iRegAB) = dipole_ab(iRegA,iRegB,3) ! y-dipole
                            lambda_ab(4,iRegAB) = dipole_ab(iRegA,iRegB,4) ! z-dipole
                         ELSEIF (setting%scheme%pari_charge) THEN
                            lambda_ab(1,iRegAB) = charge_ab(iRegA,iRegB)
                         ENDIF ! if pari_charge or pari_dipole
                      ENDIF
                   ENDDO
                ENDDO

                ! --- Calculate the auxiliary charges and dipoles (alpha_dipole={ <alpha>, <xalpha>, <yalpha>, <zalpha>})          
                ! --- and store them in additional columns of alpha_ab_red(:,:)
                IF (setting%scheme%pari_dipole) THEN
                   call mem_alloc(alpha_dipole,nAux,4)
                   call pari_alpha_dipole(alpha_dipole,setting,molecule,atoms_a,iAtomA,iAtomB,&
                        &                             nAuxA,nAuxB,lupri,luerr)
                   alpha_ab_red(:,nRegAB)   = alpha_dipole(:,4) ! z-dipole fit
                   alpha_ab_red(:,nRegAB-1) = alpha_dipole(:,3) ! y-dipole fit 
                   alpha_ab_red(:,nRegAB-2) = alpha_dipole(:,2) ! x-dipole fit
                   alpha_ab_red(:,nRegAB-3) = alpha_dipole(:,1) ! charge fit
                ELSE IF (setting%scheme%pari_charge) THEN
                   ! Calculate the auxiliary charges (alpha_charge=<alpha>)
                   call mem_alloc(alpha_charge,nAux)
                   call mem_alloc(K_red,nRegAB)
                   call pari_alpha_charge(alpha_charge,setting,molecule,atoms_a,iAtomA,iAtomB,&
                        &                             nAuxA,nAuxB,lupri,luerr)
                   alpha_ab_red(:,nRegAB) = alpha_charge
                ENDIF ! if pari_charge or pari_dipole

                ! --- Solve the linear system to find the unconstrained density-fitting coeff.
                ! --- c_alpha^ab = (alpha|beta)^-1 (beta|ab) with alpha in AUB
                ! --- (here the alpha_ab is changed to be c_alpha^ab)
                call DPOSV('U',nAux,nRegAB,alphaBeta,nAux,alpha_ab_red,nAux,info)
                If (info.ne. 0) THEN
                   WRITE(LUPRI,'(1X,A,I5)') 'DPOSV error in getPariCoefficients 1/2. Info =',info
                   call LSQUIT('DPOSV error in getPariCoefficients 1/2',lupri)
                ENDIF

                ! --- find the Lagrange multiplier(s)
                IF (setting%scheme%pari_dipole) THEN
                   ! --- Solve the linear system to find the Langrangian multipliers
                   ! --- solve Amat*lambda_ab = B_ab (B_ab stored in lambda_ab)
                   call mem_alloc(Amat,4,4)
                   Amat = 0E0_realk
                   DO ksi1=1,4
                      DO ksi2=1,4
                         DO iAlpha=1,nAux
                            Amat(ksi1,ksi2) = Amat(ksi1,ksi2) + alpha_dipole(iAlpha,ksi1)*alpha_ab_red(iAlpha,nRegAB-4+ksi2)
                         ENDDO
                      ENDDO
                   ENDDO
                   DO iregAB=1,nregAB-4
                     DO iAux=1,nAux
                        lambda_ab(1,iRegAB) = lambda_ab(1,iRegAB) - alpha_dipole(iAux,1)*alpha_ab_red(iAux,iRegAB)
                        lambda_ab(2,iRegAB) = lambda_ab(2,iRegAB) - alpha_dipole(iAux,2)*alpha_ab_red(iAux,iRegAB)
                        lambda_ab(3,iRegAB) = lambda_ab(3,iRegAB) - alpha_dipole(iAux,3)*alpha_ab_red(iAux,iRegAB)
                        lambda_ab(4,iRegAB) = lambda_ab(4,iRegAB) - alpha_dipole(iAux,4)*alpha_ab_red(iAux,iRegAB)
                     ENDDO
                   ENDDO

                   call DPOSV('U',4,nRegAB-4,Amat,4,lambda_ab,4,info)
                   If (info.ne. 0) THEN
                      WRITE(LUPRI,'(1X,A,I5)') 'DPOSV error in getPariCoefficients 2/2. Info =',info
                      call LSQUIT('DPOSV error in getPariCoefficients 2/2',lupri)
                   ENDIF

                   ! --- rewrite the charge-/dipole-constrained fitting coefficients in 
                   ! --- terms of the unconstrained fitting coeff. and the corrections terms
                   DO iRegAB=1,nRegAB-4
                      DO iAux=1,nAux
                         alpha_ab_red(iAux,iRegAB) = alpha_ab_red(iAux,iRegAB) &
                              &                   + lambda_ab(1,iRegAB)*alpha_ab_red(iAux,nRegAB-3) &
                              &                   + lambda_ab(2,iRegAB)*alpha_ab_red(iAux,nRegAB-2) &
                              &                   + lambda_ab(3,iRegAB)*alpha_ab_red(iAux,nRegAB-1) &
                              &                   + lambda_ab(4,iRegAB)*alpha_ab_red(iAux,nRegAB)
                     ENDDO
                   ENDDO
                   charge_diff_max = 0E0_realk
                   xdipole_diff_max = 0E0_realk
                   ydipole_diff_max = 0E0_realk
                   zdipole_diff_max = 0E0_realk
                   DO iRegAB = 1,nRegAB-4
                      iRegA = indexAB(1,iRegAB)
                      iRegB = indexAB(2,iRegAB)
                      lambda_ab(1,iRegAB) = dipole_ab(iRegA,iRegB,1)
                      lambda_ab(2,iRegAB) = dipole_ab(iRegA,iRegB,2)
                      lambda_ab(3,iRegAB) = dipole_ab(iRegA,iRegB,3)
                      lambda_ab(4,iRegAB) = dipole_ab(iRegA,iRegB,4)
                      DO iAlpha=1,nAux
                         lambda_ab(1,iRegAB) = lambda_ab(1,iRegAB) - alpha_dipole(iAlpha,1)*alpha_ab_red(iAlpha,iRegAB)
                         lambda_ab(2,iRegAB) = lambda_ab(2,iRegAB) - alpha_dipole(iAlpha,2)*alpha_ab_red(iAlpha,iRegAB)
                         lambda_ab(3,iRegAB) = lambda_ab(3,iRegAB) - alpha_dipole(iAlpha,3)*alpha_ab_red(iAlpha,iRegAB)
                         lambda_ab(4,iRegAB) = lambda_ab(4,iRegAB) - alpha_dipole(iAlpha,4)*alpha_ab_red(iAlpha,iRegAB)
                      ENDDO
                      charge_diff_max  = max(charge_diff_max ,abs(lambda_ab(1,iRegAB)))
                      xdipole_diff_max = max(xdipole_diff_max,abs(lambda_ab(2,iRegAB)))
                      ydipole_diff_max = max(ydipole_diff_max,abs(lambda_ab(3,iRegAB)))
                      zdipole_diff_max = max(zdipole_diff_max,abs(lambda_ab(4,iRegAB)))
                   ENDDO
                   IF (charge_diff_max.GE. 1E-13_realk)  &
                        & write(LUPRI,*) 'Maximum error in PARI charge-constraint: ',charge_diff_max
                   IF (xdipole_diff_max.GE. 1E-13_realk) & 
                        & write(LUPRI,*) 'Maximum error in PARI x-dipole-constraint: ',xdipole_diff_max
                   IF (ydipole_diff_max.GE. 1E-13_realk) & 
                        & write(LUPRI,*) 'Maximum error in PARI y-dipole-constraint: ',ydipole_diff_max
                   IF (zdipole_diff_max.GE. 1E-13_realk) & 
                        & write(LUPRI,*) 'Maximum error in PARI z-dipole-constraint: ',zdipole_diff_max
                   !
                   call mem_dealloc(dipole_ab)
                   call mem_dealloc(lambda_ab)
                   call mem_dealloc(Amat)
                   call mem_dealloc(alpha_dipole)
                   nRegAB = nRegAB-4

                ELSEIF (setting%scheme%pari_charge) THEN
                   !           -- K_red contrains sum_I <I>c_I^{ab} and sum_I <I> k_I (needed to obtain lambda_ab
                   call dgemm('N','N',1,nRegAB,nAux,1E0_realk,alpha_charge,1,alpha_ab_red,nAux,&
                        &                 0E0_realk,K_red,1)
                   lambda_ab(1,:) = lambda_ab(1,:) - K_red(1:nRegAB-1)
                   lambda_ab(1,:) = lambda_ab(1,:) / K_red(nRegAB)

                   ! --- rewrite the charge-/dipole-constrained fitting coefficients in 
                   ! --- terms of the unconstrained fitting coeff. and the corrections terms
                   DO iRegAB=1,nRegAB-1
                      alpha_ab_red(:,iRegAB) = alpha_ab_red(:,iRegAB) + lambda_ab(1,iRegAB)*alpha_ab_red(:,nRegAB)
                   ENDDO
                   charge_diff_max = 0E0_realk
                   DO iRegAB = 1,nRegAB-1
                      iRegA = indexAB(1,iRegAB)
                      iRegB = indexAB(2,iRegAB)
                      lambda_ab(1,iRegAB) = charge_ab(iRegA,iRegB)
                      DO iAlpha=1,nAux
                         lambda_ab(1,iRegAB) = lambda_ab(1,iRegAB) - alpha_charge(iAlpha)*alpha_ab_red(iAlpha,iRegAB)
                      ENDDO
                      charge_diff_max = max(charge_diff_max,abs(lambda_ab(1,iRegAB)))
                   ENDDO
                   IF (charge_diff_max.GE. 1E-13_realk) THEN
                      !              write(*,*) 'error in charge-constraint:charge-diff',lambda_ab
                   ENDIF
                   call mem_dealloc(charge_ab)
                   call mem_dealloc(lambda_ab)
                   call mem_dealloc(K_red)
                   call mem_dealloc(alpha_charge)
                   nRegAB = nRegAB-1
                ENDIF ! if pari_charge or pari_dipole

                DO iRegAB=1,nRegAB
                   iRegA=indexAB(1,iRegAB)
                   iRegB=indexAB(2,iRegAB)
                   iRegAfull = startRegA+iRegA-1
                   iRegBfull = startRegB+iRegB-1
                   calpha_ab(iAtomA)%elements(:,iRegA,iRegBfull) = alpha_ab_red(1:nAuxA,iRegAB)
                   IF (iAtomA.NE.iAtomB) calpha_ab(iAtomB)%elements(:,iRegB,iRegAfull) = alpha_ab_red(nAuxA+1:nAuxA+nAuxB,iRegAB)
                ENDDO

                !          ! debug
                !           ireg = 0
                !           DO iauxa=1,naux
                ! !             write(*,'(A20,I3)') 'Calpha_ab_AB iaux=',iauxa
                !              DO irega=1,nregab
                !                 ireg=ireg+1
                !                 nline='no'
                !                 if(ireg .eq. nregb) THEN
                !                    ireg=0
                !                    nline='yes'
                !                 ENDIF
                ! !                write(*,'(E12.4)',advance=nline)  alpha_ab_red(iauxa,irega)
                !              ENDDO
                !           ENDDO
                call mem_dealloc(indexAB)
                call mem_dealloc(alpha_ab_red)
             ENDIF ! IF nRegAB >0


             call mem_dealloc(alpha_ab)
             call mem_dealloc(alphaBeta)
          ENDDO !iAtomB
       ENDDO !iAtomA
       ! Reset settings to default

#if 0
       write(lupri,*) "(alpha|beta): minEigV of all (alpha|beta) local matrices: ",minEigV
       write(*,*)     "(alpha|beta): minEigV of all (alpha|beta) local matrices: ",minEigV
       write(lupri,*) "(alpha|beta): maxEigV of all (alpha|beta) local matrices: ",maxEigV
       write(*,*)     "(alpha|beta): maxEigV of all (alpha|beta) local matrices: ",maxEigV
       write(lupri,*) "(alpha|beta) full: Condition Number (abs(max)/abs(min): ",conditionNum
       write(*,*)     "(alpha|beta) full: Condition Number (abs(max)/abs(min): ",conditionNum
#endif

       call typedef_setMolecules(setting,molecule,1,2,3,4)

       call pari_free_atomic_fragments(ATOMS_A,nAtoms)
       deallocate(ATOMS_A)  
       CALL io_add_filename(setting%io,C_filename,lupri)
       CALL io_write_mat3d(calpha_ab,nAtoms,C_filename,setting%io,lupri,luerr)
    ENDIF ! if file CALPHA exists
    !call LSQUIT('Testing eigenvalues of the (alpha|beta) matrices - quitting getPariCoefficients()',-1)
  END SUBROUTINE getPariCoefficients

SUBROUTINE freePariCoefficients(calpha_ab,nAtoms)
implicit none
Integer,intent(in)        :: nAtoms
type(MAT3D),intent(inout) :: calpha_ab(nAtoms)
!
Integer iAtom
!
DO iAtom=1,nAtoms
   call free_MAT3D(calpha_ab(iAtom))
ENDDO
END SUBROUTINE freePariCoefficients


!> \brief This subroutine constructs the density coefficients centered on one atom A for PARI-K
!> \latexonly
!>   $d_\alpha^{ac} = \sum_b c_\alpha^{ab} D_{bc}, \quad a,\alpha\in A, \quad b,c\in{\cal M}$
!> \endlatexonly
!> \author P. Merlot, S. Reine
!> \date 2010-05-04
!> \param orbitalInfo Information about the molecular orbitals (dimensions)
!> \param calpha_ab Array of extracted PARI coefficients centered on each atom
!> \param Dfull the full density matrix
!> \param dalpha_ab Array of PARI density coefficients (alpha,a in A, and b in mol.)
SUBROUTINE getDensityCoefficients(iAtomA,orbitalInfo,calpha_ab,Dfull,dalpha_ab)
IMPLICIT NONE
Integer,INTENT(IN)                     :: iAtomA
TYPE(MOLECULARORBITALINFO),INTENT(IN)  :: orbitalInfo
TYPE(MAT3D),intent(IN)                 :: calpha_ab(orbitalInfo%nAtoms)
Real(realk),intent(IN)                 :: Dfull(orbitalInfo%nBastReg,orbitalInfo%nBastReg,1)
Real(realk),intent(INOUT)                :: dalpha_ab(orbitalInfo%numAtomicOrbitalsAux(iAtomA),&
     &                                              orbitalInfo%numAtomicOrbitalsReg(iAtomA),&
     &                                              orbitalInfo%nBastReg)
!
Integer               :: nRegA,startRegA,endRegA
Integer               :: nAuxA,startAuxA,endAuxA
Integer               :: nBasis
Integer               :: iAlpha

! --- for now the number of regular basis functions for B is set to its maximum
nBasis     = orbitalInfo%nBastReg
call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)

call DGEMM('N','N',nAuxA*nRegA,nBasis,nBasis,1E0_realk,calpha_ab(iAtomA)%elements,&
     &     nAuxA*nRegA,Dfull,nBasis,0E0_realk,dalpha_ab,nAuxA*nRegA)

END SUBROUTINE getDensityCoefficients

!> \brief The three-center PARI contributions to the exchange matrix
!> \latexonly
!>            $K_{ab} = \sum_{\alpha,cd} c_\alpha^{ac}(\alpha|bd) D_{cd} + (ac|\alpha)c_\alpha^{bd} D_{cd}$
!> \endlatexonly
!> \author S. Reine
!> \date   2010-07-22
!> \param Kfull      
!> \param Dfull      
!> \param dalpha     
!> \param alpha_beta 
!> \param alpha_ab   
!> \param iAtomA     
!> \param orbitalInfo
!> \param nBastReg   
!> \param nAtoms     
!> \param lupri      
SUBROUTINE pariK_ThreeCenter(Kfull,Dfull,setting,calpha_ab,dalpha_ad,&
     &                       iAtomA,orbitalInfo,nBastReg,nAtoms,regCSfull,auxCSfull,lupri,luerr)
implicit none
Real(realk),pointer                   :: Dfull(:,:,:)
Real(realk),pointer                   :: Kfull(:,:,:,:,:)
type(lssetting),intent(inout)         :: setting
Integer,intent(in)                    :: lupri,luerr
Integer,intent(in)                    :: iAtomA, nAtoms,nBastReg
TYPE(MOLECULARORBITALINFO),intent(in) :: orbitalInfo
Real(realk),pointer                   :: dalpha_ad(:,:,:)
TYPE(MAT3D),pointer                   :: calpha_ab(:)
TYPE(LSTENSOR),pointer                :: regCSfull,auxCSfull
!
Integer                    :: iAtomB,iAtomC,iAtomD
Real(realk),pointer        :: dalpha_ca(:,:,:)
Real(realk),pointer        :: alpha_cd(:,:,:)
Real(realk),pointer        :: alpha_cd_full(:,:,:)
Real(realk),pointer        :: alpha_dc(:,:,:)
Real(realk),pointer        :: extracted_alpha_gamma(:,:)
Real(realk),pointer        :: Xalpha_ac(:,:,:)
Real(realk),pointer        :: Xalpha_ad(:,:,:)
Real(realk),pointer        :: Xalpha_ca(:,:,:)
Real(realk),pointer        :: Xalpha_da(:,:,:)

Integer                    :: nRegA,nRegB,nRegC,nRegD,nAuxA,nAuxB,nAuxC,nAuxD
Integer                    :: startRegA,startRegB,startRegC,startRegD,startAuxA,startAuxB,startAuxC,startAuxD
Integer                    :: endRegA,endRegB,endRegC,endRegD,endAuxA,endAuxB,endAuxC,endAuxD
Integer                    :: iAlpha,iRegA,iRegB,iRegC,iRegD
!
TYPE(AtomSparseMat)        :: alphaBeta
TYPE(MoleculeInfo),pointer :: molecule
Real(realk)                :: ts,te
TYPE(MOLECULEINFO),pointer :: ATOMS(:)
Integer                    :: iAtom
Character(len=22)          :: FRAGMENTNAME
TYPE(LSTENSOR),pointer     :: auxCSatomA

call lstimer('START ',ts,te,lupri)

molecule => SETTING%MOLECULE(1)%p      ! --- prepare for re-initialization of SETTING
allocate(ATOMS(nAtoms))  ! --- each molecule made of only one atom

CALL pari_set_atomic_fragments(molecule,ATOMS,nAtoms,lupri)
call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)

call mem_alloc(Xalpha_ad,nAuxA,nRegA,nBastReg)
call ls_DZERO(Xalpha_ad,nAuxA*nRegA*nBastReg)

call mem_alloc(alpha_cd_full,nAuxA,nBastReg,nBastReg)
alpha_cd_full=0E0_realk
IF (setting%scheme%CS_SCREEN) THEN
  call ls_subScreenAtomic(auxCSatomA,auxCSfull,iAtomA,1,nAuxA,1,.FALSE.)
  call ls_attach_gab_to_setting(setting,auxCSatomA,regCSfull)
ENDIF
call pari_alphacd_full(alpha_cd_full,setting,molecule,atoms,iAtomA,nAuxA,nBastReg,nBastReg,lupri,luerr)
IF (setting%scheme%CS_SCREEN) THEN
  call ls_free_gab_from_setting(setting,lupri)
  call lstensor_free(auxCSatomA)
  deallocate(auxCSatomA)
ENDIF

! --- Loop over atoms C
DO iAtomC=1,nAtoms
  call getAtomicOrbitalInfo(orbitalInfo,iAtomC,nRegC,startRegC,endRegC,nAuxC,startAuxC,endAuxC)
  ! --- Loop over atoms D, starting from C since using (alpha|cd) symmetries
  call mem_alloc(Xalpha_ac,nAuxA,nRegA,nRegC)
  call ls_DZERO(Xalpha_ac,nAuxA*nRegA*nRegC)
!  DO iAtomD=iAtomC,nAtoms
  DO iAtomD=1,iAtomC
    call getAtomicOrbitalInfo(orbitalInfo,iAtomD,nRegD,startRegD,endRegD,nAuxD,startAuxD,endAuxD)
    ! --- Extract (alpha|cd) with alpha in A, c in C, d in D
    call mem_alloc(alpha_cd,nAuxA,nRegC,nRegD)
    call mem_alloc(alpha_dc,nAuxA,nRegD,nRegC)

    alpha_cd = alpha_cd_full(:,startRegC:startRegC+nRegC-1,startRegD:startRegD+nRegD-1)
!   alpha_cd = 0E0_realk
!   call pari_alphacd(alpha_cd,setting,atoms,iAtomA,iAtomC,iAtomD,nAuxA,nRegC,nRegD,lupri,luerr)
!   Make transposition of C and D

    call pariK_transpose(alpha_dc,alpha_cd,nAuxA,nRegC,nRegD)

    call mem_alloc(Xalpha_ca,nAuxA,nRegC,nRegA)
    ! --- Contract the density coeff. dalpha_ad and add contributions to K
    ! --- Construct the density-weighted integrals Xalpha_ac and Xalpha_ad
    
!   Generate K(A,C) += dalpha_ad*alpha_cd
    DO iRegD=1,nRegD
      call sub_dgemm('T','N',nRegA,nRegC,nAuxA,1E0_realk,dalpha_ad(1,1,startRegD+iRegD-1),nAuxA,&
   &             alpha_cd(1,1,iRegD),nAuxA,1E0_realk,Kfull(startRegA,startRegC,1,1,1),nBastReg)
    ENDDO ! --- end loop iRegD     

!   Generate Xalpha_ca
    call sub_DGEMM('N','T',nAuxA*nRegC,nRegA,nRegD,1E0_realk,alpha_cd,nAuxA*nRegC,&
   &           Dfull(startRegA,startRegD,1),nBastReg,0E0_realk,Xalpha_ca,nAuxA*nRegC)

!   Make transposition of Xalpha_ca to form Xalpha_ac                                                                                                         
    call parik_add_transpose(Xalpha_ac,Xalpha_ca,nAuxA,nRegA,nRegC)
    call mem_dealloc(Xalpha_ca)

    IF (iAtomC.NE.iAtomD) THEN
!    Make transpose dalpha_ca of dalpha_ac
     call mem_alloc(dalpha_ca,nAuxA,nRegC,nRegA)
     CALL parik_transpose(dalpha_ca,dalpha_ad(:,:,startRegC:startRegC+nRegC-1),nAuxA,nRegA,nRegC)

!    Contract K(A,D) = K(A,D) + dalpha_ca * alpha_cd
     call sub_dgemm('T','N',nRegA,nRegD,nAuxA*nRegC,1E0_realk,dalpha_ca,nAuxA*nRegC,&
   &            alpha_cd,nAuxA*nRegC,1E0_realk,Kfull(startRegA,startRegD,1,1,1),nBastReg)
!    
     call mem_dealloc(dalpha_ca)
!
!    call mem_alloc(Xalpha_da,nAuxA,nRegD,nRegA)
     DO iRegD=1,nRegD
      CALL sub_DGEMM('N','T',nAuxA,nRegA,nRegC,1E0_realk,alpha_cd(1,1,iRegD),nAuxA,&
   &             Dfull(startRegA,startRegC,1),nBastReg,1E0_realk,Xalpha_ad(1,1,iRegD+startRegD-1),nAuxA)
     ENDDO ! --- end loop iRegD
    END IF
!    call dgemm('T','N',nAux)
    call mem_dealloc(alpha_dc)
    call mem_dealloc(alpha_cd)
  ENDDO ! --- end loop D

  ! --- Loop over atoms B
  !> \todo (Patrick) loop should be later include only atoms in the vicinity of A
  DO iAtomB=1,nAtoms
   IF (iAtomA.NE.iAtomB) THEN
    call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
    ! --- Contract the "density-weighted integrals" Xalpha_ac with the PARI coeff., and add contributions to K
    call sub_DGEMM('T','N',nRegB,nRegC,nAuxA*nRegA,&
     &             1E0_realk,calpha_ab(iAtomA)%elements(1,1,startRegB),nAuxA*nRegA,&
     &             Xalpha_ac,nAuxA*nRegA,1E0_realk,Kfull(startRegB,startRegC,1,1,1),nBastReg)

   ENDIF
  ENDDO ! --- end loop B
    !
  call mem_dealloc(Xalpha_ac)

ENDDO ! --- end loop C  
DO iAtomB=1,nAtoms
 IF (iAtomA.NE.iAtomB) THEN
  call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
  ! --- Contract the "density-weighted integrals" Xalpha_ad with the PARI coeff., and add contributions to K
  call parik_xalpha_cont(Kfull,calpha_ab(iAtomA)%elements,Xalpha_ad,&
   &                     nAuxA,nRegA,nRegB,nBastReg,startRegB)
 ENDIF
ENDDO ! --- end loop B
call lstimer('3-loop ',ts,te,lupri)
call mem_dealloc(alpha_cd_full)
call mem_dealloc(Xalpha_ad)
CALL pari_free_atomic_fragments(ATOMS,nAtoms)
deallocate(ATOMS)  
! Return setting to default
call typedef_setMolecules(setting,molecule,1,2,3,4)
END SUBROUTINE pariK_ThreeCenter

SUBROUTINE pariK_TwoCenter(Kfull,Dfull,calpha_ab,dalpha_ad,alpha_beta,&
     &                 iAtomA,orbitalInfo,nBastReg,nAtoms,lupri)
implicit none
Real(realk),pointer                   :: Dfull(:,:,:)
Real(realk),pointer                   :: Kfull(:,:,:,:,:)
Integer,intent(in)                    :: lupri
Integer,intent(in)                    :: iAtomA, nAtoms,nBastReg
TYPE(MOLECULARORBITALINFO),intent(in) :: orbitalInfo
Real(realk),pointer                   :: dalpha_ad(:,:,:)
Real(realk),pointer                   :: alpha_beta(:,:,:,:,:)
TYPE(MAT3D),pointer                   :: calpha_ab(:)
!
Integer                    :: iAtomB,iAtomC,iAtomD
Real(realk),pointer        :: dalpha_ca(:,:,:)
Real(realk),pointer        :: extracted_alpha_cd_onC(:,:,:)
Real(realk),pointer        :: extracted_alpha_dc_onC(:,:,:)
Real(realk),pointer        :: extracted_alpha_gamma(:,:)
Real(realk),pointer        :: fXalpha_ac(:,:,:)
Real(realk),pointer        :: fXalpha_ad(:,:,:)
Real(realk),pointer        :: fXalpha_ca(:,:,:)
Real(realk),pointer        :: fXalpha_da(:,:,:)
Real(realk),pointer        :: elms(:),alphaCD(:),Dalpha(:)
Real(realk),pointer        :: dalpha_da(:,:,:)

Integer                    :: nRegA,nRegB,nRegC,nRegD,nAuxA,nAuxB,nAuxC,nAuxD
Integer                    :: startRegA,startRegB,startRegC,startRegD,startAuxA,startAuxB,startAuxC,startAuxD
Integer                    :: endRegA,endRegB,endRegC,endRegD,endAuxA,endAuxB,endAuxC,endAuxD
Integer                    :: iAlpha,iRegA,iRegB,iRegC,iRegD,iAuxA,iAuxC
!
TYPE(AtomSparseMat)        :: alphaBeta
TYPE(MoleculeInfo),pointer :: molecule
Real(realk)                :: factor,maxElm,tmp
Real(realk)                :: ts,te
!
!
call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)

call mem_alloc(fXalpha_ad,nAuxA,nRegA,nBastReg)
call mem_alloc(fXalpha_da,nAuxA,nBastReg,nRegA)
call mem_alloc(dalpha_da,nAuxA,nBastReg,nRegA)

fXalpha_da = 0E0_realk

! Make transposition of dalpha_ad
DO iRegD=1,nBastReg
  DO iRegA=1,nRegA
    DO iAuxA=1,nAuxA
      dalpha_da(iAuxA,iRegD,iRegA) = dalpha_ad(iAuxA,iRegA,iRegD)
    ENDDO
  ENDDO
ENDDO
  
! --- Loop over atoms C, starting from A since using (alpha|gamma) symmetries
DO iAtomC=1,iAtomA
  IF (iAtomC.EQ.iAtomA) THEN
    factor = -0.5E0_realk
  ELSE
    factor = -1E0_realk
  ENDIF
  call getAtomicOrbitalInfo(orbitalInfo,iAtomC,nRegC,startRegC,endRegC,nAuxC,startAuxC,endAuxC)

  call mem_alloc(extracted_alpha_gamma,nAuxA,nAuxC)

  ! --- Extract (alpha|gamma) with alpha in A, gamma in C
  extracted_alpha_gamma(1:nAuxA,1:nAuxC) = alpha_beta(startAuxA:endAuxA,1,startAuxC:endAuxC,1,1)

  call mem_alloc(extracted_alpha_cd_onC,nAuxA,nRegC,nBastReg)
  call mem_alloc(extracted_alpha_dc_onC,nAuxA,nBastReg,nRegC)

  ! Calculate the (alpha|fitted cd on c) = (alpha|gamma) c_gamma^cd
  call sub_DGEMM('N','N',nAuxA,nRegC*nBastReg,nAuxC,1E0_realk,extracted_alpha_gamma,nAuxA,&
  &              calpha_ab(iAtomC)%elements,nAuxC,0E0_realk,extracted_alpha_cd_onC,nAuxA)

! Make transposition (alpha|fitted dc on c)
  DO iRegD=1,nBastReg
    DO iRegC=1,nRegC
      DO iAuxA=1,nAuxA
        extracted_alpha_dc_onC(iAuxA,iRegD,iRegC) = extracted_alpha_cd_onC(iAuxA,iRegC,iRegD)
      ENDDO
    ENDDO
  ENDDO

! K(A,C) = K(A,C) - d_alpha^{ad} (alpha|cd)
  call sub_DGEMM('T','N',nRegA,nRegC,nAuxA*nBastReg,factor,dalpha_da,nAuxA*nBastReg,&
  &              extracted_alpha_dc_onC,nAuxA*nBastReg,1E0_realk,Kfull(startRegA,startRegC,1,1,1),nBastReg)

  call mem_alloc(fXalpha_ca,nAuxA,nRegC,nRegA)
! Xalpha_ca = (alpha|cd) D_da
  call sub_DGEMM('N','T',nAuxA*nRegC,nRegA,nBastReg,1E0_realk,extracted_alpha_cd_onC,nAuxA*nRegC,&
  &               Dfull(startRegA,1,1),nBastReg,0E0_realk,fXalpha_ca,nAuxA*nRegC)

  call mem_alloc(fXalpha_ac,nAuxA,nRegA,nRegC)
! Make Xalpha_ac by taking the transpose of Xalpha_ca with respect to a and c
  DO iRegC=1,nRegC
   DO iRegA=1,nRegA
    DO iAlpha=1,nAuxA
     fXalpha_ac(iAlpha,iRegA,iRegC) = fXalpha_ca(iAlpha,iRegC,iRegA)
    ENDDO ! --- end loop iAlpha
   ENDDO ! --- end loop iRegA
  ENDDO ! --- end loop iRegC
  call mem_dealloc(fXalpha_ca)

! Remove case iAtomC = iAtomD by zeroing out extracted_alpha_cd_onC and extracted_alpha_dc_onC
  extracted_alpha_cd_onC(:,:,startRegC:endRegC) = 0E0_realk
  extracted_alpha_dc_onC(:,startRegC:endRegC,:) = 0E0_realk

! Contract K(A,D) = K(A,D) + dalpha_ca * alpha_cd
  call sub_DGEMM('T','N',nRegA,nBastReg,nAuxA*nRegC,factor,dalpha_da(1,startRegC,1),nAuxA*nBastReg,&
   &             extracted_alpha_cd_onC,nAuxA*nRegC,1E0_realk,Kfull(startRegA,1,1,1,1),nBastReg)


! X(alpha,A,D) = X(alpha,A,D) + (alpha|CD)*D(A,C)
  call sub_DGEMM('N','T',nAuxA*nBastReg,nRegA,nRegC,factor,extracted_alpha_dc_onC,nAuxA*nBastReg,&
  &              Dfull(startRegA,startRegC,1),nBastReg,1E0_realk,fXalpha_da,nAuxA*nBastReg)

  call mem_dealloc(extracted_alpha_cd_onC)
  call mem_dealloc(extracted_alpha_dc_onC)

  DO iAtomB=1,nAtoms
   IF (iAtomB.NE.iAtomA) THEN
    call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
    ! --- Contract the "fitted density-weighted integrals" fXalpha_ac and fXalpha_ad with the PARI coeff., and add contributions to K

!   K(B,C) = K(B,C) - calpha_ab * Xalpha_ac
    call sub_DGEMM('T','N',nRegB,nRegC,nAuxA*nRegA,factor,calpha_ab(iAtomA)%elements(1,1,startRegB),nAuxA*nRegA,&
   &               fXalpha_ac,nAuxA*nRegA,1E0_realk,Kfull(startRegB,startRegC,1,1,1),nBastReg)

   END IF
  ENDDO ! --- end loop B	
  call mem_dealloc(fXalpha_ac)

  call mem_dealloc(extracted_alpha_gamma)   
ENDDO ! --- end loop C  

! Make transposition of fXalpha_da
DO iRegD=1,nBastReg
  DO iRegA=1,nRegA
    DO iAuxA=1,nAuxA
      fXalpha_ad(iAuxA,iRegA,iRegD) = fXalpha_da(iAuxA,iRegD,iRegA)
    ENDDO
  ENDDO
ENDDO
call mem_dealloc(fXalpha_da)

DO iAtomB=1,nAtoms
  IF (iAtomB.NE.iAtomA) THEN
    call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
    ! K(B,D) = K(B,D) - calpha_ab * Xalpha_ad
    call sub_DGEMM('T','N',nRegB,nBastReg,nAuxA*nRegA,1E0_realk,calpha_ab(iAtomA)%elements(1,1,startRegB),nAuxA*nRegA,&
     &             fXalpha_ad,nAuxA*nRegA,1E0_realk,Kfull(startRegB,1,1,1,1),nBastReg)
  END IF
ENDDO ! --- end loop B	
call mem_dealloc(fXalpha_ad)
call mem_dealloc(dalpha_da)
END SUBROUTINE pariK_TwoCenter

SUBROUTINE pari_alphaab(alpha_ab,setting,molecule,atoms,iAtomA,iAtomB,&
     &                  nAuxA,nAuxB,nRegA,nRegB,regCSfull,auxCSfull,lupri,luerr)
implicit none
TYPE(lssetting),intent(inout) :: setting
Integer,intent(in)            :: iAtomA,iAtomB,nAuxA,nAuxB,nRegA,nRegB,lupri,luerr
TYPE(moleculeinfo),pointer    :: atoms(:)
TYPE(moleculeinfo),intent(in) :: molecule
Real(realk),pointer           :: alpha_ab(:,:,:) !nAux,nRegA,nRegB
TYPE(LSTENSOR),pointer        :: regCSfull,auxCSfull
!
Logical                    :: sameAtoms
TYPE(moleculeinfo),pointer :: AB
TYPE(moleculeinfo),target  :: ABtarget
Integer                    :: nAux
TYPE(LSTENSOR),pointer     :: auxCSab,regCSab
INTEGER                    :: atomsAB(2),dummyAtoms(1)

!set threshold
SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%PARI_THRESHOLD

CALL pariSetPairFragment(AB,ABtarget,setting%basis(1)%p,molecule,atoms,molecule%nAtoms,&
     &                   iAtomA,iAtomB,nAuxA,nAuxB,nAux,lupri)

IF (setting%scheme%CS_SCREEN) THEN
  NULLIFY(auxCSab)
  ALLOCATE(auxCSab)
  NULLIFY(regCSab)
  ALLOCATE(regCSab)
  call ls_subScreenAtomic(regCSab,regCSfull,iAtomA,iAtomB,nRegA,nRegB,.FALSE.)
  IF (iAtomA.Eq.iAtomB) THEN
    call ls_subScreenAtomic(auxCSab,auxCSfull,iAtomA,1,nAuxA,1,.FALSE.)
  ELSE
    atomsAB(1)    = iAtomA
    atomsAB(2)    = iAtomB
    dummyAtoms(1) = 1
    call ls_subScreenFromList(auxCSab,auxCSfull,atomsAB,dummyAtoms,2,1,nAux,1,.FALSE.)
  ENDIF
ENDIF

call typedef_setMolecules(setting,AB,1,atoms(iAtomA),3,atoms(iAtomB),4)
call initIntegralOutputDims(setting%Output,nAux,1,nRegA,nRegB,1)
IF (setting%scheme%CS_SCREEN) call ls_attach_gab_to_setting(setting,auxCSab,regCSab)
call ls_getIntegrals(AODFdefault,AOEmpty,AORdefault,AORdefault,CoulombOperator,RegularSpec,&
     &               Contractedinttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,alpha_ab,.FALSE.)
IF (setting%scheme%CS_SCREEN) THEN
  call ls_free_gab_from_setting(setting,lupri)
  call lstensor_free(regCSab)
  call lstensor_free(auxCSab)
  DEALLOCATE(auxCSab)
  DEALLOCATE(regCSab)
ENDIF

CALL pariFreePairFragment(AB,iAtomA,iAtomB)

END SUBROUTINE pari_alphaab

! SUBROUTINE pari_alphacd(alpha_cd,setting,atoms,iAtomA,iAtomC,iAtomD,&
!      &                  nAuxA,nRegC,nRegD,lupri,luerr)
! implicit none
! TYPE(lssetting),intent(inout) :: setting
! Integer,intent(in)            :: iAtomA,iAtomC,iAtomD,nAuxA,nRegC,nRegD,lupri,luerr
! TYPE(moleculeinfo),pointer    :: atoms(:)
! Real(realk),intent(in)        :: alpha_cd(nAuxA,nRegC,nRegD)

! !set threshold
! SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
! call typedef_setMolecules(setting,atoms(iAtomA),1,atoms(iAtomC),3,atoms(iAtomD),4)
! call initIntegralOutputDims(setting%Output,nAuxA,1,nRegC,nRegD,1)
! call ls_getIntegrals(AODFdefault,AOEmpty,AORdefault,AORdefault,CoulombOperator,RegularSpec,&
!      &               Contractedinttype,SETTING,LUPRI,LUERR)
! CALL retrieve_Output(lupri,setting,alpha_cd,.FALSE.)
! END SUBROUTINE pari_alphacd

SUBROUTINE pari_alphacd_full(alpha_cd,setting,molecule,atoms,iAtomA,&
     &                       nAuxA,nRegC,nRegD,lupri,luerr)
implicit none
TYPE(lssetting),intent(inout) :: setting
Integer,intent(in)            :: iAtomA,nAuxA,nRegC,nRegD,lupri,luerr
TYPE(moleculeinfo),intent(in) :: molecule
TYPE(moleculeinfo),pointer    :: atoms(:)
Real(realk),intent(inout)     :: alpha_cd(nAuxA,nRegC,nRegD)

!set threshold
SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
call typedef_setMolecules(setting,molecule,3,4,atoms(iAtomA),1)
call initIntegralOutputDims(setting%Output,nAuxA,1,nRegC,nRegD,1)
call ls_getIntegrals(AODFdefault,AOEmpty,AORdefault,AORdefault,CoulombOperator,RegularSpec,&
     &               Contractedinttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,alpha_cd,.FALSE.)

END SUBROUTINE pari_alphacd_full

SUBROUTINE pari_alphaBeta(alphaBeta,setting,molecule,atoms,iAtomA,iAtomB,&
     &                    nAuxA,nAuxB,regCSfull,auxCSfull,lupri,luerr)
implicit none
TYPE(lssetting),intent(inout) :: setting
Integer,intent(in)            :: iAtomA,iAtomB,nAuxA,nAuxB,lupri,luerr
TYPE(moleculeinfo),pointer    :: atoms(:)
TYPE(moleculeinfo),intent(in) :: molecule
Real(realk),pointer           :: alphaBeta(:,:) !nAux,nAux
TYPE(LSTENSOR),pointer        :: regCSfull,auxCSfull

TYPE(moleculeinfo),pointer :: AB
TYPE(moleculeinfo),target  :: ABtarget
Integer                    :: nAux

TYPE(LSTENSOR),pointer     :: auxCSab
INTEGER                    :: atomsAB(2),dummyAtoms(1)
!set threshold
SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%PARI_THRESHOLD

CALL pariSetPairFragment(AB,ABtarget,setting%basis(1)%p,molecule,atoms,molecule%nAtoms,&
     &                   iAtomA,iAtomB,nAuxA,nAuxB,nAux,lupri)

IF (setting%scheme%CS_SCREEN) THEN
  NULLIFY(auxCSab)
  ALLOCATE(auxCSab)
  IF (iAtomA.EQ.iAtomB) THEN
    call ls_subScreenAtomic(auxCSab,auxCSfull,iAtomA,1,nAuxA,1,.FALSE.)
  ELSE
    atomsAB(1)    = iAtomA
    atomsAB(2)    = iAtomB
    dummyAtoms(1) = 1
    call ls_subScreenFromList(auxCSab,auxCSfull,atomsAB,dummyAtoms,2,1,nAux,1,.FALSE.)
  ENDIF
ENDIF

call typedef_setMolecules(setting,AB,1,3)
call initIntegralOutputDims(setting%Output,nAux,1,nAux,1,1)
IF (setting%scheme%CS_SCREEN) call ls_attach_gab_to_setting(setting,auxCSab,auxCSab)
call ls_getIntegrals(AODFdefault,AOEmpty,AODFdefault,AOEmpty,CoulombOperator,RegularSpec,&
     &               Contractedinttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,alphaBeta,.FALSE.)
IF (setting%scheme%CS_SCREEN) THEN
  call ls_free_gab_from_setting(setting,lupri)
  call lstensor_free(auxCSab)
  DEALLOCATE(auxCSab)
ENDIF
CALL pariFreePairFragment(AB,iAtomA,iAtomB)

END SUBROUTINE pari_alphaBeta

SUBROUTINE pari_charge_ab(charge_ab,setting,atomA,atomB,nRegA,nRegB,lupri,luerr)
implicit none
real(realk),intent(inout)       :: charge_ab(nRegA,nRegB)
TYPE(lssetting),intent(inout) :: setting
TYPE(moleculeinfo),intent(in) :: atomA,atomB
integer, intent(in)           :: nRegA,nRegB,lupri,luerr

!set threshold
SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%PARI_THRESHOLD
charge_ab = 0E0_realk
call typedef_setMolecules(setting,atomA,1,atomB,2)
call initIntegralOutputDims(setting%Output,nRegA,nRegB,1,1,1)
!call putScreeningSubblockToSetting(setting,screeningMatrices)
call ls_getIntegrals(AORdefault,AORdefault,AOEmpty,AOEmpty,OverlapOperator,RegularSpec,&
     &               Contractedinttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,charge_ab,.FALSE.)


END SUBROUTINE pari_charge_ab


SUBROUTINE pari_dipole_ab(dipole_ab,setting,atomA,atomB,nRegA,nRegB,lupri,luerr)
implicit none
real(realk),intent(inout)     :: dipole_ab(nRegA,nRegB,4)
TYPE(lssetting),intent(inout) :: setting
TYPE(moleculeinfo),intent(in) :: atomA,atomB
integer, intent(in)           :: nRegA,nRegB,lupri,luerr
integer                       :: nsphmat,nderiv,nbast
!
!set threshold
SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%PARI_THRESHOLD
call typedef_setMolecules(setting,atomA,1,atomB,2)
dipole_ab = 0E0_realk
nsphmat = 4 ! charge + {x,y,z} components of the dipole
call initIntegralOutputDims(setting%output,nRegA,nRegB,1,1,nsphmat)
setting%scheme%CMORDER = 1 ! First derivative to get up to the dipole moment
call ls_getIntegrals(AORdefault,AORdefault,AOEmpty,AOEmpty,CarmomOperator,RegularSpec,Contractedinttype,&
     &                SETTING,LUPRI,LUERR)
setting%scheme%CMORDER = 0
call retrieve_Output(lupri,setting,dipole_ab,.FALSE.)
END SUBROUTINE pari_dipole_ab

SUBROUTINE pari_alpha_charge(alpha_charge,setting,molecule,atoms,iAtomA,iAtomB,nAuxA,nAuxB,&
     &                       lupri,luerr)
implicit none
real(realk),intent(inout)       :: alpha_charge(:)
TYPE(lssetting),intent(inout) :: setting
TYPE(moleculeinfo),intent(in) :: molecule
TYPE(moleculeinfo),pointer    :: atoms(:)
integer, intent(in)           :: nAuxA,nAuxB,lupri,luerr,iAtomA,iAtomB

TYPE(moleculeinfo),pointer :: AB
TYPE(moleculeinfo),target  :: ABtarget
Integer                    :: nAux

!set threshold
SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%PARI_THRESHOLD
alpha_charge = 0E0_realk
CALL pariSetPairFragment(AB,ABtarget,setting%basis(1)%p,molecule,atoms,molecule%nAtoms,&
     &                   iAtomA,iAtomB,nAuxA,nAuxB,nAux,lupri)
call typedef_setMolecules(setting,AB,1)
call initIntegralOutputDims(setting%Output,nAux,1,1,1,1)
!call putScreeningSubblockToSetting(setting,screeningMatrices)
call ls_getIntegrals(AODFdefault,AOEmpty,AOEmpty,AOEmpty,OverlapOperator,RegularSpec,&
     &               Contractedinttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,alpha_charge,.FALSE.)
CALL pariFreePairFragment(AB,iAtomA,iAtomB)

END SUBROUTINE pari_alpha_charge


SUBROUTINE pari_alpha_dipole(alpha_dipole,setting,molecule,atoms,iAtomA,iAtomB,nAuxA,nAuxB,&
     &                       lupri,luerr)
implicit none
real(realk),intent(inout)       :: alpha_dipole(:,:)
TYPE(lssetting),intent(inout) :: setting
TYPE(moleculeinfo),intent(in) :: molecule
TYPE(moleculeinfo),pointer    :: atoms(:)
integer, intent(in)           :: nAuxA,nAuxB,lupri,luerr,iAtomA,iAtomB
!
TYPE(moleculeinfo),pointer    :: AB
TYPE(moleculeinfo),target     :: ABtarget
Integer                       :: nAux, nsphmat, nderiv
!
Integer :: i,j
!
!set threshold
SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%PARI_THRESHOLD
alpha_dipole = 0E0_realk
nderiv  = 1 ! First derivative to get up to the dipole moment
nsphmat = 4 ! charge + {x,y,z} components of the dipole
call pariSetPairFragment(AB,ABtarget,setting%basis(1)%p,molecule,atoms,molecule%nAtoms,&
     &                   iAtomA,iAtomB,nAuxA,nAuxB,nAux,lupri)
call typedef_setMolecules(setting,AB,1)
call initIntegralOutputDims(setting%output,nAux,1,1,1,nsphmat)
setting%scheme%cmorder = nderiv
call ls_getIntegrals(AODFdefault,AOEmpty,AOEmpty,AOEmpty,CarmomOperator,RegularSpec,Contractedinttype,&
     &                SETTING,LUPRI,LUERR)
setting%scheme%cmorder = 0
call retrieve_Output(lupri,setting,alpha_dipole,.FALSE.)

 call pariFreePairFragment(AB,iAtomA,iAtomB)
END SUBROUTINE pari_alpha_dipole

SUBROUTINE pariK_transpose(alpha_dc,alpha_cd,nAuxA,nRegC,nRegD)
implicit none
Integer, intent(in)     :: nAuxA,nRegC,nRegD
Real(realk),intent(IN)  :: alpha_cd(nAuxA,nRegC,nRegD)
Real(realk),intent(INOUT) :: alpha_dc(nAuxA,nRegD,nRegC)
!
Integer :: iRegC,iRegD
!
DO iRegC=1,nRegC
  DO iRegD=1,nRegD
    alpha_dc(:,iRegD,iRegC) = alpha_cd(:,iRegC,iRegD)
  ENDDO
ENDDO
END SUBROUTINE pariK_transpose

SUBROUTINE parik_add_transpose(Xalpha_ac,Xalpha_ca,nAuxA,nRegA,nRegC)
implicit none
Integer, intent(in)     :: nAuxA,nRegA,nRegC
Real(realk),intent(IN)  :: Xalpha_ca(nAuxA,nRegC,nRegA)
Real(realk),intent(INOUT) :: Xalpha_ac(nAuxA,nRegA,nRegC)
!
Integer :: iRegA,iRegC,iAlpha
!
DO iRegC=1,nRegC
  Xalpha_ac(:,:,iRegC) = Xalpha_ac(:,:,iRegC) + Xalpha_ca(:,iRegC,:)
ENDDO
!O iRegC=1,nRegC
!DO iRegA=1,nRegA
! DO iAlpha=1,nAuxA
!  Xalpha_ac(iAlpha,iRegA,iRegC) = Xalpha_ac(iAlpha,iRegA,iRegC) + Xalpha_ca(iAlpha,iRegC,iRegA)
! ENDDO
!ENDDO
!ENDDO
END SUBROUTINE parik_add_transpose

SUBROUTINE pari_set_atomic_fragments(molecule,ATOMS,nAtoms,lupri)
implicit none
Integer,intent(in)         :: lupri
Integer,intent(in)         :: nAtoms
TYPE(MoleculeInfo),pointer :: molecule
TYPE(MOLECULEINFO),pointer :: ATOMS(:)
Integer                    :: iAtom
Character(len=22)          :: FRAGMENTNAME
DO iAtom = 1,nAtoms
  IF (iAtom.GT. 999999999) CALL LSQUIT('Error in pari_set_atomic_fragments -> FRAGMENTNAME',-1)
  write(FRAGMENTNAME,'(A4,I18)') 'Atom',iAtom
  call init_MoleculeInfo(ATOMS(iAtom),1,FRAGMENTNAME)
  call copy_atom(molecule, iAtom, ATOMS(iAtom),1,LUPRI)
ENDDO
END SUBROUTINE pari_set_atomic_fragments

SUBROUTINE pari_free_atomic_fragments(ATOMS,nAtoms)
implicit none
Integer,intent(in)         :: nAtoms
TYPE(MOLECULEINFO),pointer :: ATOMS(:)
!
Integer                    :: iAtom
!
DO iAtom = 1,nAtoms 
  call free_MoleculeInfo(ATOMS(iAtom))
ENDDO
END SUBROUTINE pari_free_atomic_fragments

SUBROUTINE parik_xalpha_cont(Kfull,calpha_ab,Xalpha_ad,nAuxA,nRegA,nRegB,nBastReg,startRegB)
implicit none
Integer,intent(IN) :: nAuxA,nRegA,nRegB,nBastReg,startRegB
Real(realk),intent(INOUT) :: Kfull(nBastReg,nBastReg,1,1,1)
Real(realk),intent(IN)    :: calpha_ab(nAuxA,nRegA,nBastReg)
Real(realk),intent(IN)    :: Xalpha_ad(nAuxA,nRegA,nBastReg)
!
Integer :: dimA
Integer :: iAlpha,iRegA,iRegB,iRegD
Real(realk) :: tmp

dimA = nAuxA*nRegA

call sub_DGEMM('T','N',nRegB,nBastReg,dimA,1E0_realk,calpha_ab(1,1,startRegB),dimA,&
     &         Xalpha_ad,dimA,1E0_realk,Kfull(startRegB,1,1,1,1),nBastReg)
!call mem_alloc(Kred,nRegB,nBastReg)
!DO iRegD=1,nBastReg
!  DO iRegB=1,nRegB
!    call sub_DGEMM()
!    tmp = 0E0_realk
!    DO iRegA=1,nRegA
!      DO iAlpha=1,nAuxA
!        tmp = tmp + calpha_ab(iAlpha,iRegA,iRegB+startRegB-1)*Xalpha_ad(iAlpha,iRegA,iRegD)
!        Kfull(iRegB+startRegB-1,iRegD,1,1,1)=Kfull(iRegB+startRegB-1,iRegD,1,1,1)&
!      & + calpha_ab(iAlpha,iRegA,iRegB+startRegB-1)*Xalpha_ad(iAlpha,iRegA,iRegD)
!      ENDDO ! --- end loop iRegC
!    ENDDO ! --- end loop iAlpha
!  ENDDO ! --- end loop iRegA
!ENDDO ! --- end loop iRegB
END SUBROUTINE parik_xalpha_cont

SUBROUTINE pariSetPairFragment(AB,ABtarget,basis,molecule,atoms,nAtoms,iAtomA,iAtomB,&
     &                         nAuxA,nAuxB,nAux,lupri)
TYPE(basisinfo),intent(inout) :: basis
TYPE(moleculeinfo),intent(in) :: molecule
TYPE(moleculeinfo),pointer    :: AB
TYPE(moleculeinfo),target     :: ABtarget
Integer,intent(in)            :: nAtoms,iAtomA,iAtomB,nAuxA,nAuxB,lupri
Integer,intent(inout)           :: nAux
TYPE(moleculeinfo),target     :: atoms(nAtoms)

Logical                    :: sameAtoms
Integer                    :: iAtoms(2)
Character(len=22)          :: FRAGMENTNAME

sameAtoms = iAtomA.EQ.iAtomB
IF (sameAtoms) THEN
  AB => atoms(iAtomA)
  nAux = nAuxA
ELSE
  AB => ABtarget
  IF ((iAtomA.GT. 999999999).OR.(iAtomB.GT. 999999999)) CALL LSQUIT('Error in pari_alphaab -> FRAGMENTNAME',-1)
  write(FRAGMENTNAME,'(A4,2I9)') 'Atom',iAtomA,iAtomB
!  call init_MoleculeInfo(AB,2,FRAGMENTNAME)
  iAtoms(1) = iAtomA
  iAtoms(2) = iAtomB
  call build_fragment(molecule,AB,basis,.TRUE.,.FALSE.,.FALSE.,iAtoms,2,lupri)
  nAux = nAuxA + nAuxB
ENDIF
END SUBROUTINE pariSetPairFragment

SUBROUTINE pariFreePairFragment(AB,iAtomA,iAtomB)
implicit none
TYPE(moleculeinfo),pointer    :: AB
Integer,intent(IN)            :: iAtomA,iAtomB
IF (iAtomA.NE.iAtomB) THEN
  call free_MoleculeInfo(AB)
ENDIF
END SUBROUTINE pariFreePairFragment

!> \brief Calculate eigenvalues(W)/eigenvectors(A) of a square (symmetric or not) matrix A with rank rankA
!> \author P. Merlot
!> \date 2012-09-17
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param A The square matrix to diagonalize as input, its eigenvectors in output
!> \param W The eigenvalues of matrix A
!> \param rankA The rank of the square matrix A
SUBROUTINE II_get_eigv_square_mat(LUPRI,LUERR,A,W,nbast)
  implicit none
  INTEGER,INTENT(IN)        :: LUPRI,LUERR,nbast
  REAL(REALK),INTENT(INOUT) :: A(nbast,nbast)
  REAL(REALK),INTENT(INOUT) :: W(nbast) ! eigenvalues
  !
  REAL(REALK), pointer  :: WORK(:)
  INTEGER               :: INFO,LWORK,LWMAX
  REAL(realk) :: ts,te
  CHARACTER             :: JOBVS, SORT
  INTEGER               :: LDA, LDVS, N, SDIM, nb,  nbZeroEig, nbComplexEig, iReg
  LOGICAL,pointer       :: BWORK(:)
  REAL(REALK),pointer   :: VS(:,:), WI(:)
  INTEGER               :: iNegImagEigV
  REAL(REALK),pointer   :: negImagEigVal(:,:)
  INFO=0
  !
  !--------------------- DSYEV METHOD -------------------------
  !--------------------- Symmetric matrix A -------------------
#if 1
  LWORK = -1 ! returns the optimal value of LWORK in WORK(1)
  INFO = 0
  call mem_alloc(WORK,2)
  ! write(*,*) '   --- Query the optimal workspace.'
  CALL DSYEV( 'V', 'U', nbast, A,nbast, W, WORK, LWORK, INFO )
  LWORK = NINT(WORK(1))
  call mem_dealloc(WORK)
  call mem_alloc(WORK,LWORK)
  !  write(*,*) '   --- Solve eigenproblem'
  ! DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
  CALL DSYEV( 'V', 'U', nbast, A,nbast, W, WORK, LWORK, INFO )
  ! Check for convergence.                                             
  IF( INFO.GT.0 ) THEN
     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
     STOP
  END IF
  !  write(*,*) 'Eigenvalues:'
  !  write(*,*) W
  call mem_dealloc(WORK)
  !--------------------- DSYEV END -------------------------

  !--------------------- DGEES METHOD -------------------------
  !--------------------- non-symmetric matrix A ---------------
#else
  call mem_alloc(BWORK,nbast)
  call mem_alloc(WI,nbast)
  call mem_alloc(VS,nbast,nbast)
  call ls_dzero(WI,nbast)
  call ls_dzero(W,nbast)
  nb = 64
  LWORK = max((2+nb)*nbast,999998 )
  call mem_alloc(WORK,LWORK)
  SDIM = 0
  call DGEES( 'N', 'N', SELECTEIGV, nbast, A, nbast, SDIM, W, WI,VS, nbast, WORK, LWORK, BWORK, INFO )

  nbZeroEig = 0
  nbComplexEig = 0 
  DO iReg=1,nbast
     IF (abs(WI(iReg)) .GT. 1E-10_realk) THEN
        nbComplexEig = nbComplexEig +1
     ELSE
        nbZeroEig = nbZeroEig + 1
     ENDIF
  ENDDO

  IF (nbComplexEig .GT. 0) THEN
     write (*,*) '  ---> Nb. complex eigenvalues:',nbComplexEig , ' out of ',nbast
  ENDIF
  IF (nbZeroEig .GT. 0) THEN
     write (*,*) '  ---> Nb. real eigenvalues:    ',nbZeroEig, ' out of ',nbast
  ENDIF

  ! print complex eigenvalues if imaginary part non-zero
  IF (nbComplexEig .GT. 0) THEN
     iNegImagEigV = 0
     call mem_alloc(negImagEigVal,nbComplexEig,2)
     DO iReg=1,nbast
        IF (abs(WI(iReg)) .GT. 1E-10_realk) THEN
           iNegImagEigV = iNegImagEigV +1
           negImagEigVal(iNegImagEigV,1) = W(ireg)
           negImagEigVal(iNegImagEigV,2) = WI(ireg)
        ENDIF
     ENDDO
          write (*,*) "The negative COMPLEX eigenvalues:"
          DO ireg=1,nbComplexEig
             write(*,*) "(",negImagEigVal(ireg,1),";",negImagEigVal(ireg,2),")"
          ENDDO
     call mem_dealloc(negImagEigVal)
  ENDIF
  IF( INFO.NE.0 ) THEN
     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
     STOP
  END IF
  !write(*,*) 'Eigenvalues:'
  call mem_dealloc(WORK)
  call mem_dealloc(BWORK)
  call mem_dealloc(WI)
  call mem_dealloc(VS)
  !--------------------- DGEES END -------------------------
#endif
  !  CALL lstimer('get_eigv_mat',ts,te,6)

END SUBROUTINE II_get_eigv_square_mat

Logical FUNCTION selecteigv(ar, ai)
  !     .. Scalar Arguments ..
  !     Logical function SELECT for use with DGEES
  !     Returns the value .TRUE. if the imaginary part of the eigenvalue
  !     (AR + AI*i) is zero, i.e. the eigenvalue is real
  real(realk) :: ar
  real(realk) :: ai
  !     .. Local Scalars ..
  Logical :: d
  !     .. Executable Statements ..
  IF (ai .LT. -1E-10_realk) THEN
     d = .TRUE.
  ELSE
     d = .FALSE.
  END IF
  selecteigv = d
  RETURN
END FUNCTION selecteigv


!> \brief Calculate the integrals (alpha|beta) with aux. functions over the whole molecule
!> \author P. Merlot
!> \date 2012-09-17
!> \param alphaBeta The matrix of 2-center 2-electron integrals
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Default error print unit
SUBROUTINE pari_alphaBetaFull(alphaBetaFull,setting,molecule,lupri,luerr)
  implicit none
  Real(realk),pointer           :: alphaBetaFull(:,:)
  TYPE(lssetting),intent(inout) :: setting
  TYPE(moleculeinfo),intent(in) :: molecule
  Integer,intent(in)            :: lupri,luerr
  TYPE(LSTENSOR),pointer        :: auxCSfull
  !
  Integer                    :: nAux
  TYPE(LSTENSOR),pointer     :: auxCSab
  !
  nAux = molecule%nbastAux
  !set threshold
  SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%PARI_THRESHOLD
  call typedef_setMolecules(setting,molecule,1,3)
  call initIntegralOutputDims(setting%Output,nAux,1,nAux,1,1)
  call ls_getIntegrals(AODFdefault,AOEmpty,AODFdefault,AOEmpty,CoulombOperator,RegularSpec,&
       &               Contractedinttype,SETTING,LUPRI,LUERR)
  call retrieve_Output(lupri,setting,alphaBetaFull,.FALSE.)
END SUBROUTINE pari_alphaBetaFull

!> \brief Identifies the max, min values of an array, and their ratio (condition number), typically used for eigenvalues
!> \author P. Merlot
!> \date 2012-09-17
!> \param eigVal The array to investigate
!> \param nb_eigVal The length of the array
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param minEigv The minimum number of the array
!> \param maxEigv The maximum number of the array
!> \param conditionNum The ratio abs(maxEigv)/abs(minEigv)
SUBROUTINE check_min_max_Array_elem(minEigv,maxEigv,conditionNum,eigVal,nb_eigVal,lupri,luerr)
  IMPLICIT NONE
  INTEGER,INTENT(IN)        :: lupri,luerr,nb_eigVal
  Real(realk),INTENT(IN)    :: eigval(nb_eigVal)
  Real(realk),INTENT(INOUT) :: minEigv,maxEigv,conditionNum
  !                                    
  INTEGER                   :: i
  minEigV = eigval(1)
  maxEigV = eigval(1)
  conditionNum = abs(maxEigV)/abs(minEigV)
  DO i=1,nb_eigVal
     minEigV = min(minEigV,eigval(i))
     maxEigV = max(maxEigV,eigval(i))
  ENDDO
  conditionNum = abs(maxEigV)/abs(minEigV)
  write(lupri,*) "minEigV= ",minEigv
  write(*,*) "minEigV= ",minEigv
  write(lupri,*) "maxEigV= ",maxEigv
  write(*,*) "maxEigV= ",maxEigv
  write(lupri,*) "condition nb.= ",conditionNum
  write(*,*) "condition nb.= ",conditionNum
END SUBROUTINE check_min_max_Array_elem

END MODULE pari_mod

SUBROUTINE sub_DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
use memory_handling
use precision
implicit none
CHARACTER * 1  :: TRANSA, TRANSB
INTEGER        :: M, N, K, LDA, LDB, LDC
Real(realk)    :: ALPHA, BETA, Amax, Bmax
Real(realk)    :: A(LDA,*), B(LDB,*), C(LDC,*)
Real(realk)    :: maxkA(M), maxkB(N), maxiA(K),maxjB(K)
Real(realk),pointer :: Ared(:,:),Bred(:,:),Cred(:,:)

INTEGER        :: indexM(M), indexN(N), indexK(K)
INTEGER        :: i, j, kk
INTEGER        :: Mred, Nred, Kred, iMred, iNred, iKred,LDAred,LDBred,LDCred
INTEGER        :: iM, iN, iK
Real(realk)    :: threshold,red_fac,tmp
Logical        :: DO_REGULAR

#if 0
CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)

#else
! Set up an arbitrary threshold
!ToDo Make both as input to sub_DGEMM
threshold = 1.0E-15_realk
red_fac   = 0.8E0_realk


! find Amax,Bmax: the max element of A and B in absolute value
! find maxkA: the M-dim column vector of max. values of A for each row i
! find maxkB: the N-dim row    vector of max. values of B for each column j
! find maxiA: the K-dim column vector of max. values of A for each column kk
! find maxjB: the K-dim column vector of max. values of B for each column kk
Amax  = 0E0_realk
maxkA = 0E0_realk
Bmax  = 0E0_realk
maxkB = 0E0_realk
maxiA = 0E0_realk
maxjB = 0E0_realk
IF (TRANSA.EQ.'T') THEN
   DO i=1,M
      DO kk=1,K
         tmp = ABS(A(kk,i))
         maxkA(i)  = MAX(tmp,maxkA(i))  
         maxiA(kk) = MAX(tmp,maxiA(kk)) 
         Amax = MAX(maxkA(i),Amax)
      ENDDO
   ENDDO
ELSE
   DO kk=1,K
      DO i=1,M
         tmp = ABS(A(i,kk))
         maxkA(i)  = MAX(tmp,maxkA(i))  
         maxiA(kk) = MAX(tmp,maxiA(kk))
         Amax = MAX(maxkA(i),Amax)
      ENDDO
   ENDDO
ENDIF

IF (TRANSB.EQ.'T') THEN
   DO kk=1,K
      DO j=1,N
         tmp = ABS(B(j,kk))
         maxkB(j)  = MAX(tmp,maxkB(j))
         maxjB(kk) = MAX(tmp,maxjB(kk))
         Bmax = MAX(maxkB(j),Bmax)
      ENDDO
   ENDDO
ELSE
   DO j=1,N
      DO kk=1,K
         tmp = ABS(B(kk,j))
         maxkB(j)  = MAX(tmp,maxkB(j))
         maxjB(kk) = MAX(tmp,maxjB(kk))
         Bmax = MAX(maxkB(j),Bmax)
      ENDDO
   ENDDO
ENDIF

! Which rows of A can we keep?
Mred = 0
indexM   = 0
DO i=1,M
   IF (maxkA(i)*Bmax .GT. threshold) THEN
      Mred = Mred + 1
      indexM(Mred)=i
   ENDIF
ENDDO

! Which columns of B can we remove?
Nred = 0
indexN = 0
DO j=1,N
   IF (Amax*maxkB(j) .GT. threshold) THEN
      Nred = Nred +1
      indexN(Nred) = j
   ENDIF
ENDDO

! Which rows (and columns) of B (and A) can we remove?
Kred = 0
indexK = 0

DO kk = 1,K
   IF (maxiA(kk)*maxjB(kk) .GT. threshold) THEN
      Kred = Kred + 1
      indexK(Kred) = kk
   ENDIF
ENDDO

!CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)

#if 1
! Perform DGEMM only when there is a signficant contribution
IF ((NRED.GT. 0).AND.(MRED.GT. 0).AND.(KRED.GT. 0)) THEN
  
 DO_REGULAR = (1E0_realk*MRED/M.GT.red_fac).AND.(1E0_realk*NRED/N.GT.red_fac).AND.(1E0_realk*KRED/K.GT.red_fac)
 !If reduction is insignificant do it the traditional way
 IF (DO_REGULAR) THEN
  CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
 !otherwise do a reduction to save time
 ELSE
  ! Create the reduced matrices
  IF (TRANSA.EQ.'T') THEN
     CALL mem_alloc(Ared,Kred,Mred)
     LDAred = Kred
     DO iM = 1,Mred
        iMred = indexM(iM)
        DO iK = 1,Kred
           iKred = indexK(iK)
           Ared(iK,iM) = A(iKred,iMred)
        ENDDO
     ENDDO
  ELSE
     CALL mem_alloc(Ared,Mred,Kred)
     LDAred = Mred
     DO iK = 1,Kred
        iKred = indexK(iK)
        DO iM = 1,Mred
           iMred = indexM(iM)
           Ared(iM,iK) = A(iMred,iKred)
        ENDDO
     ENDDO
  ENDIF
  
  IF (TRANSB.EQ.'T') THEN
     CALL mem_alloc(Bred,Nred,Kred)
     LDBred = Nred
     DO iK = 1,Kred
        iKred = indexK(iK)
        DO iN = 1,Nred
           iNred = indexN(iN)
           Bred(iN,iK) = B(iNred,iKred)
        ENDDO
     ENDDO
  ELSE
     CALL mem_alloc(Bred,Kred,Nred)
     LDBred = Kred
     DO iN = 1,Nred
        iNred = indexN(iN)
        DO iK = 1,Kred
           iKred = indexK(iK)
           Bred(iK,iN) = B(iKred,iNred)
        ENDDO
     ENDDO
  ENDIF
  
  ! Compute the multiplication of the reduced matrices Cred = Ared * Bred
  CALL mem_alloc(Cred,Mred,Nred)
  LDCred = Mred
  CALL DGEMM(TRANSA, TRANSB, Mred, Nred, Kred, ALPHA, Ared, LDAred, Bred, LDBred, 0E0_realk, Cred, LDCred)
  
  
  
  IF (beta.EQ. 0E0_realk) THEN
  ! Zeros in the resulting matrix C
    C(1:M,1:N) = 0E0_realk
    DO iN = 1,Nred
       iNred = indexN(iN)
       DO iM = 1,Mred
          iMred = indexM(iM)
          C(iMred,iNred) = Cred(iM,iN)
       ENDDO
    ENDDO
  ELSE
    ! Map back the reduced Cred matrix to the full matrix C
    DO iN = 1,Nred
       iNred = indexN(iN)
       DO iM = 1,Mred
          iMred = indexM(iM)
          C(iMred,iNred) = C(iMred,iNred) + beta*Cred(iM,iN)
       ENDDO
    ENDDO
  ENDIF
  
  
  CALL mem_dealloc(Ared)
  CALL mem_dealloc(Bred)
  CALL mem_dealloc(Cred)
 ENDIF
  
ELSE
  IF (BETA.EQ. 0E0_realk) THEN
    C(1:M,1:N) = 0E0_realk
  ELSE IF (BETA.NE. 1E0_realk) THEN
    C(1:M,1:N) = BETA*C(1:M,1:N)
  ENDIF
ENDIF

#endif
#endif

END SUBROUTINE sub_DGEMM
