!> @file
!> Full molecular calculation of RI-MP2-F12

module fullrimp2f12 

  use precision
  use typedeftype!,only:lsitem
  use typedef
  use dec_typedef_module
  use matrix_module
  use matrix_operations
  use memory_handling
  use dec_fragment_utils
  use CABS_operations


  !WARNING FOR TESTING
  use f12_routines_module


  public :: full_canonical_rimp2_f12

  private

contains
  !> \brief Calculate canonical MP2 energy for full molecular system
  !> keeping full AO integrals in memory. Only for testing.
  !> \author Thomas Kjaergaard
  !> \date 2015
  subroutine full_canonical_rimp2_f12(MyMolecule,MyLsitem,Dmat,mp2f12_energy)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: Dmat
    !> MP2-F12 correlation energy
    real(realk),intent(inout) :: mp2f12_energy
    !> Canonical MP2 correlation energy
    real(realk) :: mp2_energy
    !local variables
    integer :: nbasis,nocc,nvirt,ncabsAO,ncabs,noccfull
    real :: E21
    !========================================================
    !WARNING THESE SHOULD NOT BE HERE - ONLY FOR TESTING
    !========================================================
    Real(realk),pointer :: Vijij(:,:)
    Real(realk)  :: Ripjq(:,:,:,:) !nocc,nbasis,nocc,nbasis
    Real(realk)  :: Gipjq(:,:,:,:) !nocc,nbasis,nocc,nbasis
    Real(realk)  :: Fijkl(:,:,:,:) !nocc,nocc,nocc,nocc
    Real(realk)  :: Fijkl(:,:,:,:) !nocc,nocc,nocc,nocc
    Real(realk)  :: gao(:,:,:,:) !nbasis,nbasis,nbasis,nbasis
    !========================================================
    logical :: Test 

    Test =.TRUE.
    ! Init stuff
    ! **********
    nbasis = MyMolecule%nbasis
    nocc   = MyMolecule%nocc
    nvirt  = MyMolecule%nunocc
    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)
    noccfull = nocc
    IF(DECinfo%frozencore)call lsquit('RI-MP2-F12 frozen core not implemented',-1)

    !==========================================================
    !=                                                        =
    !=             Step 1  Fijkl                              =
    !=                                                        =
    !==========================================================
    call II_get_CoulombEcont(DECinfo%output,DECinfo%output,mylsitem%setting,[
    CoulombF12 = Econt(1) 
    Econt(1) = 0.0E0_realk
    call II_get_exchangeEcont(DECinfo%output,DECinfo%output,mylsitem%setting,
    ExchangeF12 = Econt(1)       
    E21 = -0.5E0_realk*((5.0E0_realk/4.0E0_realk)*CoulombF12+ExchangeF12*0.5E
    WRITE(LUPRI,*)'E(Fijkl,LS) = ',E21

    IF(Test)THEN
       call mem_alloc(Vijij,nocc,nocc)
       call mem_alloc(Vjiij,nocc,nocc)
       call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
       call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRF')
       call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
            &                          MyMolecule%Co, MyMolecule%Cv,'iiii',gAO,Fijkl)
       call mem_dealloc(gao)
       do j=1,nocc
          do i=1,nocc
             Vijij(i,j) = Fijkl(i,i,j,j)
          enddo
       enddo
       do j=1,nocc
          do i=1,nocc
             Vjiij(i,j) = Fijkl(i,j,j,i)
          enddo
       enddo
       call mem_dealloc(Fijkl)
       E21 = 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)
       call mem_dealloc(Vijij)
       call mem_dealloc(Vjiij)
       WRITE(LUPRI,*)'E(Fijkl,FullIntegral) = ',E21
    ENDIF
    !==========================================================
    !=                                                        =
    !=             Step 2  Ripjq*Gipjq                        =
    !=                                                        =
    !==========================================================
    
!!$    IF(Test)THEN
!!$       do q=1,nbasis
!!$          do j=1,nocc
!!$             do p=1,nbasis
!!$                do i=1,nocc
!!$                   Vijij(i,j) = Vijij(i,j) - Ripjq(i,p,j,q)*Gipjq(i,p,j,q)
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo
!!$       
!!$       do q=1,nbasis
!!$          do j=1,nocc
!!$             do p=1,nbasis
!!$                do i=1,nocc
!!$                   Vjiij(i,j) = Vjiij(i,j) - Ripjq(i,p,j,q)*Gipjq(j,p,i,q)
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo
       
!       E21 = 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)
!       WRITE(LUPRI,*)'E(Ripjq*Gipjq,FullIntegral) = ',E21
    ENDIF
  end subroutine full_canonical_mp2_f12

  !> Function for finding the E21 energy  
  function mp2f12_E21(Vijij,Vjiij,nocc) result(energy)
    implicit none
    Integer,intent(IN)     :: nocc
    Real(realk),intent(IN) :: Vijij(nocc,nocc),Vjiij(nocc,nocc)
    Real(realk) :: energy
    !
    Integer     :: i,j
    Real(realk) :: tmp
    
    tmp = 0E0_realk
    DO i=1,nocc
       tmp = tmp + Vijij(i,i)
    ENDDO
    
    energy = -0.5E0_realk*tmp
    tmp = 0E0_realk
    
    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp + 5E0_realk * Vijij(i,j) - Vjiij(i,j)
       ENDDO
    ENDDO
    energy = energy - 0.25E0_realk*tmp
  end function mp2f12_E21
  
end module fullrimp2f12

