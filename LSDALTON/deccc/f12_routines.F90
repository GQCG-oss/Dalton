!> @file
!> Full calculation
!> This file is mainly a playground for new developments, not intended to be included in a release.

module f12_routines_module
  use fundamental
  use precision
  use typedeftype!, only:lsitem
  use typedef
  use dec_typedef_module
  use matrix_module
  use matrix_operations
  use memory_handling

  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use dec_fragment_utils
  use CABS_operations

#ifdef MOD_UNRELEASED
  use full_f12contractions
#endif 
  use ccintegrals!,only: get_full_AO_integrals,get_AO_hJ,get_AO_K,get_AO_Fock

  public :: MO_transform_AOMatrix, get_F12_mixed_MO_Matrices_real, get_F12_mixed_MO_Matrices, free_F12_mixed_MO_Matrices, &
       & free_F12_mixed_MO_Matrices_real

  private

contains

  !> Needs documentation
  subroutine get_F12_mixed_MO_Matrices_real(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
       & nocc,noccfull,nvirt,ncabs,HJir_real,Krr_real,Frr_real,Fac_real,Fpp_real,Fii_real,Fmm_real,Frm_real,Fcp_real)

    implicit none
    !> Fragmet molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: MyLsitem
    integer :: nbasis,nocc,nvirt,noccfull,ncabsAO,ncabs
    type(matrix), intent(in) :: Dmat
  
    real(realk), intent(inout) :: HJir_real(nocc,ncabsAO) 
    real(realk), intent(inout) :: Krr_real(ncabsAO,ncabsAO)
    real(realk), intent(inout) :: Frr_real(ncabsAO,ncabsAO)
    real(realk), intent(inout) :: Fac_real(nvirt,ncabs)
    real(realk), intent(inout) :: Fpp_real(nbasis,nbasis)
    real(realk), intent(inout) :: Fii_real(nocc,nocc)
    real(realk), intent(inout) :: Fmm_real(nocc,nocc)
    real(realk), intent(inout) :: Frm_real(ncabsAO,nocc)
    real(realk), intent(inout) :: Fcp_real(ncabs,nbasis)

    type(matrix) :: HJir
    type(matrix) :: Krr
    type(matrix) :: Frr
    type(matrix) :: Fac
    type(matrix) :: Fpp
    type(matrix) :: Fii
    type(matrix) :: Fmm
    type(matrix) :: Frm
    type(matrix) :: Fcp

    !> Temp
    type(matrix) :: HJrc
    type(matrix) :: Kcc
    type(matrix) :: Fcc
    type(matrix) :: Frc  
   
    !> Mixed regular/CABS one-electron and Coulomb matrix (h+J) combination in AO basis
    !> hJir
    call mat_init(HJrc,nbasis,ncabsAO)
    call get_AO_hJ(nbasis,ncabsAO,HJrc,Dmat,MyLsitem,'RCRRC')
    call mat_init(HJir,nocc,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'ir',HJrc,HJir)
    call mat_to_full(HJir,1.0E0_realk,HJir_real)
    call mat_free(HJrc)
    call mat_free(HJir)

    !> Mixed CABS/CABS exchange matrix
    !> Krr
    call mat_init(Kcc,ncabsAO,ncabsAO)
    call get_AO_K(nbasis,ncabsAO,Kcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Krr,ncabsAO,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'rr',Kcc,Krr)
    call mat_free(Kcc)
    call mat_to_full(Krr,1.0E0_realk,Krr_real)
    call mat_free(Krr)
    
    !> Mixed CABS/CABS Fock matrix
    !> Frr
    call mat_init(Fcc,ncabsAO,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Frr,ncabsAO,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'rr',Fcc,Frr)
    call mat_free(Fcc)   
    call mat_to_full(Frr,1.0E0_realk,Frr_real)
    call mat_free(Frr)
    
    !> Mixed AO/CABS Fock matrix
    !> Fac
    call mat_init(Frc,nbasis,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Frc,Dmat,MyLsitem,'RCRRC')
    call mat_init(Fac,nvirt,ncabs)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'ac',Frc,Fac)    
    call mat_free(Frc)
    call mat_to_full(Fac,1.0E0_realk,Fac_real)
    call mat_free(Fac)

    !> Mixed AO/AO full MO Fock matrix
    !> Temp Fcc
    call mat_init(Fcc,nbasis,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'RRRRC')

    !> Fpp
    call mat_init(Fpp,nbasis,nbasis)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'pp',Fcc,Fpp)
    call mat_to_full(Fpp,1.0E0_realk,Fpp_real)
    call mat_free(Fpp)
    
    !> Fii 
    call mat_init(Fii,nocc,nocc)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'ii',Fcc,Fii)
    call mat_to_full(Fii,1.0E0_realk,Fii_real)
    call mat_free(Fii)

    !> Fmm
    call mat_init(Fmm,noccfull,noccfull)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'mm',Fcc,Fmm)
    call mat_to_full(Fmm,1.0E0_realk,Fmm_real)
    call mat_free(Fmm)

    !> Free temp Fcc
    call mat_free(Fcc)

    !> Mixed CABS/AO MO Fock matrix
    !> Temp Fcc
    call mat_init(Fcc,ncabsAO,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CRRRC')
  
    !> Frm
    call mat_init(Frm,ncabsAO,noccfull)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'rm',Fcc,Frm)
    call mat_to_full(Frm,1.0E0_realk,Frm_real)
    call mat_free(Frm)
    
    !> Fcp
    call mat_init(Fcp,ncabs,nbasis)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'cp',Fcc,Fcp)
    call mat_to_full(Fcp,1.0E0_realk,Fcp_real)
    call mat_free(Fcp)

    !> Free Temp Fcc
    call mat_free(Fcc) 

!!$    print *, '****************************************'
!!$    print *, '(Norm of HJir_real):', norm2(HJir_real)       
!!$    print *, '(Norm of HJir):', sqrt(mat_sqnorm2(HJir))
!!$    print *, '****************************************'
        
  end subroutine get_F12_mixed_MO_Matrices_real

  !> Need documentation...
  subroutine get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
       & nocc,noccfull,nvirt,ncabs,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: Mylsitem
    integer :: nbasis,nocc,nvirt,noccfull,ncabsAO,ncabs
    type(matrix) :: Dmat,K
    type(matrix) :: HJir
    type(matrix) :: Krr
    type(matrix) :: Frr
    type(matrix) :: Frc
    type(matrix) :: Fpp
    type(matrix) :: Fmm
    type(matrix) :: Frm
    type(matrix) :: Fcp
    type(matrix) :: Fii
    type(matrix) :: Fac

    ! Temp
    type(matrix) :: HJrc
    type(matrix) :: Kcc
    type(matrix) :: Fcc

    ! Mixed regular/CABS one-electron and Coulomb matrix (h+J) combination in AO basis
    call mat_init(HJrc,nbasis,ncabsAO)
    call get_AO_hJ(nbasis,ncabsAO,HJrc,Dmat,MyLsitem,'RCRRC')
    call mat_init(HJir,nocc,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'ir',HJrc,HJir)
  
!!$    print *, '---------------------------------------------------'
!!$    print *, '---------------------------------------------------'
!!$    print *, '         Inside get_F12_mixed_MO_Matrices          '   
!!$    print *, '---------------------------------------------------'
!!$    print *, '---------------------------------------------------'
!!$    print *, 'nbabsis:', nbasis
!!$    print *, 'ncabsAO:', ncabsAO
!!$    print *, 'nocc:   ', nocc
!!$    print *, 'noccfull:', noccfull
!!$    print *, 'nvirt:  ', nvirt
!!$    print *, 'sqrt(mat_sqnorm2(Dmat)):', sqrt(mat_sqnorm2(Dmat))
!!$    print *, 'sqrt(mat_sqnorm2(HJrc)):', sqrt(mat_sqnorm2(HJrc))
!!$    print *, 'sqrt(mat_sqnorm2(HJir)):', sqrt(mat_sqnorm2(HJir))
    
    call mat_free(HJrc)

    ! Mixed CABS/CABS exchange matrix
    call mat_init(Kcc,ncabsAO,ncabsAO)
    call get_AO_K(nbasis,ncabsAO,Kcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Krr,ncabsAO,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'rr',Kcc,Krr)
    call mat_free(Kcc)

    ! Mixed CABS/CABS Fock matrix
    call mat_init(Fcc,ncabsAO,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Frr,ncabsAO,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'rr',Fcc,Frr)
    call mat_free(Fcc)

    ! Mixed AO/CABS Fock matrix
    call mat_init(Frc,nbasis,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Frc,Dmat,MyLsitem,'RCRRC')
    call mat_init(Fac,nvirt,ncabs)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'ac',Frc,Fac)
    call mat_free(Frc)

    ! Mixed AO/AO full MO Fock matrix
    call mat_init(Fcc,nbasis,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'RRRRC')
    !Fpp
    call mat_init(Fpp,nbasis,nbasis)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'pp',Fcc,Fpp)
    !Fii
    call mat_init(Fii,nocc,nocc)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'ii',Fcc,Fii)
    !Fmm
    call mat_init(Fmm,noccfull,noccfull)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'mm',Fcc,Fmm)
    call mat_free(Fcc)

    ! Mixed CABS/AO MO Fock matrix
    call mat_init(Fcc,ncabsAO,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CRRRC')
    !Frm
    call mat_init(Frm,ncabsAO,noccfull)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'rm',Fcc,Frm)
    !Fcc
    call mat_init(Fcp,ncabs,nbasis)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'cp',Fcc,Fcp)
    call mat_free(Fcc)
  end subroutine get_F12_mixed_MO_Matrices

  subroutine free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

    implicit none
    type(matrix) :: HJir
    type(matrix) :: Krr
    type(matrix) :: Frr
    type(matrix) :: Frc
    type(matrix) :: Fpp
    type(matrix) :: Fmm
    type(matrix) :: Frm
    type(matrix) :: Fcp
    type(matrix) :: Fii
    type(matrix) :: Fac

    call mat_free(HJir)
    call mat_free(Krr)
    call mat_free(Frr)
    call mat_free(Fac)
    call mat_free(Fpp)
    call mat_free(Fii)
    call mat_free(Fmm)
    call mat_free(Frm)
    call mat_free(Fcp)

  end subroutine free_F12_mixed_MO_Matrices

  subroutine free_F12_mixed_MO_Matrices_real(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

    implicit none  
    real(realk), pointer :: HJir(:,:) 
    real(realk), pointer :: Krr(:,:)
    real(realk), pointer :: Frr(:,:)
    real(realk), pointer :: Fac(:,:)
    real(realk), pointer :: Fpp(:,:)
    real(realk), pointer :: Fii(:,:)
    real(realk), pointer :: Fmm(:,:)
    real(realk), pointer :: Frm(:,:)
    real(realk), pointer :: Fcp(:,:)

    call mem_dealloc(HJir)
    call mem_dealloc(Krr)
    call mem_dealloc(Frr)
    call mem_dealloc(Fac)
    call mem_dealloc(Fpp)
    call mem_dealloc(Fii)
    call mem_dealloc(Fmm)
    call mem_dealloc(Frm)
    call mem_dealloc(Fcp)

  end subroutine free_F12_mixed_MO_Matrices_real

  !> Need documentation...
  subroutine MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & Cocc,Cvirt,inputstring,matAO,matMO)
    implicit none
    !> Lsitem structure
    integer :: nocc,noccfull,nvirt,nCabsAO,nCabs,nbasis
    type(lsitem), intent(inout) :: mylsitem
    integer :: ndim2(2),ndim1(2)
    type(matrix) :: matAO,matMO
    real(realk),pointer :: elms(:)
    type(matrix) :: CMO(2)
    real(realk),dimension(nbasis,nocc),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk),dimension(nbasis,nvirt),intent(in) :: Cvirt
    type(matrix) :: CMO_cabs,CMO_ri,tmp
    character(len=2) :: inputstring
    logical :: doCABS,doRI
    integer :: i,lupri
    character :: string(2)
    string(1)=inputstring(1:1) 
    string(2)=inputstring(2:2) 
    lupri=6
    doCABS = .FALSE.
    do i=1,2
       if(string(i).EQ.'c')then !cabs
          doCABS = .TRUE.
       endif
    enddo
    doRI = .FALSE.
    do i=1,2
       if(string(i).EQ.'r')then !RI
          doRI = .TRUE.
       endif
    enddo
    call determine_CABS_nbast(nCabsAO,nCabs,mylsitem%SETTING,lupri)
    IF(doCABS)THEN
       call mat_init(CMO_cabs,nCabsAO,nCabs)
       call build_CABS_MO(CMO_cabs,nCabsAO,mylsitem%SETTING,lupri)
    ENDIF
    IF(doRI)THEN
       call mat_init(CMO_ri,nCabsAO,nCabsAO)
       call build_RI_MO(CMO_ri,nCabsAO,mylsitem%SETTING,lupri)
    ENDIF
    do i=1,2
       if(string(i).EQ.'i')then !occupied active
          ndim1(i) = nbasis
          ndim2(i) = nocc
       elseif(string(i).EQ.'m')then !all occupied
          ndim1(i) = nbasis
          ndim2(i) = noccfull
       elseif(string(i).EQ.'p')then !all occupied + virtual
          ndim1(i) = nbasis
          ndim2(i) = nbasis
       elseif(string(i).EQ.'a')then !virtual
          ndim1(i) = nbasis
          ndim2(i) = nvirt
       elseif(string(i).EQ.'c')then !cabs
          ndim1(i) = ncabsAO
          ndim2(i) = ncabs
       elseif(string(i).EQ.'r')then !ri - MOs
          ndim1(i) = ncabsAO
          ndim2(i) = ncabsAO
       endif
       call mat_init(CMO(i),ndim1(i),ndim2(i))
       if(string(i).EQ.'i')then !occupied active
          call dcopy(ndim2(i)*ndim1(i),Cocc,1,CMO(i)%elms,1)
       elseif(string(i).EQ.'m')then !all occupied
          call dcopy(ndim2(i)*ndim1(i),Cocc,1,CMO(i)%elms,1)
       elseif(string(i).EQ.'p')then !all occupied + virtual
          call dcopy(noccfull*nbasis,Cocc,1,CMO(i)%elms,1)
          call dcopy(nvirt*nbasis,Cvirt,1,CMO(i)%elms(noccfull*nbasis+1:nbasis*nbasis),1)
       elseif(string(i).EQ.'a')then !virtual
          call dcopy(ndim2(i)*ndim1(i),Cvirt,1,CMO(i)%elms,1)
       elseif(string(i).EQ.'c')then !cabs
          call dcopy(ndim2(i)*ndim1(i),CMO_cabs%elms,1,CMO(i)%elms,1)
       elseif(string(i).EQ.'r')then !ri - MOs
          call dcopy(ndim2(i)*ndim1(i),CMO_RI%elms,1,CMO(i)%elms,1)
       endif
    enddo
    IF(doCABS)THEN
       call mat_free(CMO_cabs)
    ENDIF
    IF(doRI)THEN
       call mat_free(CMO_ri)
    ENDIF
    call mat_init(tmp,CMO(1)%ncol,matAO%ncol)
    call mat_mul(CMO(1),matAO,'t','n',1E0_realk,0E0_realk,tmp)
    call mat_mul(tmp,CMO(2),'n','n',1E0_realk,0E0_realk,matMO)
    call mat_free(tmp)
    do i=1,2
       call mat_free(CMO(i))
    enddo

  end subroutine MO_transform_AOMatrix



end module f12_routines_module
