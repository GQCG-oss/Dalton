!> @file
!> Full calculation
!> This file is mainly a playground for new developments, not intended to be included in a release.

module f12_routines_module
  use fundamental
  use precision
  use typedeftype!,only:lsitem
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
  
  public :: MO_transform_AOMatrix
  private

contains

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
          call dcopy(nvirt*nbasis,Cvirt,1,CMO(i)%elms(noccfull*nbasis+1),1)
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
