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
  use CABS_operations

#ifdef MOD_UNRELEASED
  use full_f12contractions
#endif 
  use ccintegrals  
 

  public :: MO_transform_AOMatrix, get_F12_mixed_MO_Matrices_real, get_F12_mixed_MO_Matrices, free_F12_mixed_MO_Matrices, &
       & free_F12_mixed_MO_Matrices_real, norm1D, norm2D, norm4D,&
       & F12_RI_transform_realMat, F12_CABS_transform_realMat ! atomic_fragment_free_f12, atomic_fragment_init_f12

  private
  
contains
  function norm1D(A)
    implicit none
    real(realk), intent(in) :: A(:)
    integer :: m,i
    real(realk) :: norm1D   
    
    norm1D = 0.0E0_realk
    
    m = size(A,1)
    do i=1,m
       norm1D = norm1D + A(i)
    enddo
  end function norm1D


  !> Takes the norm of a matrix 
  function norm2D(A)
    implicit none
    real(realk), intent(in) :: A(:,:)
    integer :: m,n,i,j
    real(realk) :: norm2D   
    
    norm2D = 0.0E0_realk
    
    m = size(A,1)
    n = size(A,2)
    
    do j=1,n
       do i=1,m
          norm2D = norm2D + A(i,j)
       enddo
    enddo
    
  end function norm2D

  function norm4D(A)
    implicit none
    real(realk), intent(in) :: A(:,:,:,:)
    integer :: m,n,p,q,i,j,k,l
    real(realk) :: norm4D   
    
    norm4D = 0.0E0_realk
    
    m = size(A,1)
    n = size(A,2)
    p = size(A,3)
    q = size(A,4)
  
    do l=1,q
       do k=1,p
          do j=1,n
             do i=1,m 
                norm4D = norm4D + A(i,j,k,l)*A(i,j,k,l)
             enddo
          enddo
       enddo
    enddo

   norm4D = sqrt(norm4D)

  end function norm4D


  !> Needs documentation
!WARNING THIS IS NOT THE SAME AS get_F12_mixed_MO_Matrices
!THIS SUBROUTINE ONLY DOES HALF TRANSFORMATIONS 
!ALL TRANSFORMATIONS RELATED TO CABS AND RI IS NOT DONE YET
!THIS IS DONE AT THE FRAGMENT LEVEL. 
  subroutine get_F12_mixed_MO_Matrices_real(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
       & nocc,noccfull,nvirt,HJir_real,Krr_real,Frr_real,Fac_real,Fii_real,Frm_real,Fcp_real)

    implicit none
    !> Fragmet molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: MyLsitem
    integer :: nbasis,nocc,nvirt,noccfull,ncabsAO!,ncabs
    type(matrix), intent(in) :: Dmat
  
    real(realk), intent(inout) :: HJir_real(nocc,ncabsAO)  !ONLY HALF TRANSFORMED
    real(realk), intent(inout) :: Krr_real(ncabsAO,ncabsAO)!NOT TRANSFORMED
    real(realk), intent(inout) :: Frr_real(ncabsAO,ncabsAO)!NOT TRANSFORMED
    real(realk), intent(inout) :: Fac_real(nvirt,ncabsAO)  !HACK not (nvirt,ncabsMO)
    real(realk), intent(inout) :: Fii_real(nocc,nocc)
    real(realk), intent(inout) :: Frm_real(ncabsAO,nocc)   !ONLY HALF TRANSFORMED
    real(realk), intent(inout) :: Fcp_real(ncabsAO,nbasis) !HACK not (ncabsMO,nbasis)

    type(matrix) :: HJir
    type(matrix) :: Krr
    type(matrix) :: Frr
    type(matrix) :: Fac
    type(matrix) :: Fii
    type(matrix) :: Frm
    type(matrix) :: Fcp

    !> Temp
    type(matrix) :: HJrc
    type(matrix) :: Kcc
    type(matrix) :: Fcc
    type(matrix) :: Frc  
    type(matrix) :: Fcr  
   
    !> Mixed regular/CABS one-electron and Coulomb matrix (h+J) combination in AO basis
    !> hJir
    call mat_init(HJrc,nbasis,ncabsAO)
    call get_AO_hJ(nbasis,ncabsAO,HJrc,Dmat,MyLsitem,'RCRRC')
    call mat_init(HJir,nocc,ncabsAO)
    call MO_halftransform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'i',HJrc,HJir,1)
    call mat_free(HJrc)
    call mat_to_full(HJir,1.0E0_realk,HJir_real)
    call mat_free(HJir)

    !> Mixed CABS/CABS exchange matrix 
    !> Krr
    call mat_init(Kcc,ncabsAO,ncabsAO)
    call get_AO_K(nbasis,ncabsAO,Kcc,Dmat,MyLsitem,'CCRRC')
    call mat_to_full(Kcc,1.0E0_realk,Krr_real)
    call mat_free(Kcc)
    
    !> Mixed CABS/CABS Fock matrix
    !> Frr
    call mat_init(Fcc,ncabsAO,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CCRRC')
    call mat_to_full(Fcc,1.0E0_realk,Frr_real)
    call mat_free(Fcc)
    
    !> Mixed AO/CABS Fock matrix
    !> Fac
    call mat_init(Frc,nbasis,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Frc,Dmat,MyLsitem,'RCRRC')
    call mat_init(Fac,nvirt,ncabsAO)
    call MO_halftransform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'a',Frc,Fac,1)    
    call mat_free(Frc)
    call mat_to_full(Fac,1.0E0_realk,Fac_real)
    call mat_free(Fac)

    !> Mixed AO/AO full MO Fock matrix
    !> Fii 
    call mat_init(Frr,nbasis,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Frr,Dmat,MyLsitem,'RRRRC')
    call mat_init(Fii,nocc,nocc)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'ii',Frr,Fii)
    call mat_free(Frr)
    call mat_to_full(Fii,1.0E0_realk,Fii_real)
    call mat_free(Fii)


    !> Mixed CABS/AO MO Fock matrix
    call mat_init(Fcr,ncabsAO,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcr,Dmat,MyLsitem,'CRRRC')

    call mat_init(Frm,ncabsAO,noccfull)
    call MO_halftransform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'m',Fcr,Frm,2)
    call mat_to_full(Frm,1.0E0_realk,Frm_real)
    call mat_free(Frm)
    
    call mat_init(Fcp,ncabsAO,nbasis)
    call MO_halftransform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'p',Fcr,Fcp,2)
    call mat_to_full(Fcp,1.0E0_realk,Fcp_real)
    call mat_free(Fcp)

    call mat_free(Fcr) 

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
    type(matrix) :: Frm
    type(matrix) :: Fcp
    type(matrix) :: Fpp
    type(matrix) :: Fmm
    type(matrix) :: Fii
    type(matrix) :: Fac

    call mat_free(HJir)
    call mat_free(Krr)
    call mat_free(Frr)
    call mat_free(Fac)
    call mat_free(Fii)
    call mat_free(Fpp)
    call mat_free(Fmm)
    call mat_free(Frm)
    call mat_free(Fcp)

  end subroutine free_F12_mixed_MO_Matrices

  subroutine free_F12_mixed_MO_Matrices_real(HJir,Krr,Frr,Fac,Fii,Frm,Fcp)

    implicit none  
    real(realk), pointer :: HJir(:,:) 
    real(realk), pointer :: Krr(:,:)
    real(realk), pointer :: Frr(:,:)
    real(realk), pointer :: Fac(:,:)
    real(realk), pointer :: Fii(:,:)
    real(realk), pointer :: Frm(:,:)
    real(realk), pointer :: Fcp(:,:)

    call mem_dealloc(HJir)
    call mem_dealloc(Krr)
    call mem_dealloc(Frr)
    call mem_dealloc(Fac)
    call mem_dealloc(Fii)
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

  !> Need documentation...
  !> Half transform a AO matrix to MO matrix
  subroutine MO_halftransform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & Cocc,Cvirt,inputstring,matAO,matMO,option)
    implicit none
    !> Lsitem structure
    integer :: nocc,noccfull,nvirt,nCabsAO,nbasis,option
    type(lsitem), intent(inout) :: mylsitem
    integer :: ndim2(1),ndim1(1)
    type(matrix) :: matAO,matMO
    real(realk),pointer :: elms(:)
    type(matrix) :: CMO(1)
    real(realk),dimension(nbasis,nocc),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk),dimension(nbasis,nvirt),intent(in) :: Cvirt
    type(matrix) :: CMO_cabs,CMO_ri,tmp
    character(len=1) :: inputstring
    logical :: doCABS,doRI
    integer :: i,lupri
    character :: string(1)
    string(1)=inputstring(1:1) 
    lupri=6
    doCABS = .FALSE.
       if(string(1).EQ.'i')then !occupied active
          ndim1(1) = nbasis
          ndim2(1) = nocc
       elseif(string(1).EQ.'m')then !all occupied
          ndim1(1) = nbasis
          ndim2(1) = noccfull
       elseif(string(1).EQ.'p')then !all occupied + virtual
          ndim1(1) = nbasis
          ndim2(1) = nbasis
       elseif(string(1).EQ.'a')then !virtual
          ndim1(1) = nbasis
          ndim2(1) = nvirt
       elseif(string(1).EQ.'c')then !cabs
          call lsquit('input error MO_transform1_AOMatrix',-1)
       elseif(string(1).EQ.'r')then !ri - MOs
          call lsquit('input error MO_transform1_AOMatrix',-1)
       endif
       call mat_init(CMO(1),ndim1(1),ndim2(1))
       if(string(1).EQ.'i')then !occupied active
          call dcopy(ndim2(1)*ndim1(1),Cocc,1,CMO(1)%elms,1)
       elseif(string(1).EQ.'m')then !all occupied
          call dcopy(ndim2(1)*ndim1(1),Cocc,1,CMO(1)%elms,1)
       elseif(string(1).EQ.'p')then !all occupied + virtual
          call dcopy(noccfull*nbasis,Cocc,1,CMO(1)%elms,1)
          call dcopy(nvirt*nbasis,Cvirt,1,CMO(1)%elms(noccfull*nbasis+1:nbasis*nbasis),1)
       elseif(string(1).EQ.'a')then !virtual
          call dcopy(ndim2(1)*ndim1(1),Cvirt,1,CMO(1)%elms,1)
       endif
       IF(option.EQ.1)THEN
          call mat_mul(CMO(1),matAO,'t','n',1E0_realk,0E0_realk,matMO)
       ELSE
          call mat_mul(matAO,CMO(1),'n','n',1E0_realk,0E0_realk,matMO)
       ENDIF
    call mat_free(CMO(1))
  end subroutine MO_halftransform_AOMatrix

  subroutine F12_RI_transform_realMat(Mat,ndim1,ndim2,Cri,nCabsAO)
    implicit none
    integer,intent(in) :: ndim1,ndim2,nCabsAO
    real(realk),intent(inout) :: Mat(ndim1,ndim2)
    real(realk),intent(in) :: Cri(nCabsAO,nCabsAO)
    !local variables
    real(realk),pointer :: tmp(:,:)
    call mem_alloc(tmp,ndim1,ndim2)
    IF((ndim1.EQ.nCabsAO).AND.(ndim2.EQ.nCabsAO))THEN
       !       call DGEMM('t','n',Cri%ncol,Mat%ncol,Cri%nrow,1E0_realk,&
       !            &Cri,Cri%nrow,Mat,Mat%nrow,0E0_realk,tmp,tmp%nrow)
       !       call DGEMM('n','n',Tmp%nrow,Cri%ncol,Tmp%ncol,1E0_realk,&
       !            &Tmp,Tmp%nrow,Cri,Cri%nrow,0E0_realk,Mat,Mat%nrow)
       call DGEMM('t','n',nCabsAO,ndim2,nCabsAO,1E0_realk,&
            &Cri,nCabsAO,Mat,ndim1,0E0_realk,tmp,ndim1)
       call DGEMM('n','n',ndim1,nCabsAO,ndim2,1E0_realk,&
            &Tmp,ndim1,Cri,nCabsAO,0E0_realk,Mat,ndim1)
    ELSEIF(ndim1.EQ.nCabsAO)THEN
       call DGEMM('t','n',nCabsAO,ndim2,nCabsAO,1E0_realk,&
            &Cri,nCabsAO,Mat,ndim1,0E0_realk,tmp,ndim1)
       call DCOPY(ndim1*ndim2,Tmp,1,Mat,1)
    ELSEIF(ndim2.EQ.nCabsAO)THEN
       call DGEMM('n','n',ndim1,nCabsAO,ndim2,1E0_realk,&
            &Mat,ndim1,Cri,nCabsAO,0E0_realk,Tmp,ndim1)
       call DCOPY(ndim1*ndim2,Tmp,1,Mat,1)
    ELSE
       call lsquit('unknown option in F12_RI_transform_realMat',-1)
    ENDIF
    call mem_dealloc(tmp)

  end subroutine F12_RI_transform_realMat

  subroutine F12_CABS_transform_realMat(MOmat,AOmat,ndim1,ndim2,Ccabs,nCabsAO,nCabsMO)
    implicit none
    integer,intent(in) :: ndim1,ndim2,nCabsAO,nCabsMO
    real(realk),intent(in) :: AOMat(ndim1,ndim2)
    real(realk),intent(inout) :: MOMat(nCabsMO,ndim2)
    real(realk),intent(in) :: Ccabs(nCabsAO,nCabsMO)

    IF(ndim1.EQ.nCabsAO)THEN
       !       call DGEMM('t','n',Cri%ncol,Mat%ncol,Cri%nrow,1E0_realk,&
       !            &Cri,Cri%nrow,Mat,Mat%nrow,0E0_realk,tmp,tmp%nrow)
       call DGEMM('t','n',nCabsMO,ndim2,nCabsAO,1E0_realk,&
            &Ccabs,nCabsAO,AOMat,ndim1,0E0_realk,MOMat,nCabsMO)
    ELSE
       call lsquit('unknown option in F12_RI_transform_realMat',-1)
    ENDIF

  end subroutine F12_CABS_transform_realMat

!!$  !> \brief Initialize atomic fragment based on a list of specific AOS orbitals for F12-related matrices
!!$  !> \author Yang Min Wang
!!$  subroutine atomic_fragment_init_f12(fragment, MyMolecule)
!!$    type(fullmolecule), intent(in) :: MyMolecule
!!$    type(decfrag), intent(inout) :: fragment
!!$    
!!$    !> F12 Specific Variables
!!$    integer :: nbasis, noccEOS, nunoccEOS, noccfull, nocvAOS, nvirtAOS, ncabsAO, ncabsMO
!!$    integer :: noccAOS, nunoccAOS
!!$    integer :: ix, iy
!!$
!!$    ncabsAO = size(MyMolecule%Ccabs,1)
!!$    ncabsMO = size(MyMolecule%Ccabs,2)
!!$
!!$    call mem_alloc(fragment%Ccabs,ncabsAO,ncabsMO)
!!$    call dcopy(ncabsAO*ncabsMO,Mymolecule%Ccabs,1,fragment%Ccabs,1)
!!$    
!!$    call mem_alloc(fragment%Cri,ncabsAO,ncabsAO)
!!$    call dcopy(ncabsAO*ncabsAO,Mymolecule%Cri,1,fragment%Cri,1)
!!$    
!!$    ! ============================================================
!!$    !                        F12-Specific                        !
!!$    ! ============================================================
!!$    !> F12 Specific Variables
!!$    nbasis   = fragment%nbasis
!!$    noccEOS  = fragment%noccEOS
!!$    nunoccEOS = fragment%nunoccEOS
!!$ 
!!$    noccAOS  = fragment%noccAOS
!!$    nunoccAOS = fragment%nunoccAOS  
!!$    nocvAOS  = fragment%noccAOS + fragment%nunoccAOS
!!$    nvirtAOS = fragment%nunoccAOS
!!$    ncabsAO = size(fragment%Ccabs,1)    
!!$    ncabsMO = size(fragment%Ccabs,2)    
!!$    
!!$      
!!$       ! ************************************************************
!!$       ! Creating a CocvAOS matrix 
!!$       ! ************************************1***********************
!!$       !do i=1, Fragment%noccAOS
!!$       !   CocvAOS(:,i) = Fragment%Co(:,i)
!!$       !end do
!!$       
!!$       print *,"nbasis:   ", nbasis
!!$       print *,"noccEOS:  ", noccEOS
!!$       print *,"nunoccEOS:", nunoccEOS
!!$       
!!$       print *,"nocvAOS:  ", nocvAOS
!!$       print *,"noccAOS: ",  noccAOS
!!$       print *,"nunoccAOS: ", nvirtAOS
!!$       print *,"ncabsAO:  ", ncabsAO
!!$       print *,"ncabsMO:  ", ncabsMO
!!$
!!$       ! Be carefull about defining what is EOS and AOS space
!!$       ! At the moment we have a EOS partitioning scheme not a CABS, they will 
!!$       ! eventually be dependent on the EOS i, j, because we will project out 
!!$       ! from the i and j and create the CABS 
!!$       ! We have partitioned the EOS space
!!$       
!!$       ! hJir
!!$       call mem_alloc(fragment%hJir, noccEOS, ncabsAO)
!!$       do j=1,ncabsAO
!!$          do i=1, fragment%noccEOS
!!$             ix = fragment%idxo(i)
!!$             fragment%hJir(i,j) = MyMolecule%hJir(ix,j)
!!$          enddo
!!$       enddo
!!$       print *,"norm2(MyMolecule%hJir): ",norm2(MyMolecule%hJir)
!!$       print *,"norm2(fragment%hJir):   ",norm2(fragment%hJir)
!!$
!!$  end subroutine atomic_fragment_init_f12
!!$  
!!$  !> \brief Free the f12 fragment free matrices
!!$  !> \author Yang Min Wang
!!$  subroutine atomic_fragment_free_f12(fragment)
!!$    
!!$    implicit none
!!$    !> Atomic fragment to be freed
!!$    type(decfrag),intent(inout) :: fragment
!!$    
!!$    print *, "------------atomic_fragment_free_f12(fragment)----------"
!!$    
!!$    if(associated(fragment%hJir)) then
!!$       call mem_dealloc(fragment%hJir)
!!$       fragment%hJir => null()
!!$    end if
!!$    
!!$    if(associated(fragment%Krs)) then
!!$       call mem_dealloc(fragment%Krs)
!!$       fragment%Krs => null()
!!$    end if
!!$    
!!$    if(associated(fragment%Frs)) then
!!$       call mem_dealloc(fragment%Frs)
!!$       fragment%Frs => null()
!!$    end if
!!$
!!$    if(associated(fragment%Fac)) then
!!$       call mem_dealloc(fragment%Fac)
!!$       fragment%Fac => null()
!!$    end if
!!$    
!!$    if(associated(fragment%Fpq)) then
!!$       call mem_dealloc(fragment%Fpq)
!!$       fragment%Fpq => null()
!!$    end if
!!$    
!!$    if(associated(fragment%Fij)) then
!!$       call mem_dealloc(fragment%Fij)
!!$       fragment%Fij => null()
!!$    end if
!!$    
!!$    if(associated(fragment%Fmn)) then
!!$       call mem_dealloc(fragment%Fmn)
!!$       fragment%Fmn => null()
!!$    end if
!!$    
!!$  end subroutine atomic_fragment_free_f12






end module f12_routines_module
