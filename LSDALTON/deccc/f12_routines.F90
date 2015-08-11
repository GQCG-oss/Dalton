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
  use IchorErimoduleHost

  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use CABS_operations

  use full_f12contractions
  use ccintegrals  
 
  public :: MO_transform_AOMatrix, get_F12_mixed_MO_Matrices_real, get_F12_mixed_MO_Matrices, free_F12_mixed_MO_Matrices, &
       & free_F12_mixed_MO_Matrices_real, norm1D, norm2D, norm4D, &
       & F12_RI_transform_realMat, F12_CABS_transform_realMat, get_mp2f12_MO, & ! atomic_fragment_free_f12, atomic_fragment_init_f12
       & get_4Center_MO_integrals, get_4Center_F12_integrals, free_4Center_F12_integrals, &
       & mp2f12_Xijij_term3, mp2f12_Xjiij_term3, mp2f12_Xijij_term4, mp2f12_Xjiij_term4, &
       & get_ES2,get_ES2_from_dec_main,dec_get_RI_orbitals,dec_get_CABS_orbitals, get_mp2f12_MO_PDM, &
       & mp2f12_Bijij_term2, mp2f12_Bijij_term3, mp2f12_Bijij_term4, mp2f12_Bijij_term5, mp2f12_Bijij_term6, mp2f12_Bijij_term7, &
       & mp2f12_Bijij_term8, mp2f12_Bijij_term9, Contractocccalpha, &
       & ContractOne4CenterF12IntegralsRI, ContractOne4CenterF12IntegralsRI2, Contractone4centerf12integralsrib23, & 
       & Contractone4centerf12integralsrobustri, Contracttwo4centerf12integralsri2v3v4, Contracttwo4centerf12integralsriX3X4, &
       & Contracttwo4centerf12integralsri, Contracttwo4centerf12integralsrib4,Contracttwo4centerf12integralsrib5, &
       & Contracttwo4centerf12integralsrib6,Contracttwo4centerf12integralsrib7, Contracttwo4centerf12integralsrib8, &
       & Contracttwo4centerf12integralsrib9, Contracttwo4centerf12integralsric, Contracttwo4centerf12integralsriX, &
       & ContractOne4CenterF12IntegralsRI2_nc, ContractTwo4CenterF12IntegralsRIC_p, ContractOne4CenterF12IntegralsRI2_p, &
       & Contracttwo4centerf12integralsri2v3v4_p, Contracttwo4centerf12integralsriX3X4_nc, ContractTwo4CenterF12IntegralsRIX_nc, &
       & Contracttwo4centerf12integralsriX3X4_nc2
  private

  !> Coefficient Type
  TYPE ctype
     real(realk), pointer :: cmat(:,:)
     integer :: n1
     integer :: n2
  END TYPE ctype
  
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
   ! real(realk), intent(inout) :: Fcd_real(ncabsAO,ncabsAO)      
    
    !> Fock matrices for singles correction
    !real(realk), intent(inout), optional :: Fic_real(nocc,ncabsAO) !HACK not (ncabsMO,nbasis)
    !real(realk), intent(inout), optional :: Fcc_real(ncabsAO,ncabsAO) !HACK not (ncabsMO,nbasis)
        
    type(matrix) :: HJir
    type(matrix) :: Krr
    type(matrix) :: Frr
    type(matrix) :: Fac
    type(matrix) :: Fii
    type(matrix) :: Frm
    type(matrix) :: Fcp
 !   type(matrix) :: Fcd
    
    !> Temp
    type(matrix) :: HJrc
    type(matrix) :: Kcc
    type(matrix) :: Fcc
    type(matrix) :: Frc  
    type(matrix) :: Fcr  

    if( MyMolecule%mem_distributed )then
       call lsquit("ERROR(get_F12_mixed_MO_Matrices_real): this routine does not work&
       & with distributed arrays in the fullmolecule type, yet",-1)
    endif
   
    !> Mixed regular/CABS one-electron and Coulomb matrix (h+J) combination in AO basis
    !> hJir
    call mat_init(HJrc,nbasis,ncabsAO)
    call get_AO_hJ(nbasis,ncabsAO,HJrc,Dmat,MyLsitem,'RCRRC')
    call mat_init(HJir,nocc,ncabsAO)
    call MO_halftransform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'i',HJrc,HJir,1)
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
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'a',Frc,Fac,1)    
    call mat_free(Frc)
    call mat_to_full(Fac,1.0E0_realk,Fac_real)
    call mat_free(Fac)

    !> Mixed AO/AO full MO Fock matrix
    !> Fii 
    call mat_init(Frr,nbasis,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Frr,Dmat,MyLsitem,'RRRRC')
    call mat_init(Fii,nocc,nocc)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ii',Frr,Fii)
    call mat_free(Frr)
    call mat_to_full(Fii,1.0E0_realk,Fii_real)
    call mat_free(Fii)

    !> Mixed CABS/AO MO Fock matrix
    call mat_init(Fcr,ncabsAO,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcr,Dmat,MyLsitem,'CRRRC')

    call mat_init(Frm,ncabsAO,noccfull)
    call MO_halftransform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'m',Fcr,Frm,2)
    call mat_to_full(Frm,1.0E0_realk,Frm_real)
    call mat_free(Frm)
    
    call mat_init(Fcp,ncabsAO,nbasis)
    call MO_halftransform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'p',Fcr,Fcp,2)
    call mat_to_full(Fcp,1.0E0_realk,Fcp_real)
    call mat_free(Fcp)
    call mat_free(Fcr) 

!!$    !Fcd
!!$    call mat_init(Fcc,ncabsAO,ncabsAO)
!!$    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CCRRC')
!!$    call mat_init(Fcd,ncabsAO,ncabsAO)
!!$    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
!!$         & MyMolecule%Co, MyMolecule%Cv,'cc',Fcc,Fcd)
!!$    call mat_to_full(Fcd,1.0E0_realk,Fcd_real)
!!$    call mat_free(Fcd)
!!$    call mat_free(Fcc)
      
  end subroutine get_F12_mixed_MO_Matrices_real

  !> Need documentation...
  subroutine get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
       & nocc,noccfull,nvirt,ncabs,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)

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

    !> Singles contribution
    type(matrix) :: Fic
    type(matrix) :: Fcd
    
    ! Temp
    type(matrix) :: HJrc
    type(matrix) :: Kcc
    type(matrix) :: Fcc

    if( MyMolecule%mem_distributed )then
       call lsquit("ERROR(get_F12_mixed_MO_Matrices): this routine does not work&
       & with distributed arrays in the fullmolecule type, yet",-1)
    endif

    ! Mixed regular/CABS one-electron and Coulomb matrix (h+J) combination in AO basis
    !hJir
    call mat_init(HJrc,nbasis,ncabsAO)
    call get_AO_hJ(nbasis,ncabsAO,HJrc,Dmat,MyLsitem,'RCRRC')
    call mat_init(HJir,nocc,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ir',HJrc,HJir)    
    call mat_free(HJrc)

    ! Mixed CABS/CABS exchange matrix
    !Krr
    call mat_init(Kcc,ncabsAO,ncabsAO)
    call get_AO_K(nbasis,ncabsAO,Kcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Krr,ncabsAO,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'rr',Kcc,Krr)
    call mat_free(Kcc)

    ! Mixed CABS/CABS Fock matrix
    !Frr 
    call mat_init(Fcc,ncabsAO,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Frr,ncabsAO,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'rr',Fcc,Frr)
    call mat_free(Fcc)

    !Fcd
    call mat_init(Fcc,ncabsAO,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Fcd,ncabs,ncabs)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'cc',Fcc,Fcd)
    call mat_free(Fcc)
  
    ! Mixed AO/CABS Fock matrix
    !Fac
    call mat_init(Frc,nbasis,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Frc,Dmat,MyLsitem,'RCRRC')
    call mat_init(Fac,nvirt,ncabs)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ac',Frc,Fac)
    call mat_free(Frc)
       
    !Fic
    call mat_init(Frc,nbasis,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Frc,Dmat,MyLsitem,'RCRRC')
    call mat_init(Fic,nocc,ncabs)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ic',Frc,Fic)
    call mat_free(Frc)
    
    ! Mixed AO/AO full MO Fock matrix 
    call mat_init(Fcc,nbasis,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'RRRRC')
    !Fpp
    call mat_init(Fpp,nbasis,nbasis)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'pp',Fcc,Fpp)
    !Fii
    call mat_init(Fii,nocc,nocc)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ii',Fcc,Fii)
    !Fmm
    call mat_init(Fmm,noccfull,noccfull)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'mm',Fcc,Fmm)
    call mat_free(Fcc)

    ! Mixed CABS/AO MO Fock matrix
    call mat_init(Fcc,ncabsAO,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CRRRC')
    !Frm
    call mat_init(Frm,ncabsAO,noccfull)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'rm',Fcc,Frm)
    !Fcc
    call mat_init(Fcp,ncabs,nbasis)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'cp',Fcc,Fcp)
    call mat_free(Fcc)
    
  end subroutine get_F12_mixed_MO_Matrices

  subroutine free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)

    implicit none
    type(matrix),intent(inout) :: HJir
    type(matrix),intent(inout) :: Krr
    type(matrix),intent(inout) :: Frr
    type(matrix),intent(inout) :: Fac
    type(matrix),intent(inout) :: Frm
    type(matrix),intent(inout) :: Fcp
    type(matrix),intent(inout) :: Fpp
    type(matrix),intent(inout) :: Fmm
    type(matrix),intent(inout) :: Fii
    type(matrix),intent(inout) :: Fic
    type(matrix),intent(inout) :: Fcd

    call mat_free(HJir)
    call mat_free(Krr)
    call mat_free(Frr)
    call mat_free(Fac)
    call mat_free(Fii)
    call mat_free(Fpp)
    call mat_free(Fmm)
    call mat_free(Frm)
    call mat_free(Fcp)
    call mat_free(Fic)
    call mat_free(Fcd)
        
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
    real(realk),dimension(nbasis,noccfull),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk),dimension(nbasis,nvirt),intent(in) :: Cvirt
    type(matrix) :: CMO_cabs,CMO_ri,tmp
    character(len=2) :: inputstring
    logical :: doCABS,doRI
    integer :: i,lupri,offset
    character :: string(2)

    ! Offset:   Frozen core    : ncore
    !           Not frozen core: 0
    offset = noccfull - nocc

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
          call dcopy(ndim2(i)*ndim1(i),Cocc(1:nbasis,offset+1:noccfull),1,CMO(i)%elms,1)
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


  !> Brief: Get <1,2|INTSPEC|3,4> MO integrals wrapper.
  !> Author: Yang M. Wang
  !> Data: Nov 2013
  subroutine get_mp2f12_MO(MyFragment,MySetting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,INTTYPE,INTSPEC,transformed_mo)
    implicit none

    !> Atomic fragment to be determined  (NOT pair fragment)
    type(decfrag), intent(inout) :: MyFragment
    !> Integrals settings   
    type(lssetting), intent(inout) :: Mysetting
    !> Number of basis functions AO
    integer :: nbasis
    !> Number of occupied orbitals MO in EOS space
    integer :: noccEOS
    !> Number of occupied orbitals MO in AOS space (includes core, also for frozen core)
    integer :: noccAOStot
    !> Number of virtupied (virtual) orbitals MO in EOS space
    integer :: nvirtEOS
    !> Number of occupied + virtual MO in AOS space (includes core, also for frozen core)
    integer :: nocvAOStot
    !> Number of CABS AO orbitals
    integer :: ncabsAO
    !> Number of CABS MO orbitals
    integer :: ncabsMO
    !> Number of nvirt MO orbitals in AOS Space
    integer :: nvirtAOS
    !> Integral Orbital Type 
    Character, intent(in) :: intType(4) ! NB! Intent in because its read as a string!
    !> Integral Operator Type 
    Character, intent(in) :: intSpec(5) ! NB! Intent in because its read as a string!
    ! <n1,n2|INTSPEC|n3,n4> integrals stored in the order (n1,n2,n3,n4)
    real(realk), intent(inout) :: transformed_mo(:,:,:,:)

    !> MO trans coefficient for orbitals in <1,2|INTSPEC|3,4>
    type(ctype), dimension(4) :: C

    !> Dummy integer variables 
    integer :: i

    !> MO trans coefficient dimensions
    integer :: n11,n12,n21,n22,n31,n32,n41,n42

    !> MO coefficient matrix for the occupied EOS
    real(realk), target, intent(in) :: CoccEOS(:,:) !CoccEOS(nbasis,noccEOS)
    !> MO coefficient matrix for the occupied AOS
    real(realk), target, intent(in) :: CoccAOStot(:,:) !CoccEOS(nbasis,noccAOStot)
    !> MO coefficient matrix for the occupied + virtual EOS
    real(realk), target, intent(in) :: CocvAOStot(:,:) !CocvAOStot(nbasis, nocvAOS)
    !> MO coefficient matrix for the CABS 
    real(realk), target, intent(in) :: Ccabs(:,:) !Ccabs(ncabsAO, ncabsMO)
    !> MO coefficient matrix for the RI 
    real(realk), target, intent(in) :: Cri(:,:) !Cri(ncabsAO,ncabsAO)
    !> MO coefficient matrix for the Virtual AOS
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)
    integer :: offset

    nbasis   =  MyFragment%nbasis
    noccEOS  =  MyFragment%noccEOS
    noccAOStot  =  MyFragment%nocctot
    nvirtEOS = MyFragment%nvirtEOS
    nvirtAOS = MyFragment%nvirtAOS
    nocvAOStot =   noccAOStot + nvirtAOS
    ncabsAO = size(MyFragment%Ccabs,1)    
    ncabsMO = size(MyFragment%Ccabs,2)
    if(DECinfo%frozencore) then
       offset = MyFragment%ncore
    else
       offset = 0
    end if

    do i=1,4
       if(intType(i).EQ.'i') then ! occupied EOS
          C(i)%cmat => CoccEOS
          C(i)%n1 = nbasis
          C(i)%n2 = noccEOS               
       elseif(intType(i).EQ.'m') then ! all occupied (core+valence)
          C(i)%cmat => CoccAOStot
          C(i)%n1 = nbasis
          C(i)%n2 = noccAOStot 
       elseif(intType(i).EQ.'v') then 
          ! Frozen core: Only valence
          ! Not frozen core: All occupied orbitals (same as 'm')
          C(i)%cmat => CoccAOStot(1:nbasis,offset+1:noccAOStot)
          C(i)%n1 = nbasis
          C(i)%n2 = noccAOStot-offset 
       elseif(intType(i).EQ.'a') then ! virtual AOS
          C(i)%cmat => CvirtAOS
          C(i)%n1 = nbasis
          C(i)%n2 = nvirtAOS
       elseif(intType(i).EQ.'p') then !all occupied + virtual AOS
          C(i)%cmat => CocvAOStot
          C(i)%n1 = nbasis
          C(i)%n2 = nocvAOStot 
       elseif(intType(i).EQ.'c') then !cabs
          C(i)%cmat => Ccabs
          C(i)%n1 = ncabsAO
          C(i)%n2 = ncabsMO
       elseif(intType(i).EQ.'r') then !ri - MOs
          C(i)%cmat => Cri
          C(i)%n1 = ncabsAO
          C(i)%n2 = ncabsAO 
       endif
    enddo 

    !> Consistency check   
    if(size(transformed_mo,1) .NE. C(1)%n2) then
       print *, "Error: Wrong dim transformed_mo C(1)"
    end if

    if(size(transformed_mo,2) .NE. C(2)%n2) then
       print *, "Error: Wrong dim transformed_mo C(2)"
    end if

    if(size(transformed_mo,3) .NE. C(3)%n2) then
       print *, "Error: Wrong dim transformed_mo C(3)"
    end if

    if(size(transformed_mo,4) .NE. C(4)%n2) then
       print *, "Error: Wrong dim transformed_mo C(4)"
    end if

    call get_mp2f12_AO_transform_MO(MySetting,transformed_mo, C(1)%n1,C(1)%n2,C(2)%n1,C(2)%n2,C(3)%n1, &
         & C(3)%n2,C(4)%n1,C(4)%n2, C(1)%cmat,C(2)%cmat,C(3)%cmat,C(4)%cmat,intType,intSpec) 

  end subroutine get_mp2f12_MO


  !> Brief: Get <1,2|INTSPEC|3,4> MO integrals stored in the order (1,2,3,4).
  !> Author: Yang M. Wang
  !> Data: Nov 2013
  subroutine get_mp2f12_AO_transform_MO(MySetting,transformed_mo,n11,n12,n21,n22,n31,n32,n41,n42, &
       & C1,C2,C3,C4,INTTYPE,INTSPEC) 
    implicit none

    !> Integrals settings
    type(lssetting), intent(inout) :: Mysetting
    !> Integral Operator Type
    Character, intent(in) :: INTSPEC(5)
    !> Integral Orbital Type 
    Character, intent(in) :: INTTYPE(4)
    !> Orbital Type for Batching
    Character :: BatchType(4)
    !> MO trans coefficient dimensions
    integer,intent(in) :: n11,n12,n21,n22,n31,n32,n41,n42
    !> MO coefficients
    real(realk),intent(in),dimension(n11,n12) :: C1
    real(realk),intent(in),dimension(n21,n22) :: C2
    real(realk),intent(in),dimension(n31,n32) :: C3
    real(realk),intent(in),dimension(n41,n42) :: C4
    ! <n1,n2|INTSPEC|n3,n4> integrals stored in the order (n1,n2,n3,n4)
    real(realk), intent(inout) :: transformed_mo(n12,n22,n32,n42)  
    !> Dummy integral stored in the order (n3,n2,n4,n1)
    real(realk), pointer :: kjli(:,:,:,:)  
    !> Dummy MO coefficients
    real(realk), pointer :: C4T(:,:)  
    real(realk), pointer :: C3T(:,:)  
    real(realk), pointer :: C1T(:,:) 

    !> Variables for BATCH
    integer :: alphaB,gammaB,dimAlpha,dimGamma,GammaStart, GammaEnd, AlphaStart, AlphaEnd
    real(realk),pointer :: tmp1(:),tmp2(:),CoccEOST(:,:),CocvAOST(:,:)
    integer(kind=long) :: dim1,dim2
    integer :: i,m,k,n,idx,j,l
    logical :: FullRHS,doscreen
    integer :: MaxActualDimAlpha,nbatchesAlpha,nbatches,MaxActualDimGamma,nbatchesGamma,iorb
    type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
    type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
    integer :: iAO,nAObatches(4),AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
    logical :: MoTrans, NoSymmetry,SameMol
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
    type(batchtoorb), pointer :: batch2orbAlpha(:),batch2orbGamma(:)
#ifdef VAR_OMP
    integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
    TYPE(DECscreenITEM)   :: DecScreen

    integer :: MinAObatchSize, MaxAObatchSize, GammaBatchSize, AlphaBatchSize

    integer :: MinAlphaBatchSize, MinGammaBatchSize

    integer :: MaxdimAlpha, MaxdimGamma

    ! ****************************************************************
    ! Allocate mem space for a temporary array that will be reordered
    ! ****************************************************************
    call mem_alloc(kjli,n32,n22,n42,n12)

    ! ****************************************************************
    ! Allocate mem space for a temporary array that will be reordered
    ! ****************************************************************
    call mem_alloc(C4T,n42,n41)
    call mem_alloc(C3T,n32,n31)
    call mem_alloc(C1T,n12,n11)
    call mat_transpose(n41,n42, 1.0E0_realk,C4, 0.0E0_realk,C4T)
    call mat_transpose(n31,n32, 1.0E0_realk,C3, 0.0E0_realk,C3T)
    call mat_transpose(n11,n12, 1.0E0_realk,C1, 0.0E0_realk,C1T)

    ! ***********************************
    ! Determine batch Types ('R' or 'C')
    ! ***********************************
    do i=1,4
       if(intType(i).EQ.'i') then !occupied active
          BatchType(i) = 'R'          
       elseif(intType(i).EQ.'m') then !all occupied
          BatchType(i) = 'R'
       elseif(intType(i).EQ.'a') then !all virtual
          BatchType(i) = 'R'
       elseif(intType(i).EQ.'p') then !all occupied + virtual
          BatchType(i) = 'R'
       elseif(intType(i).EQ.'c') then !cabs
          BatchType(i) = 'C'          
       elseif(intType(i).EQ.'r') then !ri - MOs
          BatchType(i) = 'C'
       elseif(intType(i).EQ.'v') then !only valence
          BatchType(i) = 'R'
       else
          call lsquit('unknown option in get_mp2f12_AO_transform_MO',-1)
       endif
    enddo  
    !Determine MinGamma and MinAlpha    
    call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,MinGammaBatchSize,BatchType(3))
    call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,MinAlphaBatchSize,BatchType(1))
   
    ! ************************
    ! Determine AO batch sizes
    ! ************************
    ! NOTE: Ideally the batch sizes should be optimized according to the available memory
    ! (as is done e.g. in get_optimal_batch_sizes_for_mp2_integrals).
    ! For simplicity we simply choose the gamma batch to contain all basis functions,
    ! while we make the alpha batch as small as possible
    
    IF(DECinfo%useIchor)THEN
       iprint = 0           !print level for Ichor Integral code
       MoTrans = .FALSE.    !Do not transform to MO basis! 
       NoSymmetry = .FALSE. !Use Permutational Symmetry! 
       SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
       !Determine the full number of AO batches - not to be confused with the batches of AOs
       !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
       do i=1,4
          call determine_Ichor_nAObatches(mysetting,i,BatchType(i),nAObatches(i),DECinfo%output)
       enddo
       !Determine the minimum allowed AObatch size MinAObatch
       !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
       !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch = 6 (the 2*(Px,Py,Pz))
       !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
       !'R'  !Specifies that it is the Regular AO basis that should be batched
       !'C'  !Specifies that it is the CABS AO basis that should be batched
       iAO = 3 !the center that the batching should occur on.  
       call determine_MinimumAllowedAObatchSize(mysetting,iAO,BatchType(3),GammaBatchSize)
       iAO = 1 !the center that the batching should occur on.  
       call determine_MinimumAllowedAObatchSize(mysetting,iAO,BatchType(1),AlphaBatchSize)
    ELSE
       !> Minimum AO batch size
       call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,GammaBatchSize,BatchType(3))
       call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,AlphaBatchSize,BatchType(1))
    ENDIF
 
    !> Maximum AO batch size (all basis functions)
    !MaxAObatchSize = n31
    !> Setting MinAO to AO batch size for debug purposes
    !MinAObatchSize = n11

    !> Set alpha and gamma batch size as written above
    !GammaBatchSize = n31 ! Needs to be changed, For DEBUG purposes MaxAObatchSize
    !AlphaBatchSize = n11 ! Needs to be changes, For DEBUG purposes MinAObatchSize

    !MinAlphaBatchSize = AlphaBatchSize
    !MinGammaBatchSize = GammaBatchSize

    !if(DECinfo%F12DEBUG) then
    !   print *, "call get_max_batchsize..."
    !endif
    
    call get_max_batchsize(MaxdimAlpha,MaxdimGamma,MinAlphaBatchSize,MinGammaBatchSize,n11,n12,n21,n22,n31,n32,n41,n42)

    !if(DECinfo%F12DEBUG) then
   !    print *, "exit get_max_batchsize..."
   !    print *, "----------------------------"
   !    print *, "MaxdimAlpha: ", MaxdimAlpha
   !    print *, "MaxdimGamma: ", MaxdimGamma
   !    print *, "----------------------------"
   !    print *, "MindimAlpha: ", AlphaBatchSize
   !    print *, "MindimGamma: ", GammaBatchSize
   !    print *, "----------------------------"
   ! endif
        
    GammaBatchSize = MaxdimGamma
    AlphaBatchSize = MaxdimAlpha

    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
    IF(DECinfo%useIchor)THEN
       iAO = 3 !Gamma is the 3. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the GammaBatchSize, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mysetting,iAO,BatchType(iAO),GammaBatchSize,&
            & nbatchesGamma,DECinfo%output)
       call mem_alloc(AOGammabatchinfo,nbatchesGamma)
       !Construct the batches of AOS based on the GammaBatchSize, the requested
       !size of the AO batches - GammaBatchSize must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimGamma must be less og equal to GammaBatchSize
       call determine_Ichor_batchesofAOS(mysetting,iAO,BatchType(iAO),GammaBatchSize,&
            & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)
    ELSE
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchGamma,n31)
    
       call build_batchesofAOS(DECinfo%output,mysetting,GammaBatchSize,n31,MaxActualDimGamma,&
            & batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,orb2BatchGamma,BatchType(3))
       
       ! Batch to orbital information
       ! ----------------------------
       call mem_alloc(batch2orbGamma,nbatchesGamma)
       do idx=1,nbatchesGamma
          call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
          batch2orbGamma(idx)%orbindex = 0
          batch2orbGamma(idx)%norbindex = 0
       end do
       do iorb=1, n31
          idx = orb2batchGamma(iorb)
          batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
          K = batch2orbGamma(idx)%norbindex
          batch2orbGamma(idx)%orbindex(K) = iorb
       end do
    ENDIF

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************
    IF(DECinfo%useIchor)THEN
       iAO = 1 !Alpha is the 1. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the AlphaBatchSize, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mysetting,iAO,BatchType(iAO),AlphaBatchSize,&
            & nbatchesAlpha,DECinfo%output)
       call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
       !Construct the batches of AOS based on the AlphaBatchSize, the requested
       !size of the AO batches - AlphaBatchSize must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimAlpha must be less og equal to AlphaBatchSize
       call determine_Ichor_batchesofAOS(mysetting,iAO,BatchType(iAO),AlphaBatchSize,&
            & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
    ELSE
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchAlpha,n11)
       call build_batchesofAOS(DECinfo%output,mysetting,AlphaBatchSize,n11,&
            & MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,BatchType(1))

       ! Batch to orbital information
       ! ----------------------------
       call mem_alloc(batch2orbAlpha,nbatchesAlpha)
       do idx=1,nbatchesAlpha
          call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
          batch2orbAlpha(idx)%orbindex = 0
          batch2orbAlpha(idx)%norbindex = 0
       end do
       do iorb=1, n11
          idx = orb2batchAlpha(iorb)
          batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
          K = batch2orbAlpha(idx)%norbindex
          batch2orbAlpha(idx)%orbindex(K) = iorb
       end do
    ENDIF

    IF(DECinfo%useIchor)THEN
       !Calculate Screening integrals 
       call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mysetting,INTSPEC,SameMOL)
    ELSE
       ! Setting to FALSE for DEBUG purposes
       Mysetting%scheme%cs_screen = .FALSE.
       Mysetting%scheme%ps_screen = .FALSE.
       
       ! Integral screening stuff
       doscreen = Mysetting%scheme%cs_screen .or. Mysetting%scheme%ps_screen
       call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mysetting,&
            & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%integralthreshold)
       IF(doscreen)then
          call II_getBatchOrbitalScreen(DecScreen,mysetting,&
               & n31,nbatchesAlpha,nbatchesGamma,&
               & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
               & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
       endif
    ENDIF
    FullRHS = (nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)

#ifdef VAR_OMP
    if(DECinfo%PL>0) write(DECinfo%output,*) 'Starting VOVO integrals - OMP. Number of threads: ', &
         & OMP_GET_MAX_THREADS()
#else
    if(DECinfo%PL>0) write(DECinfo%output,*) 'Starting VOVO integrals - NO OMP!'
#endif

    ! ******************************************************************
    ! Start looping over gamma and alpha batches and calculate integrals
    ! ******************************************************************
    ! Zero output integrals to be on the safe side
    kjli = 0.0_realk

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
       IF(DECinfo%useIchor)THEN
          dimGamma = AOGammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
          GammaStart = AOGammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
          GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
          AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
          AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
       ELSE
          dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
          GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
          GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
       ENDIF
       BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
          IF(DECinfo%useIchor)THEN
             dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
             AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
             AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
             AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
             AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
          ELSE
             dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
             AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
             AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch
          ENDIF
          ! Get tmp1(beta(n21),delta(n41)|INTSPEC|alphaB(n11),gammaB(n31)) 
          ! ************************************************************************************
          dim1 = i8*n21*n41*dimAlpha*dimGamma   ! dimension for integral array tmp1
          call mem_alloc(tmp1,dim1)
          tmp1 = 0.0_realk

          IF(DECinfo%useIchor)THEN
             call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,mysetting,n21,n41,dimAlpha,dimGamma,&
                  & tmp1,INTSPEC,FULLRHS,1,nAObatches(2),1,nAObatches(4),AOAlphaStart,&
                  & AOAlphaEnd,AOGammaStart,AOGammaEnd,MoTrans,n21,n41,dimAlpha,dimGamma,NoSymmetry,DECinfo%integralthreshold)
          ELSE
             IF(doscreen) mysetting%LST_GAB_RHS => DECSCREEN%masterGabRHS
             IF(doscreen) mysetting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p
             call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
                  & mysetting, tmp1, batchindexAlpha(alphaB), batchindexGamma(gammaB), &
                  & batchsizeAlpha(alphaB), batchsizeGamma(gammaB), n21, n41, dimAlpha, dimGamma, FullRHS,&
                  & INTSPEC,DECinfo%integralthreshold)
          ENDIF
          ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)

          ! Transform beta(n21) to index "j" with C2(n21,j (n22))
          ! ***********************************************
          ! Note: ";" indicates the place where the array is transposed:
          ! tmp2(delta(n41),alphaB(n11),gammaB(n31),j) = 
          ! sum_{beta(n21)} tmp1^T(beta(n21);delta(n41),alphaB(n11),gammaB(n31)) * C2{beta(n21) j}
          m = n41*dimGamma*dimAlpha            ! first  dim of tmp1^T
          k = n21                              ! second dim of tmp1^T and first dim of C2
          n = n22                              ! second dim of C2
          dim2 = i8*n41*dimAlpha*dimGamma*n22  ! dim of tmp2 

          call mem_alloc(tmp2,dim2)
          call dec_simple_dgemm(m,k,n,tmp1,C2,tmp2, 't', 'n')
          call mem_dealloc(tmp1)

          ! Transform delta(n41) to index "l" with C4(n41,l)
          ! ************************************************
          ! tmp1(b,alphaB(n11),gammaB(n31),j) = 
          ! sum_{delta(n41)} C4^T(l,delta(n41)) tmp2(delta(n41),alphaB(n11),gammaB(n31),j)
          ! Note: We have stored the transposed C4^T matrix, so no need to transpose in
          ! the call to dgemm.

          m = n42                              ! first  dim of C4^T
          k = n41                              ! second dim of C4^T and first dim of tmp2
          n = dimAlpha*dimGamma*n22            ! second dim of tmp2 array
          dim1 = i8*n42*dimAlpha*dimGamma*n22  ! dim of tmp1 
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,C4T,tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2) 
                   
          ! Transpose to make alphaB(n11) and gammaB(n31) indices available
          ! ***************************************************************
          dim2 = dim1
          call mem_alloc(tmp2,dim2)
          ! tmp2(gammaB(n31), j, l, alphaB(n11) = tmp1^T(l, alphaB(n11); gammaB(n31), j)
          m = n42*dimAlpha       ! first  dim of tmp1 array
          n = n22*dimGamma       ! second dim of tmp1 array

          call mat_transpose(m, n, 1.0E0_realk, tmp1, 0.0E0_realk,tmp2)
          call mem_dealloc(tmp1)
          
          do i=1, dim2
             if (abs(tmp2(i)) > 10.0E-10_realk) then
                ! write(DECinfo%output,'(1X,a,1i4,f20.10)') 'tmp4:', i, tmp2(i)
             endif
          enddo
    
          ! Transform gammaB(n31) to index "k" with C3(n31,k)
          ! *************************************************
          ! tmp1(k,j,l,alphaB) = sum_{gammaBatch(n31) in gamma} C3T(k,gammaB) * tmp2(gammaB,j,l,alphaB)
          m = n32                          ! first  dim of C3T 
          k = dimGamma                     ! second dim of C3T and first dim of tmp2 
          n = n22*n42*dimAlpha             ! second dim of tmp2 
          dim1 = i8*n32*n22*n42*dimAlpha   ! dim of tmp1
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,C3T(:,GammaStart:GammaEnd),tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2)

          do i=1, dim1
             if (abs(tmp1(i)) > 10.0E-10_realk) then
               ! write(DECinfo%output,'(1X,a,1i4,f20.10)') 'tmp5:', i, tmp1(i)
             endif
          enddo

          ! Transform alphaB(n11) to index "i" with C1(n11,i)
          ! ************************************************
          ! kjli(k,j,l,i) =+ sum_{alphaB(n11) in alpha} tmp1(k,j,l,alphaB(n11))  C1T^T(alphaB(n11),i)
          m = n32*n22*n42              ! first dim of  tmp1 
          k = dimAlpha                 ! second dim of tmp1 and first of C1T^T
          n = n12                      ! second dim of C1T^T
          dim2 = i8*n32*n22*n42*n12    ! dim of kjli
          call dec_simple_dgemm_update(m,k,n,tmp1,C1T(:,AlphaStart:AlphaEnd),kjli, 'n', 't')
          call mem_dealloc(tmp1)

          ! Note: To have things consecutive in memory it is better to pass C1T to the dgemm
          ! routine and then transpose (insted of passing C1 and not transpose).

       end do BatchAlpha
    end do BatchGamma

    !> Reorder tmp(k,j,l,i) -> tmp(i,j,k,l)
    call array_reorder_4d(1.0E0_realk,kjli,n32,n22,n42,n12,[4,2,1,3],0.0E0_realk,transformed_mo)

    ! Free and nullify stuff
    ! **********************
    IF(DECinfo%useIchor)THEN
       call FREE_SCREEN_ICHORERI()
       call mem_dealloc(AOGammabatchinfo)
       call mem_dealloc(AOAlphabatchinfo)
    ELSE
       nullify(mysetting%LST_GAB_LHS)
       nullify(mysetting%LST_GAB_RHS)
       call free_decscreen(DECSCREEN)
       call free_batch(orb2batchGamma, batchdimGamma, batchsizeGamma, batchindexGamma, batch2orbGamma, &
            & orb2batchAlpha, batchdimAlpha, batchsizeAlpha, batchindexAlpha, batch2orbAlpha, &
            & nbatchesGamma, nbatchesAlpha)
    ENDIF

    ! Free F12 related pointers
    call mem_dealloc(C4T)
    call mem_dealloc(C3T)
    call mem_dealloc(C1T)
    call mem_dealloc(kjli)

  end subroutine get_mp2f12_AO_transform_MO

  !> Brief: Get <1,2|INTSPEC|3,4> MO integrals wrapper.
  !> Author: Yang M. Wang
  !> Data: Nov 2013
  subroutine get_mp2f12_MO_PDM(MyFragment,MySetting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,INTTYPE,INTSPEC,transformed_mo)
    implicit none

    !> Atomic fragment to be determined  (NOT pair fragment)
    type(decfrag), intent(inout) :: MyFragment
    !> Integrals settings   
    type(lssetting), intent(inout) :: Mysetting
    !> Number of basis functions AO
    integer :: nbasis
    !> Number of occupied orbitals MO in EOS space
    integer :: noccEOS
    !> Number of occupied orbitals MO in AOS space
    integer :: noccAOS
    !> Number of virtupied (virtual) orbitals MO in EOS space
    integer :: nvirtEOS
    !> Number of occupied + virtual MO in AOS space 
    integer :: nocvAOS
    !> Number of CABS AO orbitals
    integer :: ncabsAO
    !> Number of CABS MO orbitals
    integer :: ncabsMO
    !> Number of nvirt MO orbitals in AOS Space
    integer :: nvirtAOS
    !> Integral Orbital Type 
    Character, intent(in) :: intType(4) ! NB! Intent in because its read as a string!
    !> Integral Operator Type 
    Character, intent(in) :: intSpec(5) ! NB! Intent in because its read as a string!
    ! <n1,n2|INTSPEC|n3,n4> integrals stored in the order (n1,n2,n3,n4)
    type(tensor), intent(inout) :: transformed_mo

    !> MO trans coefficient for orbitals in <1,2|INTSPEC|3,4>
    type(ctype), dimension(4) :: C

    !> Dummy integer variables 
    integer :: i

    !> MO trans coefficient dimensions
    integer :: n11,n12,n21,n22,n31,n32,n41,n42

    !> MO coefficient matrix for the occupied EOS
    real(realk), target, intent(in) :: CoccEOS(:,:) !CoccEOS(nbasis,noccEOS)
    !> MO coefficient matrix for the occupied AOS
    real(realk), target, intent(in) :: CoccAOS(:,:) !CoccEOS(nbasis,noccAOS)
    !> MO coefficient matrix for the occupied + virtual EOS
    real(realk), target, intent(in) :: CocvAOS(:,:) !CocvAOS(nbasis, nocvAOS)
    !> MO coefficient matrix for the CABS 
    real(realk), target, intent(in) :: Ccabs(:,:) !Ccabs(ncabsAO, ncabsMO)
    !> MO coefficient matrix for the RI 
    real(realk), target, intent(in) :: Cri(:,:) !Cri(ncabsAO,ncabsAO)
    !> MO coefficient matrix for the Virtual AOS
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    nbasis   =  MyFragment%nbasis
    noccEOS  =  MyFragment%noccEOS
    noccAOS  =  MyFragment%noccAOS
    nvirtEOS = MyFragment%nvirtEOS
    nvirtAOS = MyFragment%nvirtAOS
    nocvAOS =   MyFragment%noccAOS + MyFragment%nvirtAOS
    ncabsAO = size(MyFragment%Ccabs,1)    
    ncabsMO = size(MyFragment%Ccabs,2)

    do i=1,4
       if(intType(i).EQ.'i') then ! occupied EOS
          C(i)%cmat => CoccEOS
          C(i)%n1 = nbasis
          C(i)%n2 = noccEOS               
       elseif(intType(i).EQ.'m') then ! occupied AOS
          C(i)%cmat => CoccAOS
          C(i)%n1 = nbasis
          C(i)%n2 = noccAOS 
       elseif(intType(i).EQ.'a') then ! virtual AOS
          C(i)%cmat => CvirtAOS
          C(i)%n1 = nbasis
          C(i)%n2 = nvirtAOS
       elseif(intType(i).EQ.'p') then !all occupied + virtual AOS
          C(i)%cmat => CocvAOS
          C(i)%n1 = nbasis
          C(i)%n2 = nocvAOS 
       elseif(intType(i).EQ.'c') then !cabs
          C(i)%cmat => Ccabs
          C(i)%n1 = ncabsAO
          C(i)%n2 = ncabsMO
       elseif(intType(i).EQ.'r') then !ri - MOs
          C(i)%cmat => Cri
          C(i)%n1 = ncabsAO
          C(i)%n2 = ncabsAO 
       endif
    enddo 

    !> Consistency check   
    if(transformed_mo%dims(1) .NE. C(1)%n2) then
       print *, "Error: Wrong dim transformed_mo C(1)"
    end if

    if(transformed_mo%dims(2) .NE. C(2)%n2) then
       print *, "Error: Wrong dim transformed_mo C(2)"
    end if

    if(transformed_mo%dims(3) .NE. C(3)%n2) then
       print *, "Error: Wrong dim transformed_mo C(3)"
    end if

    if(transformed_mo%dims(4) .NE. C(4)%n2) then
       print *, "Error: Wrong dim transformed_mo C(4)"
    end if

    call get_mp2f12_AO_transform_MO_PDM(MySetting,transformed_mo, C(1)%n1,C(1)%n2,C(2)%n1,C(2)%n2,C(3)%n1, &
         & C(3)%n2,C(4)%n1,C(4)%n2, C(1)%cmat,C(2)%cmat,C(3)%cmat,C(4)%cmat,intType,intSpec) 

  end subroutine get_mp2f12_MO_PDM


  !> Brief: Get <1,2|INTSPEC|3,4> MO integrals stored in the order (1,2,3,4).
  !> Author: Yang M. Wang
  !> Data: Nov 2013
  subroutine get_mp2f12_AO_transform_MO_PDM(MySetting,transformed_mo,n11,n12,n21,n22,n31,n32,n41,n42, &
       & C1,C2,C3,C4,INTTYPE,INTSPEC) 
    implicit none

    !> Integrals settings
    type(lssetting), intent(inout) :: Mysetting
    !> Integral Operator Type
    Character, intent(in) :: INTSPEC(5)
    !> Integral Orbital Type 
    Character, intent(in) :: INTTYPE(4)
    !> Orbital Type for Batching
    Character :: BatchType(4)
    !> MO trans coefficient dimensions
    integer,intent(in) :: n11,n12,n21,n22,n31,n32,n41,n42
    !> MO coefficients
    real(realk),intent(in),dimension(n11,n12) :: C1
    real(realk),intent(in),dimension(n21,n22) :: C2
    real(realk),intent(in),dimension(n31,n32) :: C3
    real(realk),intent(in),dimension(n41,n42) :: C4
    ! <n1,n2|INTSPEC|n3,n4> integrals stored in the order (n1,n2,n3,n4)
    type(tensor), intent(inout) :: transformed_mo
    !> Dummy integral stored in the order (n3,n2,n4,n1)
    real(realk), pointer :: kjli(:,:,:,:)  
    !> Dummy MO coefficients
    real(realk), pointer :: C4T(:,:)  
    real(realk), pointer :: C3T(:,:)  
    real(realk), pointer :: C1T(:,:) 

    !> Variables for BATCH
    integer :: alphaB,gammaB,dimAlpha,dimGamma,GammaStart, GammaEnd, AlphaStart, AlphaEnd
    real(realk),pointer :: tmp1(:),tmp2(:),CoccEOST(:,:),CocvAOST(:,:)
    integer(kind=long) :: dim1,dim2
    integer :: i,m,k,n,idx,j,l
    logical :: FullRHS,doscreen
    integer :: MaxActualDimAlpha,nbatchesAlpha,nbatches,MaxActualDimGamma,nbatchesGamma,iorb
    type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
    type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
    integer :: iAO,nAObatches(4),AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
    logical :: MoTrans, NoSymmetry,SameMol
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
    type(batchtoorb), pointer :: batch2orbAlpha(:),batch2orbGamma(:)
#ifdef VAR_OMP
    integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
    TYPE(DECscreenITEM)   :: DecScreen

    integer :: MinAObatchSize, MaxAObatchSize, GammaBatchSize, AlphaBatchSize

    integer :: MinAlphaBatchSize, MinGammaBatchSize

    integer :: MaxdimAlpha, MaxdimGamma

    ! ****************************************************************
    ! Allocate mem space for a temporary array that will be reordered
    ! ****************************************************************
    call mem_alloc(kjli,n32,n22,n42,n12)

    ! ****************************************************************
    ! Allocate mem space for a temporary array that will be reordered
    ! ****************************************************************
    call mem_alloc(C4T,n42,n41)
    call mem_alloc(C3T,n32,n31)
    call mem_alloc(C1T,n12,n11)
    call mat_transpose(n41,n42, 1.0E0_realk,C4, 0.0E0_realk,C4T)
    call mat_transpose(n31,n32, 1.0E0_realk,C3, 0.0E0_realk,C3T)
    call mat_transpose(n11,n12, 1.0E0_realk,C1, 0.0E0_realk,C1T)


    ! ************************
    ! Determine AO batch sizes
    ! ************************
    ! NOTE: Ideally the batch sizes should be optimized according to the available memory
    ! (as is done e.g. in get_optimal_batch_sizes_for_mp2_integrals).
    ! For simplicity we simply choose the gamma batch to contain all basis functions,
    ! while we make the alpha batch as small as possible
    
    ! ***********************************
    ! Determine batch Types ('R' or 'C')
    ! ***********************************
    do i=1,4
       if(intType(i).EQ.'i') then !occupied active
          BatchType(i) = 'R'          
       elseif(intType(i).EQ.'m') then !all occupied
          BatchType(i) = 'R'
       elseif(intType(i).EQ.'a') then !all virtual
          BatchType(i) = 'R'
       elseif(intType(i).EQ.'p') then !all occupied + virtual
          BatchType(i) = 'R'
       elseif(intType(i).EQ.'c') then !cabs
          BatchType(i) = 'C'          
       elseif(intType(i).EQ.'r') then !ri - MOs
          BatchType(i) = 'C'
       endif
    enddo

    !For debugging purposes    
    call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,MinGammaBatchSize,BatchType(3))
    call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,MinAlphaBatchSize,BatchType(1))

    IF(DECinfo%useIchor)THEN
       iprint = 0           !print level for Ichor Integral code
       MoTrans = .FALSE.    !Do not transform to MO basis! 
       NoSymmetry = .FALSE. !Use Permutational Symmetry! 
       SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
       !Determine the full number of AO batches - not to be confused with the batches of AOs
       !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
       do i=1,4
          call determine_Ichor_nAObatches(mysetting,i,BatchType(i),nAObatches(i),DECinfo%output)
       enddo
       !Determine the minimum allowed AObatch size MinAObatch
       !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
       !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch = 6 (the 2*(Px,Py,Pz))
       !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
       !'R'  !Specifies that it is the Regular AO basis that should be batched
       !'C'  !Specifies that it is the CABS AO basis that should be batched
       iAO = 3 !the center that the batching should occur on.  
       call determine_MinimumAllowedAObatchSize(mysetting,iAO,BatchType(3),GammaBatchSize)
       iAO = 1 !the center that the batching should occur on.  
       call determine_MinimumAllowedAObatchSize(mysetting,iAO,BatchType(1),AlphaBatchSize)
    ELSE
       !> Minimum AO batch size
       call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,GammaBatchSize,BatchType(3))
       call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,AlphaBatchSize,BatchType(1))
    ENDIF
 
    !> Maximum AO batch size (all basis functions)
    !MaxAObatchSize = n31
    !> Setting MinAO to AO batch size for debug purposes
    !MinAObatchSize = n11

    !> Set alpha and gamma batch size as written above
    !GammaBatchSize = n31 ! Needs to be changed, For DEBUG purposes MaxAObatchSize
    !AlphaBatchSize = n11 ! Needs to be changes, For DEBUG purposes MinAObatchSize

    !MinAlphaBatchSize = AlphaBatchSize
    !MinGammaBatchSize = GammaBatchSize
    
    call get_max_batchsize(MaxdimAlpha,MaxdimGamma,MinAlphaBatchSize,MinGammaBatchSize,n11,n12,n21,n22,n31,n32,n41,n42)

    GammaBatchSize = MaxdimGamma
    AlphaBatchSize = MaxdimAlpha

    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
    IF(DECinfo%useIchor)THEN
       iAO = 3 !Gamma is the 3. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the GammaBatchSize, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mysetting,iAO,BatchType(iAO),GammaBatchSize,&
            & nbatchesGamma,DECinfo%output)
       call mem_alloc(AOGammabatchinfo,nbatchesGamma)
       !Construct the batches of AOS based on the GammaBatchSize, the requested
       !size of the AO batches - GammaBatchSize must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimGamma must be less og equal to GammaBatchSize
       call determine_Ichor_batchesofAOS(mysetting,iAO,BatchType(iAO),GammaBatchSize,&
            & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)
    ELSE
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchGamma,n31)
       
       call build_batchesofAOS(DECinfo%output,mysetting,GammaBatchSize,n31,MaxActualDimGamma,&
            & batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,orb2BatchGamma,BatchType(3))
       
       ! Batch to orbital information
       ! ----------------------------
       call mem_alloc(batch2orbGamma,nbatchesGamma)
       do idx=1,nbatchesGamma
          call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
          batch2orbGamma(idx)%orbindex = 0
          batch2orbGamma(idx)%norbindex = 0
       end do
       do iorb=1, n31
          idx = orb2batchGamma(iorb)
          batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
          K = batch2orbGamma(idx)%norbindex
          batch2orbGamma(idx)%orbindex(K) = iorb
       end do
    ENDIF

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************
    IF(DECinfo%useIchor)THEN
       iAO = 1 !Alpha is the 1. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the AlphaBatchSize, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mysetting,iAO,BatchType(iAO),AlphaBatchSize,&
            & nbatchesAlpha,DECinfo%output)
       call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
       !Construct the batches of AOS based on the AlphaBatchSize, the requested
       !size of the AO batches - AlphaBatchSize must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimAlpha must be less og equal to AlphaBatchSize
       call determine_Ichor_batchesofAOS(mysetting,iAO,BatchType(iAO),AlphaBatchSize,&
            & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
    ELSE
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchAlpha,n11)
       call build_batchesofAOS(DECinfo%output,mysetting,AlphaBatchSize,n11,&
            & MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,BatchType(1))

       ! Batch to orbital information
       ! ----------------------------
       call mem_alloc(batch2orbAlpha,nbatchesAlpha)
       do idx=1,nbatchesAlpha
          call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
          batch2orbAlpha(idx)%orbindex = 0
          batch2orbAlpha(idx)%norbindex = 0
       end do
       do iorb=1, n11
          idx = orb2batchAlpha(iorb)
          batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
          K = batch2orbAlpha(idx)%norbindex
          batch2orbAlpha(idx)%orbindex(K) = iorb
       end do
    ENDIF

    IF(DECinfo%useIchor)THEN
       !Calculate Screening integrals 
       call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mysetting,INTSPEC,SameMOL)
    ELSE
       ! Setting to FALSE for DEBUG purposes
       Mysetting%scheme%cs_screen = .FALSE.
       Mysetting%scheme%ps_screen = .FALSE.
       
       ! Integral screening stuff
       doscreen = Mysetting%scheme%cs_screen .or. Mysetting%scheme%ps_screen
       call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mysetting,&
            & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%integralthreshold)
       IF(doscreen)then
          call II_getBatchOrbitalScreen(DecScreen,mysetting,&
               & n31,nbatchesAlpha,nbatchesGamma,&
               & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
               & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
       endif
    ENDIF
    FullRHS = (nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)

#ifdef VAR_OMP
    if(DECinfo%PL>0) write(DECinfo%output,*) 'Starting VOVO integrals - OMP. Number of threads: ', &
         & OMP_GET_MAX_THREADS()
#else
    if(DECinfo%PL>0) write(DECinfo%output,*) 'Starting VOVO integrals - NO OMP!'
#endif

    ! ******************************************************************
    ! Start looping over gamma and alpha batches and calculate integrals
    ! ******************************************************************
    ! Zero output integrals to be on the safe side
    kjli = 0.0_realk

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
       IF(DECinfo%useIchor)THEN
          dimGamma = AOGammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
          GammaStart = AOGammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
          GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
          AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
          AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
       ELSE
          dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
          GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
          GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
       ENDIF
       BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
          IF(DECinfo%useIchor)THEN
             dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
             AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
             AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
             AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
             AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
          ELSE
             dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
             AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
             AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch
          ENDIF
          ! Get tmp1(beta(n21),delta(n41)|INTSPEC|alphaB(n11),gammaB(n31)) 
          ! ************************************************************************************
          dim1 = i8*n21*n41*dimAlpha*dimGamma   ! dimension for integral array tmp1
          call mem_alloc(tmp1,dim1)
          tmp1 = 0.0_realk

          IF(DECinfo%useIchor)THEN
             call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,mysetting,n21,n41,dimAlpha,dimGamma,&
                  & tmp1,INTSPEC,FULLRHS,1,nAObatches(2),1,nAObatches(4),AOAlphaStart,&
                  & AOAlphaEnd,AOGammaStart,AOGammaEnd,MoTrans,n21,n41,dimAlpha,dimGamma,NoSymmetry,DECinfo%integralthreshold)
          ELSE
             IF(doscreen) mysetting%LST_GAB_RHS => DECSCREEN%masterGabRHS
             IF(doscreen) mysetting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p

             call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
                  & mysetting, tmp1, batchindexAlpha(alphaB), batchindexGamma(gammaB), &
                  & batchsizeAlpha(alphaB), batchsizeGamma(gammaB), n21, n41, dimAlpha, dimGamma, FullRHS,&
                  & INTSPEC,DECinfo%integralthreshold)
          ENDIF
          ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)

          ! Transform beta(n21) to index "j" with C2(n21,j (n22))
          ! ***********************************************
          ! Note: ";" indicates the place where the array is transposed:
          ! tmp2(delta(n41),alphaB(n11),gammaB(n31),j) = 
          ! sum_{beta(n21)} tmp1^T(beta(n21);delta(n41),alphaB(n11),gammaB(n31)) * C2{beta(n21) j}
          m = n41*dimGamma*dimAlpha            ! first  dim of tmp1^T
          k = n21                              ! second dim of tmp1^T and first dim of C2
          n = n22                              ! second dim of C2
          dim2 = i8*n41*dimAlpha*dimGamma*n22  ! dim of tmp2 

          call mem_alloc(tmp2,dim2)
          call dec_simple_dgemm(m,k,n,tmp1,C2,tmp2, 't', 'n')
          call mem_dealloc(tmp1)

          ! Transform delta(n41) to index "l" with C4(n41,l)
          ! ************************************************
          ! tmp1(b,alphaB(n11),gammaB(n31),j) = 
          ! sum_{delta(n41)} C4^T(l,delta(n41)) tmp2(delta(n41),alphaB(n11),gammaB(n31),j)
          ! Note: We have stored the transposed C4^T matrix, so no need to transpose in
          ! the call to dgemm.

          m = n42                              ! first  dim of C4^T
          k = n41                              ! second dim of C4^T and first dim of tmp2
          n = dimAlpha*dimGamma*n22            ! second dim of tmp2 array
          dim1 = i8*n42*dimAlpha*dimGamma*n22  ! dim of tmp1 
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,C4T,tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2) 
                   
          ! Transpose to make alphaB(n11) and gammaB(n31) indices available
          ! ***************************************************************
          dim2 = dim1
          call mem_alloc(tmp2,dim2)
          ! tmp2(gammaB(n31), j, l, alphaB(n11) = tmp1^T(l, alphaB(n11); gammaB(n31), j)
          m = n42*dimAlpha       ! first  dim of tmp1 array
          n = n22*dimGamma       ! second dim of tmp1 array

          call mat_transpose(m, n, 1.0E0_realk, tmp1, 0.0E0_realk,tmp2)
          call mem_dealloc(tmp1)
          
          do i=1, dim2
             if (abs(tmp2(i)) > 10.0E-10_realk) then
                ! write(DECinfo%output,'(1X,a,1i4,f20.10)') 'tmp4:', i, tmp2(i)
             endif
          enddo
    
          ! Transform gammaB(n31) to index "k" with C3(n31,k)
          ! *************************************************
          ! tmp1(k,j,l,alphaB) = sum_{gammaBatch(n31) in gamma} C3T(k,gammaB) * tmp2(gammaB,j,l,alphaB)
          m = n32                          ! first  dim of C3T 
          k = dimGamma                     ! second dim of C3T and first dim of tmp2 
          n = n22*n42*dimAlpha             ! second dim of tmp2 
          dim1 = i8*n32*n22*n42*dimAlpha   ! dim of tmp1
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,C3T(:,GammaStart:GammaEnd),tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2)

          do i=1, dim1
             if (abs(tmp1(i)) > 10.0E-10_realk) then
               ! write(DECinfo%output,'(1X,a,1i4,f20.10)') 'tmp5:', i, tmp1(i)
             endif
          enddo

          ! Transform alphaB(n11) to index "i" with C1(n11,i)
          ! ************************************************
          ! kjli(k,j,l,i) =+ sum_{alphaB(n11) in alpha} tmp1(k,j,l,alphaB(n11))  C1T^T(alphaB(n11),i)
          m = n32*n22*n42              ! first dim of  tmp1 
          k = dimAlpha                 ! second dim of tmp1 and first of C1T^T
          n = n12                      ! second dim of C1T^T
          dim2 = i8*n32*n22*n42*n12    ! dim of kjli
          call dec_simple_dgemm_update(m,k,n,tmp1,C1T(:,AlphaStart:AlphaEnd),kjli, 'n', 't')
          call mem_dealloc(tmp1)

          ! Note: To have things consecutive in memory it is better to pass C1T to the dgemm
          ! routine and then transpose (insted of passing C1 and not transpose).

       end do BatchAlpha
    end do BatchGamma

    !> Reorder tmp(k,j,l,i) -> tmp(i,j,k,l)
    call array_reorder_4d(1.0E0_realk,kjli,n32,n22,n42,n12,[4,2,1,3],0.0E0_realk,transformed_mo%elm1)

    ! Free and nullify stuff
    ! **********************
    IF(DECinfo%useIchor)THEN
       call FREE_SCREEN_ICHORERI()
       call mem_dealloc(AOGammabatchinfo)
       call mem_dealloc(AOAlphabatchinfo)
    ELSE
       nullify(mysetting%LST_GAB_LHS)
       nullify(mysetting%LST_GAB_RHS)
       call free_decscreen(DECSCREEN)
       call free_batch(orb2batchGamma, batchdimGamma, batchsizeGamma, batchindexGamma, batch2orbGamma, &
            & orb2batchAlpha, batchdimAlpha, batchsizeAlpha, batchindexAlpha, batch2orbAlpha, &
            & nbatchesGamma, nbatchesAlpha)
    ENDIF

    ! Free F12 related pointers
    call mem_dealloc(C4T)
    call mem_dealloc(C3T)
    call mem_dealloc(C1T)
    call mem_dealloc(kjli)

  end subroutine get_mp2f12_AO_transform_MO_PDM

  subroutine free_batch(orb2batchGamma, batchdimGamma, batchsizeGamma, batchindexGamma, batch2orbGamma, &
       & orb2batchAlpha, batchdimAlpha, batchsizeAlpha, batchindexAlpha, batch2orbAlpha, nbatchesGamma, nbatchesAlpha)
    implicit none

    integer :: idx
    integer, intent(in) :: nbatchesAlpha,nbatchesGamma
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
    type(batchtoorb), pointer :: batch2orbAlpha(:),batch2orbGamma(:)

    ! Free gamma batch stuff
    call mem_dealloc(orb2batchGamma)
    call mem_dealloc(batchdimGamma)
    call mem_dealloc(batchsizeGamma)
    call mem_dealloc(batchindexGamma)
    do idx=1,nbatchesGamma
       call mem_dealloc(batch2orbGamma(idx)%orbindex)
    end do
    call mem_dealloc(batch2orbGamma)

    ! Free alpha batch stuff
    call mem_dealloc(orb2batchAlpha)
    call mem_dealloc(batchdimAlpha)
    call mem_dealloc(batchsizeAlpha)
    call mem_dealloc(batchindexAlpha)
    do idx=1,nbatchesAlpha
       call mem_dealloc(batch2orbAlpha(idx)%orbindex)
    end do
    call mem_dealloc(batch2orbAlpha)

  end subroutine free_batch

  subroutine get_4Center_F12_integrals(mylsitem,MyMolecule,nbasis,nocc,noccfull,nvirt,ncabsAO,&
       & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
       & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    integer :: nbasis,nocc,nvirt,noccfull,ncabsAO
    real(realk),pointer :: Ripjq(:,:,:,:)
    real(realk),pointer :: Fijkl(:,:,:,:)
    real(realk),pointer :: Tijkl(:,:,:,:)
    real(realk),pointer :: Rimjc(:,:,:,:)
    real(realk),pointer :: Dijkl(:,:,:,:)
    real(realk),pointer :: Tirjk(:,:,:,:)
    real(realk),pointer :: Tijkr(:,:,:,:)
    real(realk),pointer :: Gipjq(:,:,:,:)
    real(realk),pointer :: Gimjc(:,:,:,:)
    real(realk),pointer :: Girjs(:,:,:,:)
    real(realk),pointer :: Girjm(:,:,:,:)
    real(realk),pointer :: Grimj(:,:,:,:)
    real(realk),pointer :: Gipja(:,:,:,:)
    real(realk),pointer :: Gpiaj(:,:,:,:)
    real(realk),pointer :: Gicjm(:,:,:,:)
    real(realk),pointer :: Gcimj(:,:,:,:)
    real(realk),pointer :: Gcirj(:,:,:,:)
    real(realk),pointer :: Gciaj(:,:,:,:)
    real(realk),pointer :: Giajc(:,:,:,:)
    !
    real(realk),pointer :: gao(:,:,:,:)

    if( MyMolecule%mem_distributed )then
       call lsquit("ERROR(get_4Center_F12_integrals): this routine does not work&
       & with distributed arrays in the fullmolecule type, yet",-1)
    endif

    call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')

    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ipip',gAO,Ripjq)

    !Calculate the various Gaussian geminal integrals with four regular AO indeces
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ipip',gAO,Gipjq)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRF')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'iiii',gAO,Fijkl)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRD')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                        MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'iiii',gAO,Dijkl)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRR2')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'iiii',gAO,Tijkl)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ipia',gAO,Gipja)


    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'piai',gAO,Gpiaj)


    call mem_dealloc(gao)

    !Calculate the various Gaussian geminal integrals with RRRC
    call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRC2')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'iiir',gAO,Tijkr)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCC')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'imic',gAO,Rimjc)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCG')

    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'imic',gAO,Gimjc)
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'iaic',gAO,Giajc)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRC2')

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with RCRR
    call mem_alloc(gao,nbasis,ncabsAO,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRR2')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'irii',gAO,Tirjk)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'irim',gAO,Girjm)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'icim',gAO,Gicjm)

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with RCRC
    call mem_alloc(gao,nbasis,ncabsAO,nbasis,ncabsAO)

    !Calculate the various Gaussian geminal integrals with four regular AO indeces
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRCG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'irir',gAO,Girjs)

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with CRRR
    call mem_alloc(gao,ncabsAO,nbasis,nbasis,nbasis)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'rimi',gAO,Grimj)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'cimi',gAO,Gcimj)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ciai',gAO,Gciaj)

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with CRCR
    call mem_alloc(gao,ncabsAO,nbasis,ncabsAO,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRCRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ciri',gAO,Gcirj)
    call mem_dealloc(gao)

  end subroutine get_4Center_F12_integrals

  subroutine free_4Center_F12_integrals(&
       & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
       & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)
    implicit none
    real(realk),pointer :: Ripjq(:,:,:,:)
    real(realk),pointer :: Fijkl(:,:,:,:)
    real(realk),pointer :: Tijkl(:,:,:,:)
    real(realk),pointer :: Rimjc(:,:,:,:)
    real(realk),pointer :: Dijkl(:,:,:,:)
    real(realk),pointer :: Tirjk(:,:,:,:)
    real(realk),pointer :: Tijkr(:,:,:,:)
    real(realk),pointer :: Gipjq(:,:,:,:)
    real(realk),pointer :: Gimjc(:,:,:,:)
    real(realk),pointer :: Girjs(:,:,:,:)
    real(realk),pointer :: Girjm(:,:,:,:)
    real(realk),pointer :: Grimj(:,:,:,:)
    real(realk),pointer :: Gipja(:,:,:,:)
    real(realk),pointer :: Gpiaj(:,:,:,:)
    real(realk),pointer :: Gicjm(:,:,:,:)
    real(realk),pointer :: Gcimj(:,:,:,:)
    real(realk),pointer :: Gcirj(:,:,:,:)
    real(realk),pointer :: Gciaj(:,:,:,:)
    real(realk),pointer :: Giajc(:,:,:,:)
    call mem_dealloc(Ripjq)
    call mem_dealloc(Fijkl)
    call mem_dealloc(Tijkl)
    call mem_dealloc(Rimjc)
    call mem_dealloc(Dijkl)
    call mem_dealloc(Tirjk)
    call mem_dealloc(Tijkr)
    call mem_dealloc(Gipjq)
    call mem_dealloc(Gimjc)
    call mem_dealloc(Girjs)
    call mem_dealloc(Girjm)
    call mem_dealloc(Grimj)
    call mem_dealloc(Gipja)
    call mem_dealloc(Gpiaj)
    call mem_dealloc(Gicjm)
    call mem_dealloc(Gcimj)
    call mem_dealloc(Gcirj)
    call mem_dealloc(Gciaj)
    call mem_dealloc(Giajc)
  end subroutine free_4Center_F12_integrals

  subroutine get_4Center_MO_integrals(mylsitem,lupri,nbasis,nocc,noccfull,nvirt,&
       & Cocc,Cvirt,inputstring,gAO,gMO)
    implicit none
    integer :: nocc,noccfull,nvirt,nCabsAO,nCabs,nbasis
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    character(len=4) :: inputstring
    integer :: ndim2(4),ndim1(4)
    real(realk),pointer :: gAO(:,:,:,:)
    real(realk),pointer :: gMO(:,:,:,:) ,elms(:)
    type(matrix) :: CMO(4)
    real(realk),dimension(nbasis,noccfull),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk),dimension(nbasis,nvirt),intent(in) :: Cvirt
    type(matrix) :: CMO_cabs,CMO_ri
    real(realk),pointer :: tmp(:,:,:,:)
    real(realk),pointer :: tmp2(:,:,:,:)
    character :: string(4)
    logical :: doCABS,doRI
    integer :: i,lupri,offset

    ! Offset:   Frozen core    : ncore
    !           Not frozen core: 0
    offset = noccfull - nocc

    string(1) = inputstring(1:1)
    string(2) = inputstring(2:2)
    string(3) = inputstring(3:3)
    string(4) = inputstring(4:4)
    doCABS = .FALSE.
    do i=1,4
       if(string(i).EQ.'c')then !occupied active
          doCABS = .TRUE.
       endif
    enddo
    doRI = .FALSE.
    do i=1,4
       if(string(i).EQ.'r')then !occupied active
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
    do i=1,4
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
          call dcopy(ndim2(i)*ndim1(i),Cocc(1:nbasis,offset+1:noccfull),1,CMO(i)%elms,1)
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
    call mem_alloc(tmp,ndim2(1),ndim1(2),ndim1(3),ndim1(4))
    call ls_dzero(tmp,ndim2(1)*ndim1(2)*ndim1(3)*ndim1(4)) !replace with a dgemm !!
    call sub1(gao,tmp,CMO(1)%elms,ndim2,ndim1)             !replace with a dgemm !!

    call mem_alloc(tmp2,ndim2(1),ndim2(2),ndim1(3),ndim1(4))
    call ls_dzero(tmp2,ndim2(1)*ndim2(2)*ndim1(3)*ndim1(4))
    call sub2(tmp,tmp2,CMO(2)%elms,ndim2,ndim1)
    call mem_dealloc(tmp)

    call mem_alloc(tmp,ndim2(1),ndim2(2),ndim2(3),ndim1(4))
    call ls_dzero(tmp,ndim2(1)*ndim2(2)*ndim2(3)*ndim1(4))
    call sub3(tmp2,tmp,CMO(3)%elms,ndim2,ndim1)
    call mem_dealloc(tmp2)

    call mem_alloc(gMO,ndim2(1),ndim2(2),ndim2(3),ndim2(4))
    call ls_dzero(gMO,ndim2(1)*ndim2(2)*ndim2(3)*ndim2(4))
    call sub4(tmp,gMO,CMO(4)%elms,ndim2,ndim1)
    call mem_dealloc(tmp)
    do i=1,4
       call mat_free(CMO(i))
    enddo

  contains

    subroutine sub1(gao,tmp,elms,ndim2,ndim1)
      implicit none
      integer :: ndim1(4),ndim2(4)
      real(realk),intent(in) :: elms(ndim1(1),ndim2(1))
      real(realk),intent(in) :: gao(ndim1(1),ndim1(2),ndim1(3),ndim1(4))
      real(realk) :: tmp(ndim2(1),ndim1(2),ndim1(3),ndim1(4))
      integer :: a,b,c,d,p

      do d=1,ndim1(4)
         do c=1,ndim1(3)
            do b=1,ndim1(2)
               do a=1,ndim2(1)

                  do p=1,ndim1(1)
                     tmp(a,b,c,d) = tmp(a,b,c,d) + gao(p,b,c,d)*elms(p,a)
                  end do

               end do
            end do
         end do
      end do
    end subroutine sub1

    subroutine sub2(tmp,tmp2,elms,ndim2,ndim1)
      implicit none
      integer :: ndim1(4),ndim2(4)
      real(realk),intent(in) :: elms(ndim1(2),ndim2(2))
      real(realk),intent(in) :: tmp(ndim2(1),ndim1(2),ndim1(3),ndim1(4))
      real(realk) :: tmp2(ndim2(1),ndim2(2),ndim1(3),ndim1(4))
      integer :: a,b,c,d,p

      do d=1,ndim1(4)
         do c=1,ndim1(3)
            do b=1,ndim2(2)
               do a=1,ndim2(1)

                  do p=1,ndim1(2)
                     tmp2(a,b,c,d) = tmp2(a,b,c,d) + tmp(a,p,c,d)*elms(p,b)
                  end do

               end do
            end do
         end do
      end do
    end subroutine sub2

    subroutine sub3(tmp2,tmp,elms,ndim2,ndim1)
      implicit none
      integer :: ndim1(4),ndim2(4)
      real(realk),intent(in) :: elms(ndim1(3),ndim2(3))
      real(realk),intent(in) :: tmp2(ndim2(1),ndim2(2),ndim1(3),ndim1(4))
      real(realk) :: tmp(ndim2(1),ndim2(2),ndim2(3),ndim1(4))
      integer :: a,b,c,d,p

      do d=1,ndim1(4)
         do c=1,ndim2(3)
            do b=1,ndim2(2)
               do a=1,ndim2(1)

                  do p=1,ndim1(3)
                     tmp(a,b,c,d) = tmp(a,b,c,d) + tmp2(a,b,p,d)*elms(p,c)
                  end do

               end do
            end do
         end do
      end do
    end subroutine sub3

    subroutine sub4(tmp,gMO,elms,ndim2,ndim1)
      implicit none
      
      integer :: ndim1(4),ndim2(4)
      real(realk),intent(in) :: elms(ndim1(4),ndim2(4))
      real(realk),intent(in) :: tmp(ndim2(1),ndim2(2),ndim2(3),ndim1(4))
      real(realk) :: gMO(ndim2(1),ndim2(2),ndim2(3),ndim2(4))
      integer :: a,b,c,d,p

      do d=1,ndim2(4)
         do c=1,ndim2(3)
            do b=1,ndim2(2)
               do a=1,ndim2(1)

                  do p=1,ndim1(4)
                     gMO(a,b,c,d) = gMO(a,b,c,d) + tmp(a,b,c,p)*elms(p,d)
                  end do

               end do
            end do
         end do
      end do
    end subroutine sub4

  end subroutine get_4Center_MO_integrals

  subroutine get_max_batchsize(dimAlpha,dimGamma,minAlpha,minGamma,n11,n12,n21,n22,n31,n32,n41,n42)
    implicit none
    
    !> Optimal Gamma and Alpha Batches
    integer,intent(inout) :: dimAlpha, dimGamma
    integer, intent(in) :: n11, n12, n21, n22, n31, n32, n41, n42
    integer, intent(in) :: minAlpha, minGamma
  
    !> Batching routine initializations
    real(realk) :: MemAvailable, step1, step2, step3, step4, MAXstepmem   
    real(realk) :: dim1,dim2,dim3,dim4,dim5,Ngamma,Nalpha,frac,GB,MB,KB,UNIT

    integer :: i,j,k

    !   dimAlpha = n11 !Full nbasis
    !   dimGamma = n31 !Full nbasis

    GB = 1.0E-9_realk
    MB = 1.0E-6_realk
    KB = 1.0E-3_realk
    frac = 0.9E0_realk   
    Maxstepmem = 0.0E0_realk
    UNIT = GB

    dimAlpha = n11
    dimGamma = n31
      
    MemAvailable = 0.0E0_realk   
    call get_currently_available_memory(MemAvailable)
    MemAvailable = MemAvailable*1.0E9_realk !In bytes
  
 !   call get_maxstepmem(MAXstepmem,dimAlpha,dimGamma,n11,n12,n21,n22,n31,n32,n41,n42,UNIT)
 
    if(DECinfo%F12DEBUG) then
!!$       print *, "----------------------------------"
!!$       print *, " Inside get_max_batchsize summary "
!!$       print *, "----------------------------------"
!!$       print *, "MemAvailable: ", MemAvailable*UNIT
!!$       print *, "MAXstepmem: ", MAXstepmem*UNIT
!!$       print *, "n11: ", n11
!!$       print *, "n21: ", n21
!!$       print *, "n31: ", n31
!!$       print *, "n41: ", n41
!!$       print *, "----------------------------------"
!!$       print *, "n12: ", n12
!!$       print *, "n22: ", n22
!!$       print *, "n32: ", n32
!!$       print *, "n42: ", n42
!!$       print *, "----------------------------------"
    endif
    
!!$    if(DECinfo%F12DEBUG) then
!!$       print *, "call get_maxstepmem..."
!!$       print *, "dimAlpha: ", dimAlpha
!!$       print *, "dimGamma: ", dimGamma
!!$    endif
    
    dimAlpha = minAlpha
    
    gamma: do k = n31, minGamma,-1

       call get_maxstepmem(MAXstepmem,dimAlpha,k,n11,n12,n21,n22,n31,n32,n41,n42,UNIT)

       if(Maxstepmem < frac*MemAvailable) then
          dimGamma = k
!!$          if(DECinfo%F12DEBUG) then
!!$             print *, "inside gamma loop: "
!!$             print *, "dimGamma:", dimGamma
!!$          endif
          exit gamma   
       endif
    enddo gamma

    alpha: do k = minAlpha, n11
       call get_maxstepmem(MAXstepmem,k,dimGamma,n11,n12,n21,n22,n31,n32,n41,n42,UNIT)

       if(Maxstepmem > frac*MemAvailable) then
          dimAlpha = k-1
!!$          if(DECinfo%F12DEBUG) then
!!$             print *, "inside alpha loop: "
!!$             print *, "dimAlpha:", dimAlpha
!!$          endif
          exit alpha   
       endif
       
       if(k==n11) then
          dimAlpha = n11
       endif

    enddo alpha


  end subroutine get_max_batchsize

  subroutine get_maxstepmem(MAXstepmem,dimAlpha,dimGamma,n11,n12,n21,n22,n31,n32,n41,n42,UNIT)
    implicit none

    integer, intent(in) :: n11, n12, n21, n22, n31, n32, n41, n42
    integer, intent(in) :: dimAlpha, dimGamma
    real(realk), intent(inout) :: MAXstepmem
    real(realk), intent(in) :: UNIT  

    integer(8) :: dim1,dim2,dim3,dim4
    real(realk) :: step1,step2,step3,step4
    
    dim1 = i8*n21*n41*dimAlpha*dimGamma   ! dim of tmp1
    dim2 = i8*n41*dimAlpha*dimGamma*n22   ! dim of tmp2 
    dim3 = i8*n42*dimAlpha*dimGamma*n22   ! dim of tmp1     
    dim4 = i8*n42*dimAlpha*n32*n22        ! dim of tmp1 

    step1 = (dim1 + dim2)*8 
    step2 = (dim2 + dim3)*8
    step3 = step2
    step4 = (dim3 + dim4)*8

    MAXstepmem = MAX(step1,step2,step3,step4) 

!!$    if(DECinfo%F12DEBUG) then
!!$       print *, "inside get_maxstepmem..."
!!$       print *, "step1: ", step1*UNIT
!!$       print *, "step2: ", step2*UNIT 
!!$       print *, "step3: ", step3*UNIT
!!$       print *, "step4: ", step4*UNIT
!!$       print *, "MAXstepmem: ", Maxstepmem*UNIT
!!$    endif

  end subroutine get_maxstepmem
 
  subroutine get_ES2_from_dec_main(MyMolecule,MyLsitem,Dmat,ES2)
    implicit none
    
    type(fullmolecule),intent(in) :: MyMolecule
    type(lsitem), intent(inout) :: Mylsitem
    real(realk), intent(inout) :: ES2
    type(matrix), intent(in) :: Dmat
    
    integer :: nbasis,nocc,nvirt,noccfull,ncabsAO,ncabs
    
    !> Singles contribution
    type(matrix) :: Fic
    type(matrix) :: Fcd
    type(matrix) :: Fac

    !> Fock AO
    type(matrix) :: Fcc
    type(matrix) :: Frc
 
    !> Fock matrices of type real
    real(realk), pointer :: Fic_real(:,:)
    real(realk), pointer :: Fcd_real(:,:)
    real(realk), pointer :: Fac_real(:,:)

    !Need to build Cabs
    ! Init stuff
    ! **********
    nbasis = MyMolecule%nbasis
    nocc   = MyMolecule%nocc
    nvirt  = MyMolecule%nvirt

    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)
    noccfull = nocc

    if( MyMolecule%mem_distributed )then
       call lsquit("ERROR(get_ES2_from_dec_main): this routine does not work&
       & with distributed arrays in the fullmolecule type, yet",-1)
    endif
    
    !Fcd
    call mat_init(Fcc,ncabsAO,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Fcd,ncabs,ncabs)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'cc',Fcc,Fcd)
    call mat_free(Fcc)
    
    call mem_alloc(Fcd_real,ncabs,ncabs)
    call mat_to_full(Fcd,1.0E0_realk,Fcd_real)
    call mat_free(Fcd)
    
    !Fic
    call mat_init(Frc,nbasis,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Frc,Dmat,MyLsitem,'RCRRC')
    call mat_init(Fic,nocc,ncabs)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ic',Frc,Fic)
    call mat_free(Frc)

    call mem_alloc(Fic_real,nocc,ncabs)
    call mat_to_full(Fic,1.0E0_realk,Fic_real)
    call mat_free(Fic)
    
    ! Fac
    call mat_init(Frc,nbasis,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Frc,Dmat,MyLsitem,'RCRRC')
    call mat_init(Fac,nvirt,ncabs)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ac',Frc,Fac)
    call mat_free(Frc)

    call mem_alloc(Fac_real,nvirt,ncabs)
    call mat_to_full(Fac,1.0E0_realk,Fac_real)
    call mat_free(Fac)
    
    call get_ES2(ES2,Fic_real,MyMolecule%oofock%elm2,MyMolecule%vvfock%elm2, &
         & Fcd_real,Fac_real,nocc,nvirt,ncabs,ncabsAO)
    
    call mem_dealloc(Fcd_real)
    call mem_dealloc(Fic_real)
    call mem_dealloc(Fac_real)

  end subroutine get_ES2_from_dec_main
  
  subroutine get_ES2(ES2,Fic,Fij,Fab,Fcd,Fac,nocc,nvirt,ncabs,ncabsAO)
    implicit none
    
    real(realk), target, intent(in) :: Fcd(:,:)
    real(realk), target, intent(in) :: Fab(:,:)
    real(realk), target, intent(in) :: Fij(:,:)
    real(realk), target, intent(in) :: Fac(:,:)
    real(realk), target, intent(in) :: Fic(:,:)
    
    real(realk), intent(inout) :: ES2
    real(realk) :: tmp1,tmp2,tmp3,tmp4
    real(realk) :: gamma0, gamma1, beta1, beta0, sigma1, alpha1, alpha2
    real(realk) :: denom
    real(realk) :: Ecorr
    real(realk) :: Ecorr_old

    integer, intent(inout) :: nocc,ncabs,nvirt,ncabsAO
    integer :: i,j,a,c,iter
   
    !xk
    real(realk), pointer :: x1ia(:,:) 
    real(realk), pointer :: x1ic(:,:)
    real(realk), pointer :: x2ia(:,:) 
    real(realk), pointer :: x2ic(:,:)

    !rk 
    real(realk), pointer :: r1ia(:,:) 
    real(realk), pointer :: r1ic(:,:)
    real(realk), pointer :: r2ia(:,:) 
    real(realk), pointer :: r2ic(:,:)

    !vk
    real(realk), pointer :: v1ia(:,:) 
    real(realk), pointer :: v1ic(:,:)
    
    !pk
    real(realk), pointer :: p1ia(:,:) 
    real(realk), pointer :: p1ic(:,:) 
    real(realk), pointer :: p0ia(:,:) 
    real(realk), pointer :: p0ic(:,:)

    !vectors
    call mem_alloc(x1ia,nocc,nvirt)
    call mem_alloc(x1ic,nocc,ncabs)
    call mem_alloc(x2ia,nocc,nvirt)
    call mem_alloc(x2ic,nocc,ncabs)
    
    call mem_alloc(r1ia,nocc,nvirt)
    call mem_alloc(r1ic,nocc,ncabs)
    call mem_alloc(r2ia,nocc,nvirt)
    call mem_alloc(r2ic,nocc,ncabs)

    call mem_alloc(p1ia,nocc,nvirt)
    call mem_alloc(p1ic,nocc,ncabs)
    call mem_alloc(p0ia,nocc,nvirt)
    call mem_alloc(p0ic,nocc,ncabs)

    call mem_alloc(v1ia,nocc,nvirt)
    call mem_alloc(v1ic,nocc,ncabs)

    !Initial start
    Ecorr_old = 0.0E0_realk

    r1ia = 0.0E0_realk
    r1ic = -1.0E0_realk*Fic

    x1ia = 0.0E0_realk
    x1ic = 0.0E0_realk

    gamma1 = 0.0E0_realk
    gamma0 = 1.0E0_realk

    !pk0 
    p0ia = r1ia
    p0ic = r1ic

    !CG algorithm
    CG_loop: do iter=1,1000


       gamma1 = inner_prod(r1ia,r1ic,r1ia,r1ic,nocc,nvirt,ncabs)
       beta1 = gamma1/gamma0

       !Update pk1 vector - needs to be changed to type_matrix later
       p1ia = r1ia + beta1*p0ia
       p1ic = r1ic + beta1*p0ic

       call mat_prod(Fij,Fab,Fac,Fcd,p1ia,p1ic,nocc,nvirt,ncabs,v1ia,v1ic)

       sigma1 = inner_prod(p1ia,p1ic,v1ia,v1ic,nocc,nvirt,ncabs)
         
       alpha1 = gamma1/sigma1

       x2ia = x1ia + alpha1*p1ia 
       x2ic = x1ic + alpha1*p1ic

       r2ia = r1ia - alpha1*v1ia
       r2ic = r1ic - alpha1*v1ic

       !Update variables
       r1ia = r2ia
       r1ic = r2ic
    
       gamma0 = gamma1

       p0ia = p1ia
       p0ic = p1ic
       
       x1ia = x2ia
       x1ic = x2ic
       
       Ecorr = 0.0E0_realk
       do c=1,ncabs
          do i=1,nocc
             Ecorr = Ecorr + Fic(i,c)*x2ic(i,c)
          enddo
       enddo

       !print*, "iter Ecorr Ecorr_old: ", iter, Ecorr, Ecorr_old
       if(abs(Ecorr-Ecorr_old)<1.0E-9) then
          exit 
       endif

       !Update the energy
       Ecorr_old = Ecorr  

    enddo CG_loop

    !Energy summation
    Ecorr = 0.0E0_realk 
    do c=1,ncabs 
       do i=1,nocc 
          Ecorr = Ecorr + Fic(i,c)*x2ic(i,c) 
       enddo
    enddo

    ES2 = 2.0E0_realk*Ecorr

    !Twice due to spinorbitals
    print *,"Singles Energy Contribution: ", 2.0E0_realk*Ecorr   
    print *,"Number of iterations: ", iter

    call mem_dealloc(x1ic)
    call mem_dealloc(x1ia)
    call mem_dealloc(x2ic)
    call mem_dealloc(x2ia)
    
    call mem_dealloc(p1ic)
    call mem_dealloc(p1ia)
    call mem_dealloc(p0ic)
    call mem_dealloc(p0ia)

    call mem_dealloc(r1ic)
    call mem_dealloc(r1ia)
    call mem_dealloc(r2ic)
    call mem_dealloc(r2ia)
    
    call mem_dealloc(v1ic)
    call mem_dealloc(v1ia)
    
 end subroutine get_ES2

 function inner_prod(t1ia,t1ic,t2ia,t2ic,nocc,nvirt,ncabs) 
    implicit none
    integer, intent(inout) :: nocc,ncabs,nvirt
    real(realk), target, intent(in) :: t1ia(:,:)
    real(realk), target, intent(in) :: t1ic(:,:)
    real(realk), target, intent(in) :: t2ia(:,:)
    real(realk), target, intent(in) :: t2ic(:,:)
    real(realk) :: inner_prod
 
    integer i,a,c

    inner_prod = 0.0E0_realk

    do i=1,nocc
       do a=1,nvirt
          inner_prod = inner_prod + t1ia(i,a)*t2ia(i,a)
       enddo
    enddo

    do i=1,nocc
       do c=1,ncabs
          inner_prod = inner_prod + t1ic(i,c)*t2ic(i,c)
       enddo
    enddo

 end function inner_prod

 subroutine mat_prod(Fij,Fab,Fac,Fcd,t1ia,t1ic,nocc,nvirt,ncabs,t2ia,t2ic) 
    integer, intent(inout) :: nocc,ncabs,nvirt
    real(realk), target, intent(in) :: t1ia(:,:)
    real(realk), target, intent(in) :: t1ic(:,:)

    real(realk), target, intent(in) :: Fij(:,:)
    real(realk), target, intent(in) :: Fab(:,:)
    real(realk), target, intent(in) :: Fac(:,:)
    real(realk), target, intent(in) :: Fcd(:,:)

    real(realk), target, intent(inout) :: t2ia(:,:)
    real(realk), target, intent(inout) :: t2ic(:,:)

    integer i,j,a,b,c,d
    real(realk) :: tmp1,tmp2

    do i=1,nocc
       do a=1,nvirt
          tmp1 = 0.0E0_realk
          do b=1, nvirt 
             tmp1 = tmp1 + Fab(a,b)*t1ia(i,b)
          enddo
          do d=1, ncabs
             tmp1 = tmp1 + Fac(a,d)*t1ic(i,d)
          enddo
          do j=1, nocc
             tmp1 = tmp1 - Fij(i,j)*t1ia(j,a)
          enddo
          t2ia(i,a) = tmp1
       enddo
    enddo

    do i=1,nocc
       do c=1,ncabs
          tmp2 = 0.0E0_realk
          do b=1,nvirt
             tmp2 = tmp2 + Fac(b,c)*t1ia(i,b)
          enddo
          do d=1,ncabs
             tmp2 = tmp2 + Fcd(c,d)*t1ic(i,d)
          enddo
          do j=1,nocc
             tmp2 = tmp2 - Fij(i,j)*t1ic(j,c)
          enddo
          t2ic(i,c) = tmp2
       enddo
    enddo

 end subroutine mat_prod




 subroutine  dec_get_CABS_orbitals(molecule,mylsitem)
    implicit none

    !> Full molecule structure to be initialized
    type(fullmolecule), intent(inout) :: molecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem

    type(matrix) :: CMO_cabs
    integer :: ncabsAO,ncabs

    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)
    molecule%nCabsAO = ncabsAO
    molecule%nCabsMO = ncabs
    call mat_init(CMO_cabs,nCabsAO,nCabs)

    call init_cabs(DECinfo%full_molecular_cc)
    call build_CABS_MO(CMO_cabs,ncabsAO,mylsitem%SETTING,DECinfo%output)
    IF(.NOT.DECinfo%full_molecular_cc)call free_cabs()

    ! NB! Memory leak need to be freed somewhere
    !    call mem_alloc(molecule%Ccabs,ncabsAO,nCabs)
    !    call mat_to_full(CMO_cabs,1.0E0_realk,molecule%Ccabs)
    call mat_free(CMO_cabs)

 end subroutine dec_get_CABS_orbitals

  subroutine  dec_get_RI_orbitals(molecule,mylsitem)
    implicit none

    !> Full molecule structure to be initialized
    type(fullmolecule), intent(inout) :: molecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem

    type(matrix) :: CMO_RI
    integer :: ncabsAO,ncabs

    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)
    molecule%nCabsAO = ncabsAO
    molecule%nCabsMO = ncabs
    call mat_init(CMO_RI,ncabsAO,ncabsAO)

    call init_ri(DECinfo%full_molecular_cc)
    call build_RI_MO(CMO_RI,ncabsAO,mylsitem%SETTING,DECinfo%output)
    IF(.NOT.DECinfo%full_molecular_cc)call free_cabs()

    ! NB! Memory leak need to be freed somewhere
    !    call mem_alloc(molecule%Cri,ncabsAO,ncabsAO) 
    !    call mat_to_full(CMO_RI,1.0E0_realk,molecule%Cri)
    call mat_free(CMO_RI)

  end subroutine dec_get_RI_orbitals

!==========================================================
!===           Subroutines for DEC-RIMP2-F12            ===
!==========================================================

!Perform on GPU 
subroutine ContractOne4CenterF12IntegralsRI(nBA,n,Calpha,EJ,EK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n
   real(realk),intent(in)    :: Calpha(nBA,n,n)
   real(realk),intent(inout) :: EJ,EK
   !local variables
   integer :: I,ALPHA,J
   real(realk) :: TMP,TMPV(n),TMPI
   !Dopair
   logical,intent(in),optional :: dopair_occ_in(n,n)
   logical :: dopair_occ(n,n)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif 
   
   !Exchange Fiijj
   EK = 0.0E0_realk
   EJ = 0.0E0_realk

   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(I,J,&
   !$OMP ALPHA) SHARED(Calpha,n,nba,dopair_occ) REDUCTION(+:EK,EJ)
   DO I=1,n
      DO J=1,n
         if(dopair_occ(I,J)) then
            DO ALPHA = 1, NBA
               EJ = EJ + Calpha(ALPHA,I,I)*Calpha(ALPHA,J,J)
               EK = EK + Calpha(ALPHA,I,J)*Calpha(ALPHA,J,I)
            ENDDO
         endif
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !EJK = 7.0/32.0*EJ + 1.0/32.0*EK 
   print *,"COULOMBX2:  ",  -1.0E0_realk*(5.0E0_realk*0.25E0_realk)*EJ
   print *,"EXCHANGEX2: ",   1.0E0_realk*(0.250_realk*EK) 
   print *,"COULOMBX2+EXCHANGEX2:", -1.0E0_realk*((5.0E0_realk*0.25E0_realk)*EJ-EK*0.25E0_realk)   

end subroutine ContractOne4CenterF12IntegralsRI

!Perform on GPU 
subroutine ContractOne4CenterF12IntegralsRI2(nBA,n,Calpha,Fii,EJ,EK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n
   real(realk),intent(in)    :: Calpha(nBA,n,n)
   real(realk),intent(inout) :: EJ,EK
   real(realk),intent(IN) :: Fii(n,n)
   !local variables
   integer :: i,alpha,j
   real(realk) :: tmp1,tmp2,TMPV(n),TMPI
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n,n)                            
   logical :: dopair_occ(n,n)                                                     
   if(present(dopair_occ_in)) then                                                  
      dopair_occ = dopair_occ_in                                                
   else                                                                        
      dopair_occ = .TRUE.                                                     
   endif                 
   !Exchange Fiijj
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   DO i=1,n
      DO j=1,n
         if(dopair_occ(i,j)) then
            tmp1 = 0.0E0_realk
            tmp2 = 0.0E0_realk
            DO alpha = 1, nBA
               tmp1 = tmp1 + Calpha(ALPHA,i,i)*Calpha(ALPHA,j,j)
               tmp2 = tmp2 + Calpha(ALPHA,i,j)*Calpha(ALPHA,j,i)
            ENDDO
            EJ = EJ + tmp1*(Fii(i,i)+Fii(j,j))
            EK = EK + tmp2*(Fii(i,i)+Fii(j,j))
         endif
      ENDDO
   ENDDO

end subroutine ContractOne4CenterF12IntegralsRI2

!Perform on GPU 
subroutine ContractOne4CenterF12IntegralsRI2_nc(nBA,n,Calpha,CalphaT,EJ,EK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n
   real(realk),intent(in)    :: Calpha(nBA,n,n)
   real(realk),intent(in)    :: CalphaT(nBA,n,n)
   real(realk),intent(inout) :: EJ,EK
   !local variables
   integer :: i,alpha,beta,j
   real(realk) :: tmpR1,tmpR2,tmpG1,tmpG2,TMPV(n),TMPI
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n,n)                            
   logical :: dopair_occ(n,n)                                                     
 
   if(present(dopair_occ_in)) then                                                  
      dopair_occ = dopair_occ_in                                                
   else                                                                        
      dopair_occ = .TRUE.                                                     
   endif                 

   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   DO i=1,n
      DO j=1,n
         if(dopair_occ(i,j)) then
            tmpR1 = 0.0E0_realk
            tmpR2 = 0.0E0_realk
            DO alpha = 1, nBA
               tmpR1 = tmpR1 + Calpha(ALPHA,i,i)*CalphaT(ALPHA,j,j)
               tmpR2 = tmpR2 + Calpha(ALPHA,j,j)*CalphaT(ALPHA,i,i)
            ENDDO
            tmpG1 = 0.0E0_realk
            tmpG2 = 0.0E0_realk
            DO beta = 1, nBA
               tmpG1 = tmpG1 + Calpha(BETA,j,i)*CalphaT(BETA,i,j)
               tmpG2 = tmpG2 + Calpha(BETA,i,j)*CalphaT(BETA,j,i)
            ENDDO
            EJ = EJ + tmpR1+tmpR2
            EK = EK + tmpG1+tmpG2
         endif
      ENDDO
   ENDDO

end subroutine ContractOne4CenterF12IntegralsRI2_nc

subroutine ContractOne4CenterF12IntegralsRobustRI(nBA,n,nbasis,Rtilde,CalphaR,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n,nbasis
   real(realk),intent(in)    :: Rtilde(nBA,n,n)
   real(realk),intent(in)    :: CalphaR(nBA,n,nbasis)
   real(realk),intent(inout) :: EJK
   !local variables
   integer :: I,ALPHA,J
   real(realk) :: TMP,TMP_IJIJ,TMP_JIIJ,TMPD
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n,n)
   logical :: dopair_occ(n,n)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   TMPD = 0.0E0_realk
   TMP = 0.0E0_realk
   EJK = 0.0E0_realk
   !  EJK1 = 0.0E0_realk
   !  EJK2 = 0.0E0_realk
   !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J,ALPHA,TMP_IJIJ,TMP_JIIJ) SHARED(CalphaR,Rtilde,n,&
   !$OMP nba,dopair_occ) REDUCTION(+:TMP,TMPD)
   !Diagonal
   DO J=1,n
      IF(dopair_occ(J,J)) THEN
         TMP_IJIJ = 0.0E0_realk
         DO ALPHA = 1,NBA
            TMP_IJIJ = TMP_IJIJ + CalphaR(ALPHA,J,J)*Rtilde(ALPHA,J,J) + Rtilde(ALPHA,J,J)*CalphaR(ALPHA,J,J)
         ENDDO
         TMPD = TMPD + TMP_IJIJ
      ENDIF
      !Non Diagonal
      DO I=j+1,n
         IF(dopair_occ(I,J)) THEN
            TMP_IJIJ = 0.0E0_realk
            TMP_JIIJ = 0.0E0_realk
            DO ALPHA = 1,NBA
               TMP_IJIJ = TMP_IJIJ + CalphaR(ALPHA,I,I)*Rtilde(ALPHA,J,J) + Rtilde(ALPHA,I,I)*CalphaR(ALPHA,J,J)
               TMP_JIIJ = TMP_JIIJ + CalphaR(ALPHA,I,J)*Rtilde(ALPHA,J,I) + Rtilde(ALPHA,I,J)*CalphaR(ALPHA,J,I)
            ENDDO
            TMP = TMP + 7.0E0_realk * TMP_IJIJ + TMP_JIIJ
            !        EJK1 = EJK1 + 7.0E0_realk * TMP_IJIJ
            !        EJK2 = EJK2 + TMP_JIIJ
         ENDIF
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK = 0.25E0_realk*tmpD + 0.0625E0_realk*TMP
   !  print*,'A',0.25E0_realk*tmpD
   !  print*,'B',0.0625E0_realk*EJK1
   !  print*,'C',0.0625E0_realk*EJK2
end subroutine ContractOne4CenterF12IntegralsRobustRI

!Perform on GPU 
subroutine ContractOne4CenterF12IntegralsRIB23(nBA,n1,n2,CalphaR,CalphaG,HJir,Coeff,EJK2,EJK3,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n1,n1)
   real(realk),intent(inout) :: EJK2,EJK3
   real(realk),intent(in)    :: Coeff
   real(realk)               :: ED2,EJ2,EJ3,EK2,EK3
   real(realk),intent(IN)    :: HJir(n1,n2)
   !local variables
   integer :: c,i,j,alpha,beta
   real(realk) :: tmpR2,tmpR21,tmpR22,tmpR31,tmpR32,tmp,tmpR
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   ED2 = 0.0E0_realk
   EJ2 = 0.0E0_realk
   EK2 = 0.0E0_realk
   EJ3 = 0.0E0_realk
   EK3 = 0.0E0_realk 
   DO c=1,n2
      DO j=1,n1
         IF(dopair_occ(J,J)) THEN
            !Diagonal
            tmpR2 = 0.0E0_realk
            DO alpha = 1,nBA
               tmpR2 = tmpR2 + CalphaR(alpha,J,c)*CalphaG(alpha,J,J)
            ENDDO
            ED2 = ED2 + tmpR2*hJir(J,c)
         ENDIF
         !Non Diagonal
         DO i=j+1,n1
            IF(dopair_occ(I,J)) THEN
               tmpR21 = 0.0E0_realk
               tmpR22 = 0.0E0_realk
               tmpR31 = 0.0E0_realk
               tmpR32 = 0.0E0_realk
               DO alpha = 1,nBA
                  tmpR21 = tmpR21 + CalphaR(alpha,i,c)*CalphaG(alpha,j,j)
                  tmpR22 = tmpR22 + CalphaR(alpha,j,c)*CalphaG(alpha,i,j)

                  tmpR31 = tmpR31 + CalphaR(alpha,j,c)*CalphaG(alpha,i,i)
                  tmpR32 = tmpR32 + CalphaR(alpha,i,c)*CalphaG(alpha,i,j)
               ENDDO
               EJ2 = EJ2 + tmpR21*hJir(i,c)
               EK2 = EK2 + tmpR22*hJir(i,c)
               EJ3 = EJ3 + tmpR31*hJir(j,c)
               EK3 = EK3 + tmpR32*hJir(j,c)
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   EJK2 = ED2*1.0/4.0 + (EJ2*(7.0/16.0) + EK2*(1.0/16.0))
   EJK3 = ED2*1.0/4.0 + (EJ3*(7.0/16.0) + EK3*(1.0/16.0)) 

end subroutine ContractOne4CenterF12IntegralsRIB23

subroutine ContractTwo4CenterF12IntegralsRI_p(nBA,n1,n2,CalphaR,CalphaG,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   !local variables
   integer :: q,p,i,j,alpha,beta
   real(realk) :: tmpR,tmpG1,tmpG2,tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif  
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,p,q,tmpR,&
   !$OMP tmpG1,tmpG2) SHARED(CalphaR,CalphaG,n2,n1,&
   !$OMP nba,dopair_occ) REDUCTION(+:EJ,EK,ED)
   DO q=1,n2 !nocv
      DO p=1,n2 !nocv
         DO j=1,n1 !nocc
            DO i=1,n1 !nocc
               IF(dopair_occ(I,J)) THEN
                  tmpR = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR = tmpR + CalphaR(alpha,i,p)*CalphaR(alpha,j,q)
                  ENDDO
                  tmpG1 = 0.0E0_realk
                  tmpG2 = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG1 = tmpG1 + CalphaG(beta,i,p)*CalphaG(beta,j,q)
                     tmpG2 = tmpG2 + CalphaG(beta,j,p)*CalphaG(beta,i,q)
                  ENDDO
                  EJ = EJ + tmpR*tmpG1 
                  EK = EK + tmpR*tmpG2
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK = 1.25E0_realk*EJ - 0.25E0_realk*EK 
end subroutine ContractTwo4CenterF12IntegralsRI_p

subroutine ContractTwo4CenterF12IntegralsRI(nBA,n1,n2,CalphaR,CalphaG,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n2,n1)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   !local variables
   integer :: q,p,i,j,alpha,beta
   real(realk) :: tmpR,tmpG1,tmpG2,tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif  
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,p,q,tmpR,&
   !$OMP tmpG1,tmpG2) SHARED(CalphaR,CalphaG,n2,n1,&
   !$OMP nba,dopair_occ) REDUCTION(+:EJ,EK,ED)
   DO q=1,n2 !nocv
      DO p=1,n2 !nocv
         DO j=1,n1 !nocc
            !Diagonal
            IF(dopair_occ(J,J)) THEN
               tmpR = 0.0E0_realk
               DO alpha = 1,NBA
                  tmpR = tmpR + CalphaR(alpha,j,p)*CalphaR(alpha,j,q)
               ENDDO
               tmpG1 = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG1 = tmpG1 + CalphaG(beta,p,j)*CalphaG(beta,q,j)
               ENDDO
               ED = ED + tmpR*tmpG1
            ENDIF
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpR = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR = tmpR + CalphaR(alpha,i,p)*CalphaR(alpha,j,q)
                  ENDDO
                  tmpG1 = 0.0E0_realk
                  tmpG2 = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG1 = tmpG1 + CalphaG(beta,p,i)*CalphaG(beta,q,j)
                     tmpG2 = tmpG2 + CalphaG(beta,p,j)*CalphaG(beta,q,i)
                  ENDDO
                  EJ = EJ + tmpR*tmpG1 
                  EK = EK + tmpR*tmpG2
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK = ED + 2.5E0_realk*EJ - 0.50_realk*EK 
end subroutine ContractTwo4CenterF12IntegralsRI

subroutine ContractTwo4CenterF12IntegralsRIC(nBA,n1,n2,CalphaV,CalphaD,Taibj,EJK,dopair_occ_in)
   implicit none 
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaV(nBA,n1,n2),CalphaD(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED,EJ,EK
   real(realk),pointer       :: Caibj(:,:,:,:)
   real(realk),intent(in)    :: Taibj(:,:,:,:)
   !local variables
   integer :: a,b,i,j,alpha,beta
   real(realk) :: tmp,tmpR,tmpG
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO b=1,n2
      DO a=1,n2
         DO j=1,n1
            !Diagonal
            IF(dopair_occ(J,J)) THEN
               tmp = 0.0E0_realk
               DO alpha = 1,NBA
                  tmp = tmp + CalphaV(alpha,j,a)*CalphaD(alpha,j,b) + CalphaV(alpha,j,b)*CalphaD(alpha,j,a)
               ENDDO
               ED = ED + tmp*Taibj(a,j,b,j)
            ENDIF
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpR = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR = tmpR + CalphaV(alpha,i,a)*CalphaD(alpha,j,b) + CalphaV(alpha,j,b)*CalphaD(alpha,i,a)
                  ENDDO
                  tmpG = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG = tmpG + CalphaV(beta,j,a)*CalphaD(beta,i,b) + CalphaV(beta,i,b)*CalphaD(beta,j,a)
                  ENDDO
                  EJ = EJ + tmpR*Taibj(a,i,b,j)
                  EK = EK + tmpG*Taibj(a,i,b,j)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = - ED - 2.5E0_realk*EJ + 0.5E0_realk*EK 
end subroutine ContractTwo4CenterF12IntegralsRIC

subroutine ContractTwo4CenterF12IntegralsRIC_p(nBA,n1,n2,CalphaV,CalphaD,Taibj,EJK,dopair_occ_in)
   implicit none 
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaV(nBA,n1,n2),CalphaD(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED,EJ,EK
   real(realk),pointer       :: Caibj(:,:,:,:)
   real(realk),intent(in)    :: Taibj(:,:,:,:)
   !local variables
   integer :: a,b,i,j,alpha,beta
   real(realk) :: tmp,tmpR,tmpG
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO b=1,n2
      DO a=1,n2
         DO j=1,n1
            DO i=1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpR = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR = tmpR + CalphaV(alpha,i,a)*CalphaD(alpha,j,b) + CalphaV(alpha,j,b)*CalphaD(alpha,i,a)
                  ENDDO
                  tmpG = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG = tmpG + CalphaV(beta,j,a)*CalphaD(beta,i,b) + CalphaV(beta,i,b)*CalphaD(beta,j,a)
                  ENDDO
                  EJ = EJ + tmpR*Taibj(a,i,b,j)
                  EK = EK + tmpG*Taibj(a,i,b,j)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = -1.25E0_realk*EJ + 0.25E0_realk*EK 
end subroutine ContractTwo4CenterF12IntegralsRIC_p

subroutine ContractTwo4CenterF12IntegralsRIX(nBA,n1,n2,CalphaG,Fii,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   real(realk),intent(IN)    :: Fii(n1,n1)
   !local variables
   integer :: q,p,i,j,alpha,beta
   real(realk) :: tmpG1,tmpG2,tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO q=1,n2
      DO p=1,n2
         DO j=1,n1
            !Diagonal
            IF(dopair_occ(J,J)) THEN
               tmpG1 = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG1 = tmpG1 + CalphaG(beta,j,p)*CalphaG(beta,j,q)
                  !We have a factor 2 but it's integrated into the reduction
               ENDDO
               ED = ED + tmpG1*tmpG1*Fii(j,j)
            ENDIF
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpG1 = 0.0E0_realk
                  tmpG2 = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG1 = tmpG1 + CalphaG(beta,i,p)*CalphaG(beta,j,q)
                     tmpG2 = tmpG2 + CalphaG(beta,j,p)*CalphaG(beta,i,q)
                  ENDDO
                  EJ = EJ + tmpG1*tmpG1*(Fii(i,i) + Fii(j,j)) 
                  EK = EK + tmpG1*tmpG2*(Fii(i,i) + Fii(j,j))
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = (ED*(0.5E0_realk) + (7.0/16.0*EJ + 1.0/16.0*EK)) 
end subroutine ContractTwo4CenterF12IntegralsRIX

subroutine ContractTwo4CenterF12IntegralsRIX_nc(nBA,n1,n2,CalphaC,CalphaG,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaC(nBA,n2,n1)
   real(realk),intent(in)    :: CalphaG(nBA,n2,n1)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   !local variables
   integer :: q,p,i,j,alpha,beta,gamma
   real(realk) :: tmpG1,tmpG2,tmpG3,tmpG4,tmpR1,tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO q=1,n2
      DO p=1,n2
         DO i=1,n1
            DO j=1,n1
               if(dopair_occ(i,j)) then
                  tmpR1 = 0.0E0_realk
                  DO alpha = 1, nBA
                     tmpR1 = tmpR1 + CalphaC(ALPHA,p,i)*CalphaC(ALPHA,q,j)
                  ENDDO
                  tmpG1 = 0.0E0_realk
                  DO beta = 1, nBA
                     tmpG1 = tmpG1 + CalphaC(BETA,p,i)*CalphaG(BETA,q,j)
                  ENDDO
                  tmpG2 = 0.0E0_realk
                  DO beta = 1, nBA
                     tmpG2 = tmpG2 + CalphaG(BETA,p,i)*CalphaC(BETA,q,j)
                  ENDDO
                  tmpG3 = 0.0E0_realk
                  tmpG4 = 0.0E0_realk
                  DO gamma = 1, nBA
                     tmpG3 = tmpG3 + CalphaG(GAMMA,p,j)*CalphaC(GAMMA,q,i)
                     tmpG4 = tmpG4 + CalphaC(GAMMA,p,j)*CalphaG(GAMMA,q,i)
                  ENDDO
                  EJ = EJ + tmpR1*(tmpG1 + tmpG2)
                  EK = EK + tmpR1*(tmpG3 + tmpG4)
               endif
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = 7.0/32.0*EJ + 1.0/32.0*EK 
   !print *,"COULOMBX2:  ", 7.0/32.0*EJ
   !print *,"EXCHANGEX2: ", 1.0/32.0*EK
   !print *,"COULOMBX2+EXCHANGEX2:", 7.0/32.0*EJ + 1.0/32.0*EK      

end subroutine ContractTwo4CenterF12IntegralsRIX_nc

subroutine ContractTwo4CenterF12IntegralsRIX_nc2(nBA,n1,n2,CalphaC,CalphaG,Fii,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaC(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaG(nBA,n2,n1)
   real(realk),intent(IN)    :: Fii(n1,n1)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   Real(realk)               :: Xijkl(n1,n1,n1,n1)
   !local variables
   integer :: q,p,i,j,k,l,alpha,beta,gamma
   real(realk) :: tmpG1,tmpG2,tmpG3,tmpG4,tmpR1,tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO i=1,n1
      DO j=1,n1
         DO k=1,n1
            DO l=1,n1

               tmp = 0.0E0_realk
               DO q=1,n2
                  DO p=1,n2
                     tmpR1 = 0.0E0_realk
                     DO alpha = 1, nBA
                        tmpR1 = tmpR1 + CalphaC(ALPHA,i,p)*CalphaC(ALPHA,j,q)
                     ENDDO
                     tmpG1 = 0.0E0_realk
                     DO beta = 1, nBA
                        tmpG1 = tmpG1 + CalphaC(BETA,k,p)*CalphaC(BETA,l,q)
                     ENDDO
                     tmp = tmp + tmpR1*tmpG1
                  ENDDO
               ENDDO
               Xijkl(i,j,k,l) = tmp

            ENDDO
         ENDDO
      ENDDO
   ENDDO

   DO l=1,n1
       DO k=1,n1
          DO j=1,n1
             DO i=1,n1
                if(abs(Xijkl(i,j,k,l)) > 1E-10) then
                   print *, "i j k l Xijkl: ", i, j, k, l, Xijkl(i,j,k,l)
                endif
             ENDDO
          ENDDO
       ENDDO
    ENDDO

   EJK = 7.0/32.0*EJ + 1.0/32.0*EK 
   print *,"COULOMBX2:  ", 7.0/32.0*EJ
   print *,"EXCHANGEX2: ", 1.0/32.0*EK
   print *,"COULOMBX2+EXCHANGEX2:", 7.0/32.0*EJ + 1.0/32.0*EK 

end subroutine ContractTwo4CenterF12IntegralsRIX_nc2

subroutine ContractTwo4CenterF12IntegralsRIX_nc3(nBA,n1,n2,CalphaC,Fii,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaC(nBA,n1,n2)
   real(realk),intent(in)    :: Fii(n1,n1)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   !local variables
   integer :: q,p,i,j,k,alpha,beta,gamma
   real(realk) :: tmpG1,tmpG2,tmpG3,tmpG4,tmpR1,tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO q=1,n2
      DO p=1,n2
         DO i=1,n1
            DO j=1,n1
               DO k=1,n1
                  tmpR1 = 0.0E0_realk
                  DO alpha = 1, nBA
                     tmpR1 = tmpR1 + CalphaC(ALPHA,i,p)*CalphaC(ALPHA,j,q)
                  ENDDO
                  tmpG1 = 0.0E0_realk
                  tmpG2 = 0.0E0_realk
                  DO beta = 1, nBA
                     !tmpG1 = tmpG1 + CalphaC(BETA,i,p)*CalphaC(BETA,k,q)*Fii(k,j)
                     !tmpG2 = tmpG2 + CalphaC(BETA,k,p)*CalphaC(BETA,j,q)*Fii(k,i)
                     tmpG1 = tmpG1 + CalphaC(BETA,i,p)*CalphaC(BETA,k,q)*Fii(k,j)
                     tmpG2 = tmpG2 + CalphaC(BETA,k,p)*CalphaC(BETA,j,q)*Fii(k,i)
                  ENDDO
                  tmpG3 = 0.0E0_realk
                  tmpG4 = 0.0E0_realk
                  DO gamma = 1, nBA
                     !tmpG3 = tmpG3 + CalphaC(GAMMA,k,p)*CalphaC(GAMMA,i,q)*Fii(k,j)
                     !tmpG4 = tmpG4 + CalphaC(GAMMA,j,p)*CalphaC(GAMMA,k,q)*Fii(k,i)
                     tmpG3 = tmpG3 + CalphaC(GAMMA,k,p)*CalphaC(GAMMA,i,q)*Fii(k,j)
                     tmpG4 = tmpG4 + CalphaC(GAMMA,j,p)*CalphaC(GAMMA,k,q)*Fii(k,i)
                  ENDDO
                  EJ = EJ + tmpR1*(tmpG1 + tmpG2)
                  EK = EK + tmpR1*(tmpG3 + tmpG4)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   EJK = 7.0/32.0*EJ + 1.0/32.0*EK 
   print *,"COULOMBX2:  ", 7.0/32.0*EJ
   print *,"EXCHANGEX2: ", 1.0/32.0*EK
   print *,"COULOMBX2+EXCHANGEX2:", 7.0/32.0*EJ + 1.0/32.0*EK      

end subroutine ContractTwo4CenterF12IntegralsRIX_nc3

subroutine ContractTwo4CenterF12IntegralsRI2V3V4(nBA,n1,n2,n3,nbasis,&
      & CalphaRcabs,CalphaGcabs,CalphaR,CalphaG,EJK3,EJK4,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3,nbasis
   real(realk),intent(in)    :: CalphaRcabs(nBA,n1,n3),CalphaGcabs(nBA,n1,n3)
   real(realk),intent(in)    :: CalphaR(nBA,n1,nbasis),CalphaG(nBA,nbasis,n1)
   real(realk),intent(inout) :: EJK3,EJK4
   real(realk)               :: EJ3, EJ4, EK3, EK4, ED
   !local variables
   integer :: m,c,i,j,alpha,beta
   real(realk) :: tmpR3,tmpG13,tmpG23
   real(realk) :: tmpR4,tmpG14,tmpG24
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   ED =  0.0E0_realk
   EJ3 =  0.0E0_realk
   EK3 =  0.0E0_realk
   EJK3 = 0.0E0_realk
   EJ4 =  0.0E0_realk
   EK4 =  0.0E0_realk
   EJK4 = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3,tmpR4, &
   !$OMP tmpG13,tmpG23,tmpG14,tmpG24) SHARED(CalphaRcabs,CalphaR,CalphaGcabs,CalphaG,&
   !$OMP n3,n2,n1,nba,dopair_occ) REDUCTION(+:EJ3,EK3,EJ4,EK4,ED)
   DO c=1,n3
      DO m=1,n2
         DO j=1,n1
            !Diagonal
            IF(dopair_occ(J,J)) THEN
               tmpR3 = 0.0E0_realk
               DO alpha = 1,NBA
                  tmpR3 = tmpR3 + CalphaR(alpha,j,m)*CalphaRcabs(alpha,j,c)
               ENDDO
               tmpG13 = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG13 = tmpG13 + CalphaG(beta,m,j)*CalphaGcabs(beta,j,c)
               ENDDO
               ED = ED + tmpR3*tmpG13
            ENDIF
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpR3 = 0.0E0_realk
                  tmpR4 = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR3 = tmpR3 + CalphaR(alpha,i,m)*CalphaRcabs(alpha,j,c)
                     tmpR4 = tmpR4 + CalphaR(alpha,j,m)*CalphaRcabs(alpha,i,c)
                  ENDDO              
                  tmpG13 = 0.0E0_realk
                  tmpG23 = 0.0E0_realk
                  !  tmpG14 = 0.0E0_realk
                  !  tmpG24 = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG13 = tmpG13 + CalphaG(beta,m,i)*CalphaGcabs(beta,j,c)
                     tmpG23 = tmpG23 + CalphaG(beta,m,j)*CalphaGcabs(beta,i,c)
                     ! tmpG14 = tmpG14 + CalphaGocc(beta,j,m)*CalphaG(beta,i,c)
                     ! tmpG24 = tmpG24 + CalphaGocc(beta,i,m)*CalphaG(beta,j,c)
                  ENDDO
                  EJ3 = EJ3 + tmpR3*tmpG13 
                  EK3 = EK3 + tmpR3*tmpG23
                  EJ4 = EJ4 + tmpR4*tmpG23
                  EK4 = EK4 + tmpR4*tmpG13
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK3 = ED + 2.5E0_realk*EJ3-0.5E0_realk*EK3
   EJK4 = ED + 2.5E0_realk*EJ4-0.5E0_realk*EK4
end subroutine ContractTwo4CenterF12IntegralsRI2V3V4

subroutine ContractTwo4CenterF12IntegralsRI2V3V4_p(nBA,n1,n3,n2,nbasis,&
      & CalphaRcabs,CalphaGcabs,CalphaR,CalphaG,EJK3,EJK4,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3,nbasis
   real(realk),intent(in)    :: CalphaRcabs(nBA,n1,n2),CalphaGcabs(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaR(nBA,n1,nbasis),CalphaG(nBA,n1,nbasis)
   real(realk),intent(inout) :: EJK3,EJK4
   real(realk)               :: EJ3, EJ4, EK3, EK4, ED
   !local variables
   integer :: m,c,i,j,alpha,beta
   real(realk) :: tmpR3,tmpG13,tmpG23
   real(realk) :: tmpR4,tmpG14,tmpG24
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   ED =  0.0E0_realk
   EJ3 =  0.0E0_realk
   EK3 =  0.0E0_realk
   EJK3 = 0.0E0_realk
   EJ4 =  0.0E0_realk
   EK4 =  0.0E0_realk
   EJK4 = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3,tmpR4, &
   !$OMP tmpG13,tmpG23,tmpG14,tmpG24) SHARED(CalphaRcabs,CalphaR,CalphaGcabs,CalphaG,&
   !$OMP n3,n2,n1,nba,dopair_occ) REDUCTION(+:EJ3,EK3,EJ4,EK4,ED)
   DO m=1,n3
      DO c=1,n2
         DO j=1,n1
            DO i=1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpR3 = 0.0E0_realk
                  tmpR4 = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR3 = tmpR3 + CalphaR(alpha,i,m)*CalphaRcabs(alpha,j,c)
                     tmpR4 = tmpR4 + CalphaR(alpha,j,m)*CalphaRcabs(alpha,i,c)
                  ENDDO              
                  tmpG13 = 0.0E0_realk
                  tmpG23 = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG13 = tmpG13 + CalphaG(beta,i,m)*CalphaGcabs(beta,j,c)
                     tmpG23 = tmpG23 + CalphaG(beta,j,m)*CalphaGcabs(beta,i,c)
                  ENDDO
                  EJ3 = EJ3 + tmpR3*tmpG13 
                  EK3 = EK3 + tmpR3*tmpG23
                  EJ4 = EJ4 + tmpR4*tmpG23
                  EK4 = EK4 + tmpR4*tmpG13
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK3 = 1.25E0_realk*EJ3-0.25E0_realk*EK3
   EJK4 = 1.25E0_realk*EJ4-0.25E0_realk*EK4

end subroutine ContractTwo4CenterF12IntegralsRI2V3V4_p

subroutine ContractTwo4CenterF12IntegralsRI2(nBA,n1,n3,n2,CalphaR,CalphaG,&
      & CalphaRocc,CalphaGocc,EJK3,EJK4,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaRocc(nBA,n1,n3),CalphaGocc(nBA,n1,n3)
   real(realk),intent(inout) :: EJK3,EJK4
   real(realk)               :: EJ3, EJ4, EK3, EK4, ED
   !local variables
   integer :: m,c,i,j,alpha,beta
   real(realk) :: tmpR3,tmpG13,tmpG23
   real(realk) :: tmpR4,tmpG14,tmpG24
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   ED =  0.0E0_realk
   EJ3 =  0.0E0_realk
   EK3 =  0.0E0_realk
   EJK3 = 0.0E0_realk
   EJ4 =  0.0E0_realk
   EK4 =  0.0E0_realk
   EJK4 = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3,tmpR4, &
   !$OMP tmpG13,tmpG23,tmpG14,tmpG24) SHARED(CalphaR,CalphaRocc,CalphaG,CalphaGocc,n3,n2,n1,&
   !$OMP nba,dopair_occ) REDUCTION(+:EJ3,EK3,EJ4,EK4,ED)
   DO m=1,n3
      DO c=1,n2
         DO j=1,n1
            !Diagonal
            IF(dopair_occ(J,J)) THEN
               tmpR3 = 0.0E0_realk
               DO alpha = 1,NBA
                  tmpR3 = tmpR3 + CalphaRocc(alpha,j,m)*CalphaR(alpha,j,c)
               ENDDO
               tmpG13 = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG13 = tmpG13 + CalphaGocc(beta,j,m)*CalphaG(beta,j,c)
               ENDDO
               ED = ED + tmpR3*tmpG13
            ENDIF
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpR3 = 0.0E0_realk
                  tmpR4 = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR3 = tmpR3 + CalphaRocc(alpha,i,m)*CalphaR(alpha,j,c)
                     tmpR4 = tmpR4 + CalphaRocc(alpha,j,m)*CalphaR(alpha,i,c)
                  ENDDO              
                  tmpG13 = 0.0E0_realk
                  tmpG23 = 0.0E0_realk
                  tmpG14 = 0.0E0_realk
                  tmpG24 = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG13 = tmpG13 + CalphaGocc(beta,i,m)*CalphaG(beta,j,c)
                     tmpG23 = tmpG23 + CalphaGocc(beta,j,m)*CalphaG(beta,i,c)

                     ! tmpG14 = tmpG14 + CalphaGocc(beta,j,m)*CalphaG(beta,i,c)
                     ! tmpG24 = tmpG24 + CalphaGocc(beta,i,m)*CalphaG(beta,j,c)
                  ENDDO
                  EJ3 = EJ3 + tmpR3*tmpG13 
                  EK3 = EK3 + tmpR3*tmpG23
                  EJ4 = EJ4 + tmpR4*tmpG23
                  EK4 = EK4 + tmpR4*tmpG13
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK3 = ED + 2.5E0_realk*EJ3-0.5E0_realk*EK3
   EJK4 = ED + 2.5E0_realk*EJ4-0.5E0_realk*EK4

end subroutine ContractTwo4CenterF12IntegralsRI2

subroutine ContractTwo4CenterF12IntegralsRIB4(nBA,n1,n2,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaD(nBa,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,q,r,i,j,alpha,beta,alpha1,beta1,alpha2,beta2,alpha3,beta3,alpha4,beta4
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   ED = 0.0E0_realk
   EJ = 0.0E0_realk
   EK = 0.0E0_realk
   DO q=1,n2 !ncabsAO
      DO r=1,n2 !ncabsAO
         DO j=1,n1 !nocc
            !Diagonal
            IF(dopair_occ(J,J)) THEN
               tmpR = 0.0E0_realk
               tmpG = 0.0E0_realk
               DO beta = 1,nBA
                  tmpR = tmpR + CalphaG(beta,j,r)*CalphaD(beta,j,q)
                  tmpG = tmpG + CalphaG(beta,j,r)*CalphaG(beta,j,q)
               ENDDO
               ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            ENDIF 
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  tmpGJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaG(alpha1,i,r)*CalphaD(alpha1,j,q) 
                     tmpGJ1 = tmpGJ1 + CalphaG(alpha1,i,r)*CalphaG(alpha1,j,q)
                  ENDDO
                  tmpRJ2 = 0.0E0_realk
                  tmpGJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaG(alpha2,j,r)*CalphaD(alpha2,i,q)
                     tmpGJ2 = tmpGJ2 + CalphaG(alpha2,j,r)*CalphaG(alpha2,i,q)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO   
   EJK = ED*0.5E0_realk + 7.0/16.0*EJ + 1.0/16.0*EK 
end subroutine ContractTwo4CenterF12IntegralsRIB4                

subroutine ContractTwo4CenterF12IntegralsRIB5(nBA,n1,n2,n3,nbasis,CalphaGcabs,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3,nbasis
   real(realk),intent(in)    :: CalphaG(nBA,nbasis,n1)
   real(realk),intent(in)    :: CalphaGcabs(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,q,m,i,j,alpha,beta,alpha1,beta1,alpha2,beta2
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   ED = 0.0E0_realk
   EJ = 0.0E0_realk
   EK = 0.0E0_realk
   DO q=1,n2 !ncabsAO
      DO m=1,n3 !noccAOS
         DO j=1,n1 !nocc
            !Diagonal
            IF(dopair_occ(J,J)) THEN
               tmpR = 0.0E0_realk
               tmpG = 0.0E0_realk
               DO alpha = 1,nBA
                  tmpR = tmpR + CalphaD(alpha,j,q)*CalphaG(alpha,m,j)
               ENDDO
               DO beta = 1,nBA
                  tmpG = tmpG + CalphaGcabs(beta,j,q)*CalphaG(beta,m,j)
               ENDDO
               ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            ENDIF 
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  tmpGJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaD(alpha1,i,q)*CalphaG(alpha1,m,j) 
                     tmpGJ1 = tmpGJ1 + CalphaGcabs(alpha1,i,q)*CalphaG(alpha1,m,j)
                  ENDDO
                  tmpRJ2 = 0.0E0_realk
                  tmpGJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaD(alpha2,j,q)*CalphaG(alpha2,m,i)
                     tmpGJ2 = tmpGJ2 + CalphaGcabs(alpha2,j,q)*CalphaG(alpha2,m,i)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = -1.0E0_realk*(ED*0.5E0_realk + 7.0/16.0*EJ + 1.0/16.0*EK) 
end subroutine ContractTwo4CenterF12IntegralsRIB5                

subroutine ContractTwo4CenterF12IntegralsRIB6(nBA,n1,n2,n3,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(in)    :: CalphaG(nBA,n3,n1)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: q,a,i,j,alpha,beta,alpha1,beta1,alpha2,beta2
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
   ED = 0.0E0_realk
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   DO q=1,n3 !nocv
      DO a=n1+1,n3 !nvirt
         DO j=1,n1 !nocc
            !Diagonal
            IF(dopair_occ(J,J)) THEN
               tmpR = 0.0E0_realk
               tmpG = 0.0E0_realk
               DO alpha = 1,nBA
                  tmpR = tmpR + CalphaD(alpha,j,q)*CalphaG(alpha,a,j)
                  tmpG = tmpG + CalphaG(alpha,q,j)*CalphaG(alpha,a,j)
               ENDDO
               ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            ENDIF
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  tmpGJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaD(alpha1,i,q)*CalphaG(alpha1,a,j) 
                     tmpGJ1 = tmpGJ1 + CalphaG(alpha1,q,i)*CalphaG(alpha1,a,j)
                  ENDDO
                  tmpRJ2 = 0.0E0_realk
                  tmpGJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaD(alpha2,j,q)*CalphaG(alpha2,a,i)
                     tmpGJ2 = tmpGJ2 + CalphaG(alpha2,q,j)*CalphaG(alpha2,a,i)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)           
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO   
   EJK = -1.0E0_realk*(ED*0.5E0_realk + 7.0/16.0*EJ + 1.0/16.0*EK) 
end subroutine ContractTwo4CenterF12IntegralsRIB6
            
!warning I am very unclear about the Frozen core implementation of this 
subroutine ContractOccCalpha(NBA,nocc,noccfull,nbasis,CalphaG,Fii,CalphaD,dopair_occ_in)
   implicit none
   integer,intent(in) :: NBA,nocc,noccfull,nbasis
   real(realk),intent(in) :: CalphaG(NBA,nbasis,nocc),Fii(noccfull,noccfull)
   real(realk),intent(inout) :: CalphaD(NBA,nocc,noccfull)
   integer :: i,j,alpha,m
   real(realk) :: TMP
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(nocc,nocc)
   logical :: dopair_occ(nocc,nocc)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
!   !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,alpha,m,&
!   !$OMP TMP) SHARED(NBA,nocc,noccfull,CalphaG,CalphaD,Fii)
!   !$OMP DO COLLAPSE(3)
   DO m=1,noccfull
      DO I=1,nocc
         IF(dopair_occ(I,I)) THEN
            DO ALPHA=1,NBA
               CalphaD(ALPHA,I,m) = CalphaG(ALPHA,I,1)*Fii(1,M)
            ENDDO
         ENDIF
      ENDDO
   ENDDO
!   !$OMP END DO
!   !$OMP DO COLLAPSE(2) 
   DO m=1,noccfull
      DO I=1,nocc
         DO J=2,noccfull
            IF(dopair_occ(I,J)) THEN
               TMP = Fii(J,M)
               DO ALPHA=1,NBA
                  CalphaD(ALPHA,I,m) = CalphaD(ALPHA,I,m) + CalphaG(ALPHA,I,J)*TMP
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDDO
!   !$OMP END DO
!   !$OMP END PARALLEL
end subroutine ContractOccCalpha


subroutine ContractTwo4CenterF12IntegralsRIB7(nBA,n1,n2,n3,CalphaR,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaR(nBA,n1,n3)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: c,n,i,j,alpha,beta,alpha1,beta1,alpha2,beta2
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif  
   ED = 0.0E0_realk
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   DO c=1,n3 !ncabsMO
      DO n=1,n2 !noccAOS
         DO j=1,n1 !noccEOS
            !Diagonal
            IF(dopair_occ(J,J)) THEN
               tmpR = 0.0E0_realk
               tmpG = 0.0E0_realk
               DO alpha = 1,nBA
                  tmpR = tmpR + CalphaR(alpha,j,c)*CalphaD(alpha,j,n)
               ENDDO
               DO beta = 1,nBA
                  tmpG = tmpG + CalphaR(beta,j,c)*CalphaG(beta,j,n)
               ENDDO
               ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            ENDIF
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaR(alpha1,i,c)*CalphaD(alpha1,j,n) 
                  ENDDO
                  tmpRJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaR(alpha2,j,c)*CalphaD(alpha2,i,n)
                  ENDDO
                  tmpGJ1 = 0.0E0_realk
                  DO beta1 = 1, nBA
                     tmpGJ1 = tmpGJ1 + CalphaR(beta1,i,c)*CalphaG(beta1,j,n)
                  ENDDO
                  tmpGJ2 = 0.0E0_realk
                  DO beta2 = 1, nBA
                     tmpGJ2 = tmpGJ2 + CalphaR(beta2,j,c)*CalphaG(beta2,i,n)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO   
   EJK = 1.0E0_realk*(ED*0.5E0_realk + 7.0/16.0*EJ + 1.0/16.0*EK) 

end subroutine ContractTwo4CenterF12IntegralsRIB7

subroutine ContractTwo4CenterF12IntegralsRIB8(nBA,n1,n2,nbasis,CalphaR,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,nbasis
   real(realk),intent(in)    :: CalphaG(nBA,nbasis,n1)
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2) !CalphaGcabsMO
   real(realk),intent(in)    :: CalphaD(nBA,n1,n1)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,m,i,j,alpha,beta,alpha1,beta1,alpha2,beta2
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   ED = 0.0E0_realk
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(p,m,i,j,alpha,beta,alpha1,&
   !$OMP beta1,alpha2,beta2,tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2,tmpG,tmpGJ1,tmpGJ2,&
   !$OMP tmpGK1,tmpGK2) SHARED(CalphaR,CalphaG,CalphaD,n2,n1,nba,dopair_occ) REDUCTION(+:ED,EJ,EK)
   DO p=1,n2 !ncabs
      DO m=1,n1 !noccfull
         DO j=1,n1 !nocc
            IF(dopair_occ(J,J)) THEN
               !Diagonal
               tmpR = 0.0E0_realk
               tmpG = 0.0E0_realk
               DO alpha = 1,nBA
                  tmpR = tmpR + CalphaR(alpha,j,p)*CalphaG(alpha,m,j)
               ENDDO
               DO beta = 1,nBA
                  tmpG = tmpG + CalphaR(beta,j,p)*CalphaD(beta,j,m)
               ENDDO
               ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            ENDIF
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaR(alpha1,i,p)*CalphaG(alpha1,m,j) 
                  ENDDO
                  tmpRJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaR(alpha2,j,p)*CalphaG(alpha2,m,i)
                  ENDDO
                  tmpGJ1 = 0.0E0_realk
                  DO beta1 = 1, nBA
                     tmpGJ1 = tmpGJ1 + CalphaR(beta1,i,p)*CalphaD(beta1,j,m)
                  ENDDO
                  tmpGJ2 = 0.0E0_realk
                  DO beta2 = 1, nBA
                     tmpGJ2 = tmpGJ2 + CalphaR(beta2,j,p)*CalphaD(beta2,i,m)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK = -2.0E0_realk*(ED*0.5E0_realk + 7.0_realk/16.0_realk*EJ + 1.0_realk/16.0_realk*EK)
end subroutine ContractTwo4CenterF12IntegralsRIB8

subroutine ContractTwo4CenterF12IntegralsRIB9(nBA,n1,n2,n3,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(in)    :: CalphaG(nBA,n3,n1)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,a,i,j,alpha,beta,alpha1,beta1,alpha2,beta2
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
   ED = 0.0E0_realk
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   DO p=1,n3 !ncabs
      DO a=n1+1,n3 !nvirt
         DO j=1,n1 !nocc
            IF(dopair_occ(J,J)) THEN
               !Diagonal
               tmpR = 0.0E0_realk
               DO alpha = 1,nBA
                  tmpR = tmpR + CalphaG(alpha,p,j)*CalphaG(alpha,a,j)
               ENDDO
               tmpG = 0.0E0_realk
               DO beta = 1,nBA
                  tmpG = tmpG + CalphaD(beta,j,p)*CalphaG(beta,a,j)
               ENDDO
               ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            ENDIF
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaG(alpha1,p,i)*CalphaG(alpha1,a,j) 
                  ENDDO      
                  tmpRJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaG(alpha2,p,j)*CalphaG(alpha2,a,i)
                  ENDDO
                  tmpGJ1 = 0.0E0_realk
                  DO beta1 = 1, nBA
                     tmpGJ1 = tmpGJ1 + CalphaD(beta1,i,p)*CalphaG(beta1,a,j)
                  ENDDO
                  tmpGJ2 = 0.0E0_realk
                  DO beta2 = 1, nBA
                     tmpGJ2 = tmpGJ2 + CalphaD(beta2,j,p)*CalphaG(beta2,a,i)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = -2.0E0_realk*(ED*0.5E0_realk + 7.0_realk/16.0_realk*EJ + 1.0_realk/16.0_realk*EK)
end subroutine ContractTwo4CenterF12IntegralsRIB9

subroutine ContractTwo4CenterF12IntegralsRIX3X4(nBA,n1,n3,n2,&
      & CalphaGcabs,CalphaG,Fii,EJK3,EJK4,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(in)    :: CalphaGcabs(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaG(nBA,n1,n3)
   real(realk),intent(IN)    :: Fii(n1,n1)
   real(realk),intent(inout) :: EJK3,EJK4
   real(realk)               :: EJ3, EJ4, EK3, EK4, ED
   !local variables
   integer :: m,c,i,j,alpha,beta
   real(realk) :: tmpR3,tmpG13,tmpG23
   real(realk) :: tmpR4,tmpG14,tmpG24,tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   ED =  0.0E0_realk
   EJ3 =  0.0E0_realk
   EK3 =  0.0E0_realk
   EJK3 = 0.0E0_realk
   EJ4 =  0.0E0_realk
   EK4 =  0.0E0_realk
   EJK4 = 0.0E0_realk
   !!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3,tmpR4, &
   !!$OMP tmpG13,tmpG23,tmp) SHARED(CalphaG,CalphaGcabs,n3,n2,n1,dopair_occ) &
   !!$OMP nba, Fii) REDUCTION(+:EJ3,EK3,EJ4,EK4,ED)
   DO m=1,n3
      DO c=1,n2
         DO j=1,n1
            if(dopair_occ(J,J)) then
               !Diagonal
               tmpG13 = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG13 = tmpG13 + CalphaG(beta,j,m)*CalphaGcabs(beta,j,c)
               ENDDO
               ED = ED + tmpG13*tmpG13*Fii(j,j)
            endif
            !Non Diagonal
            DO i=j+1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpG13 = 0.0E0_realk
                  tmpG23 = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG13 = tmpG13 + CalphaG(beta,i,m)*CalphaGcabs(beta,j,c)
                     tmpG23 = tmpG23 + CalphaG(beta,j,m)*CalphaGcabs(beta,i,c)
                     ! tmpG14 = tmpG14 + CalphaGocc(beta,j,m)*CalphaG(beta,i,c)
                     ! tmpG24 = tmpG24 + CalphaGocc(beta,i,m)*CalphaG(beta,j,c)
                  ENDDO
                  tmp = (Fii(i,i) + Fii(j,j)) 
                  EJ3 = EJ3 + tmpG13*tmpG13*tmp
                  EK3 = EK3 + tmpG13*tmpG23*tmp
                  EJ4 = EJ4 + tmpG23*tmpG23*tmp
                  EK4 = EK4 + tmpG23*tmpG13*tmp
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !!$OMP END PARALLEL DO

   EJK3 = ED*0.5E0_realk + (7.0_realk/16.0_realk*EJ3+1.0_realk/16.0_realk*EK3)
   EJK4 = ED*0.5E0_realk + (7.0_realk/16.0_realk*EJ4+1.0_realk/16.0_realk*EK4)
end subroutine ContractTwo4CenterF12IntegralsRIX3X4

subroutine ContractTwo4CenterF12IntegralsRIX3X4_nc(nBA,n1,n2,n3,&
      & CalphaGcabs,CalphaC,CalphaG,CalphaP,EJK3,EJK4,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(IN)    :: CalphaGcabs(nBA,n1,n3)
   real(realk),intent(IN)    :: CalphaC(nBA,n1,n2)
   real(realk),intent(IN)    :: CalphaG(nBA,n2,n1)
   real(realk),intent(IN)    :: CalphaP(nBA,n3,n1)
   real(realk),intent(inout) :: EJK3,EJK4
   real(realk)               :: EJ3, EJ4, EK3, EK4, ED
   !local variables
   integer :: m,c,i,j,alpha,beta
   real(realk) :: tmpR3,tmpG31,tmpG32,tmpG33,tmpG34
   real(realk) :: tmpR4,tmpG41,tmpG42,tmpG43,tmpG44
   real(realk) :: tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   ED =  0.0E0_realk
   EJ3 =  0.0E0_realk
   EK3 =  0.0E0_realk
   EJK3 = 0.0E0_realk
   EJ4 =  0.0E0_realk
   EK4 =  0.0E0_realk
   EJK4 = 0.0E0_realk
   !!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3,tmpR4, &
   !!$OMP tmpG13,tmpG23,tmp) SHARED(CalphaC,CalphaP,CalphaG,CalphaGcabs,n3,n2,n1,dopair_occ) &
   !!$OMP nba, Fii) REDUCTION(+:EJ3,EK3,EJ4,EK4,ED)
   DO c=1,n3
      DO m=1,n2
         DO j=1,n1
            DO i=1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpR3 = 0.0E0_realk
                  tmpR4 = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR3 = tmpR3 + CalphaC(alpha,i,m)*CalphaGcabs(alpha,j,c)
                  ENDDO
                  tmpG31 = 0.0E0_realk
                  tmpG32 = 0.0E0_realk
                  tmpG33 = 0.0E0_realk
                  tmpG34 = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG31 = tmpG31 + CalphaC(beta,m,i)*CalphaP(beta,c,j)
                     tmpG32 = tmpG32 + CalphaG(beta,m,i)*CalphaGcabs(beta,j,c)
                     tmpG33 = tmpG33 + CalphaG(beta,m,j)*CalphaGcabs(beta,i,c)
                     tmpG34 = tmpG34 + CalphaC(beta,m,j)*CalphaP(beta,c,i)
                  ENDDO
                  EJ3 = EJ3 + tmpR3*(tmpG31 + tmpG32)
                  EK3 = EK3 + tmpR3*(tmpG33 + tmpG34)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !!$OMP END PARALLEL DO
   EJK3 = 7.0/32.0_realk*EJ3+1.0_realk/32.0*EK3
   EJK4 = EJK3
!   print *,"COULOMBX2:  ", 7.0/32.0*EJ3
!   print *,"EXCHANGEX2: ", 1.0/32.0*EK3
!   print *,"COULOMBX2+EXCHANGEX2:", 7.0/32.0*EJ3 + 1.0/32.0*EK3      
end subroutine ContractTwo4CenterF12IntegralsRIX3X4_nc

end module f12_routines_module
