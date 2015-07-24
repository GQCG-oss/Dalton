!> @file
!> Author: Yang M. Wang
!> Date: April 2013
!> Calculates the expressions for the single fragment energy and pair fragment energies for MP2F12

!> General index convention:
!> *************************
!> a,b   : Virtual EOS
!> c,d,e : Virtual AOS
!> i,j   : Occupied EOS
!> k,l,m : Occupied AOS

module f12_integrals_module 

#ifdef MOD_UNRELEASED 

#ifdef VAR_MPI   
  use infpar_module
  use lsmpi_type
#endif
  use,intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
  use precision
  use lstiming!, only: lstimer
  use screen_mod!, only: DECscreenITEM
  use dec_typedef_module
  use typedeftype!, only: Lsitem,lssetting
  use BUILDAOBATCH!,only:build_batchesofaos,determine_maxbatchorbitalsize,&
  !      & determine_MaxOrbitals
  use typedef!, only: typedef_free_setting,copy_setting
  use memory_handling
  use screen_mod!,only: free_decscreen, DECSCREENITEM
  use lsparameters
  use IntegralInterfaceMod!, only: II_getBatchOrbitalInfo
  use IntegralInterfaceDEC!, only: II_precalc_DECScreenMat,&
  !       & II_getBatchOrbitalScreen

  use ccintegrals!, only: get_full_AO_integrals,get_AO_hJ,get_AO_K,get_AO_Fock

  ! Yangs F12 routines
  use f12_routines_module!, only: MO_transform_AOMatrix, matrix_print, norm4D, norm2D, get_mp2f12_MO

  ! Patricks mat_transpose routine 
  use reorder_frontend_module!, only: mat_transpose(rows,column,pref1,A,pref2,AT)

  ! Thomas free_cabs() for aa free MO_CABS_save_created, CMO_RI_save_created
  use CABS_operations

  ! MP2F12 C coupling routine
  use mp2_module

  ! *********************************************
  !   DEC DEPENDENCIES (within deccc directory) 
  ! *********************************************
#ifdef VAR_MPI
  use decmpi_module !, only: mpi_communicate_mp2_int_and_amp
#endif

  use dec_fragment_utils!,only: calculate_fragment_memory, &
  !       & dec_simple_dgemm_update,start_flop_counter,&
  !       & get_currently_available_memory, atomic_fragment_free
  use array2_simple_operations!, only: array2_free, array2_extract_EOS, &
  !       & get_mp2_integral_transformation_matrices, get_mp2_integral_transformation_matrices_fc, &
  !      & extract_occupied_eos_mo_indices, extract_virtual_EOS_MO_indices,array2_init,array2_print
  use array4_simple_operations!, only: array4_delete_file, array4_init_file, &
  !       & array4_init_standard, array4_free, array4_reorder, array4_init, &
  !       & array4_contract1, array4_open_file, array4_write_file_type2, &
  !       & array4_close_file, array4_write_file_type1, mat_transpose, &
  !     & array4_read_file_type2

  public :: get_f12_fragment_energy, matrix_print_4d, matrix_print_2d, get_mp2f12_sf_E21, get_f12_fragment_energy_slave

  private
#endif

contains
#ifdef MOD_UNRELEASED 

  !> Brief: Gives the single and pair fragment energy for V1 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EV1(Venergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V1energy

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 
    real(realk), pointer :: V1ijkl(:,:,:,:)          !V1_ijkl <ij|f12*r^-1|kl>

    integer :: noccEOS !number of occupied MO orbitals in EOS
    noccEOS  = MyFragment%noccEOS

    call mem_alloc(V1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)

    !> Get integrals <ij|f12*r^-1|kl> stored as (i,j,k,l)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,&
         & Ccabs,Cri,CvirtAOS,'iiii','RRRRF',V1ijkl)

    V1energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E21(V1ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V1energy, 1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V1ijkl, noccEOS, V1energy,  1.0E0_realk)
    endif

    Venergy(1) = V1energy

    call mem_dealloc(V1ijkl)
  end subroutine get_EV1

  !> Brief: Gives the single and pair fragment energy for V2 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EV2(Venergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V2energy

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccAOS(nbasis,nocctot)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), pointer :: V2ijkl(:,:,:,:)! sum_pq <ij|r^-1|pq> * <pq|f12|kl>   
    real(realk), pointer :: Gijpq(:,:,:,:) ! <ij|r^-1|pq>
    real(realk), pointer :: Rijpq(:,:,:,:) ! <ij|f12|pq>

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: nocvAOStot !number of occupied tot + virtual MO orbitals in EOS 

    integer :: m, k, n

    noccEOS  = MyFragment%noccEOS
    nocvAOStot  = MyFragment%nocctot + MyFragment%nvirtAOS

    call mem_alloc(V2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)  
    call mem_alloc(Gijpq,  noccEOS, noccEOS, nocvAOStot, nocvAOStot)    
    call mem_alloc(Rijpq,  noccEOS, noccEOS, nocvAOStot, nocvAOStot)

    V2energy = 0.0E0_realk

    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iipp','RRRRC',Gijpq)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iipp','RRRRG',Rijpq)

    m = noccEOS*noccEOS  ! <ij G pq> <pq R kl> = <m V2 n> 
    k = nocvAOStot*nocvAOStot  
    n = noccEOS*noccEOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','T',m,n,k,1.0E0_realk,Rijpq,m,Gijpq,n,0.0E0_realk,V2ijkl,m)

    V2energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E21(V2ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V2energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V2ijkl, noccEOS, V2energy, -1.0E0_realk)
    endif

    Venergy(2) = V2energy

    call mem_dealloc(V2ijkl)
    call mem_dealloc(Gijpq) 
    call mem_dealloc(Rijpq)

  end subroutine get_EV2

  !> Brief: Gives the single and pair fragment energy for V3 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EV3(Venergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V3energy
    real(realk) :: V4energy

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), pointer :: V3ijkl(:,:,:,:) ! sum_mc <ij|r^-1|mc> * <mc|f12|kl>   
    real(realk), pointer :: Gijmc(:,:,:,:)  ! <ij|r^-1|ma'>
    real(realk), pointer :: Rijmc(:,:,:,:)  ! <ij|f12|ma'>
    real(realk), pointer :: V4ijkl(:,:,:,:) ! sum_mc <ji|r^-1|cm> * <cm|f12|lk> 

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOStot !number of occupied orbitals in AOStot 
    integer :: ncabsMO !number of CABS MO orbitals    

    integer :: m, k, n

    noccEOS  = MyFragment%noccEOS
    noccAOStot  = MyFragment%nocctot 
    ncabsMO = size(MyFragment%Ccabs,2)

    call mem_alloc(V3ijkl, noccEOS, noccEOS, noccEOS, noccEOS)  
    call mem_alloc(Gijmc,  noccEOS, noccEOS, noccAOStot, ncabsMO)
    call mem_alloc(Rijmc,  noccEOS, noccEOS, noccAOStot, ncabsMO)

    call mem_alloc(V4ijkl, noccEOS, noccEOS, noccEOS, noccEOS) 

    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS, &
       & CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iimc','RCRRC',Gijmc)    
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS, &
       & CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iimc','RCRRG',Rijmc) 

    m = noccEOS*noccEOS  ! <ij G mc> <mc R kl> = <m V3 n> 
    k = noccAOStot*ncabsMO  ! m x k * k x n = m x n
    n = noccEOS*noccEOS

    !> dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','T',m,n,k,1.0E0_realk,Gijmc,m,Rijmc,n,0.0E0_realk,V3ijkl,m)

    !> Creating the V4ijkl = V3jilk 
    call array_reorder_4d(1.0E0_realk,V3ijkl,noccEOS,noccEOS,noccEOS,noccEOS,[2,1,4,3],0.0E0_realk,V4ijkl)

    V3energy = 0.0E0_realk
    V4energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E21(V3ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V3energy, -1.0E0_realk)
       call get_mp2f12_pf_E21(V4ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V4energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V3ijkl, noccEOS, V3energy, -1.0E0_realk)
       call get_mp2f12_sf_E21(V4ijkl, noccEOS, V4energy, -1.0E0_realk)
    endif

    Venergy(3) = V3energy
    Venergy(4) = V4energy

    call mem_dealloc(V3ijkl)
    call mem_dealloc(Gijmc)
    call mem_dealloc(Rijmc)

    call mem_dealloc(V4ijkl)

  end subroutine get_EV3

  !> Brief: Gives the single and pair fragment energy for V4 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EV4(Venergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Taibj)
    implicit none

    real(realk), intent(inout) :: Venergy(:) 
    real(realk) :: V5energy

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOStot)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Taibj(:,:,:,:)  

    real(realk), pointer :: V5ijkl(:,:,:,:) !sum_ab <ij|r^-1|ab> * <ab|t_2|kl>  
    real(realk), pointer :: Cijab(:,:,:,:)
    real(realk), pointer :: Rijac(:,:,:,:)

    real(realk), pointer :: Tabij(:,:,:,:) 

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOStot
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: i, j, m, k, n, a, b, c
    real(realk) :: tmp

    noccEOS  = MyFragment%noccEOS
    noccAOStot  = MyFragment%nocctot
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)

    call mem_alloc(V5ijkl, noccEOS, noccEOS, noccEOS, noccEOS) 
    call mem_alloc(Cijab, noccEOS,  noccEOS,  nvirtAOS,  nvirtAOS) 
    call mem_alloc(Rijac, noccEOS,  noccEOS,  nvirtAOS,  ncabsMO)
    call mem_alloc(Tabij, nvirtAOS, nvirtAOS,  noccEOS,  noccEOS)

    V5energy = 0.0E0_realk
    tmp = 0.0E0_realk

    ! **********************************************************
    !   Creating the C matrix 
    ! **********************************************************  
    !> Rijac <ij|f12|ac> stored as (i,j,a,c)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iiac','RCRRG',Rijac)
!FIXME ADD A C Couplings contribution here 
    Cijab = 0.0E0_realk
    do i=1, noccEOS
       do j=1, noccEOS      
          do a=1, nvirtAOS
             do b=1, nvirtAOS      
                ! tmp2 = 0.0E0_realk
                tmp  = 0.0E0_realk
                do c=1, ncabsMO
                   tmp =  tmp  + Rijac(i,j,a,c)*(Myfragment%Fcp(c,b+noccAOStot)) + Rijac(j,i,b,c)*(Myfragment%Fcp(c,a+noccAOStot)) 
                   ! tmp2 = tmp2 + Rijac(j,i,a,c)*(Myfragment%Fcp(c,b+noccAOStot)) + Rijac(i,j,b,c)*(Myfragment%Fcp(c,a+noccAOStot)) 
                enddo
                Cijab(i,j,a,b) = tmp 
                   ! Cijab(j,i,a,b) = tmp2 
             end do
          end do
       end do
    end do

    m = noccEOS*noccEOS    ! <ij C ab> <ab T kl> = <m V2 n> 
    k = nvirtAOS*nvirtAOS
    n = noccEOS*noccEOS

    call array_reorder_4d(1.0E0_realk,Taibj,nvirtAOS,noccEOS,nvirtAOS,noccEOS,[1,3,2,4],0.0E0_realk,Tabij)

    !> dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Cijab,m,Tabij,k,0.0E0_realk,V5ijkl,m)

    if(dopair) then
       call get_mp2f12_pf_E21(V5ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V5energy, 1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V5ijkl, noccEOS, V5energy,  1.0E0_realk)
    endif

    Venergy(5) = V5energy

    call mem_dealloc(Cijab)
    call mem_dealloc(Rijac)
    call mem_dealloc(Tabij)

    call mem_dealloc(V5ijkl)

  end subroutine get_EV4

  !> Brief: Gives the single and pair fragment energy for X1 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EX1(Xenergy,Fkj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none

    real(realk), intent(inout) :: Xenergy(:)
    real(realk) :: X1energy, tmp1, tmp2
    ! k index: occupied AOS, for frozen core: only valence!
    real(realk), target, intent(in) :: Fkj(:,:)
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOStot)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 
    real(realk), pointer :: X1ijkl(:,:,:,:)
    
    real(realk), pointer :: X1ijkn(:,:,:,:)
    real(realk), pointer :: X1ijnk(:,:,:,:)
    
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: i,j,k,n

    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS  ! (only valence for frozen core)

    call mem_alloc(X1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)

    call mem_alloc(X1ijkn, noccEOS, noccEOS, noccEOS, noccAOS)
    call mem_alloc(X1ijnk, noccEOS, noccEOS, noccAOS, noccEOS)

    ! For frozen core: index 4 (v) is only virtual
    ! Without frozen core: index 4 is core+virtual (effectively v=m)

    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iiiv','RRRR2',X1ijkn)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iivi','RRRR2',X1ijnk)

    !> Creating the X2ijkl
    do i=1,noccEOS
       do j=1,noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do k=1,noccAOS
             tmp1 = tmp1 - (X1ijkn(i,j,i,k)*Fkj(k,j) + X1ijnk(i,j,k,j)*Fkj(k,i))  
             tmp2 = tmp2 - (X1ijkn(j,i,i,k)*Fkj(k,j) + X1ijkn(i,j,j,k)*Fkj(k,i)) 
          enddo
          X1ijkl(i,j,j,i) =  tmp2
          X1ijkl(i,j,i,j) =  tmp1     !> Important ordering when i=j
       enddo 
    enddo

    X1energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E22(X1ijkl, noccEOS, Fragment1, Fragment2, MyFragment, X1energy,  1.0E0_realk)
    else    
       call get_mp2f12_sf_E22(X1ijkl, noccEOS, X1energy,  1.0E0_realk)
    endif

    Xenergy(1) = X1energy

    call mem_dealloc(X1ijkl)
    call mem_dealloc(X1ijkn)
    call mem_dealloc(X1ijnk)

  end subroutine get_EX1

  !> Brief: Gives the single and pair fragment energy for X2 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EX2(Xenergy,Fkj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none

    real(realk), intent(inout) :: Xenergy(:)
    real(realk) :: X2energy, tmp1, tmp2 

    real(realk), target, intent(in) :: Fkj(:,:)
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), pointer :: X2ijkl(:,:,:,:) !sum_pq X2ijkn * Fnl 
    real(realk), pointer :: X2ijkn(:,:,:,:) !sum_pq <ij|f12|pq> * <pq|f12|kn> 
    real(realk), pointer :: X2ijnk(:,:,:,:) 

    real(realk), pointer :: Rijpq(:,:,:,:)  !<ij|f12|pq>
    real(realk), pointer :: Rpqkn(:,:,:,:)  !<pq|f12|kn>  
    real(realk), pointer :: Rpqnk(:,:,:,:)  !<pq|f12|nk>  

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nocvAOStot
    integer :: m,k,n,i,j

    noccEOS = MyFragment%noccEOS
    noccAOS = MyFragment%noccAOS ! only valence for frozen core 
    nocvAOStot = MyFragment%nocctot + MyFragment%nvirtAOS

    call mem_alloc(X2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(X2ijkn, noccEOS, noccEOS, noccEOS, noccAOS)
    call mem_alloc(X2ijnk, noccEOS, noccEOS, noccAOS, noccEOS)

    call mem_alloc(Rijpq,  noccEOS, noccEOS, nocvAOStot, nocvAOStot)
    call mem_alloc(Rpqkn,  nocvAOStot, nocvAOStot, noccEOS, noccAOS)
    call mem_alloc(Rpqnk,  nocvAOStot, nocvAOStot, noccAOS, noccEOS)

    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iipp','RRRRG',Rijpq)

    ! For frozen core: index 4 (v) is only virtual
    ! Without frozen core: index 4 is core+virtual (effectively v=m) 
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'ppiv','RRRRG',Rpqkn)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'ppvi','RRRRG',Rpqnk)

    ! These m,k,n refer to dgemm inputs
    m = noccEOS*noccEOS   ! <ij R pq> <pq R kn> = <m X2 n>
    k = nocvAOStot*nocvAOStot
    n = noccEOS*noccAOS 

    !> Creating the X2ijkn
    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijpq,m,Rpqkn,k,0.0E0_realk,X2ijkn,m)

    m = noccEOS*noccEOS   ! <ij R pq> <pq R nk> = <ij X2 nk>
    k = nocvAOStot*nocvAOStot  
    n = noccEOS*noccAOS

    !> Creating the X2ijnk
    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijpq,m,Rpqnk,k,0.0E0_realk,X2ijnk,m)

    !> Creating the X2ijkl
    do i=1,noccEOS
       do j=1,noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do k=1,noccAOS
             tmp1 = tmp1 - (X2ijkn(i,j,i,k)*Fkj(k,j) + X2ijnk(i,j,k,j)*Fkj(k,i))   ! (Eq. 159)
             tmp2 = tmp2 - (X2ijkn(j,i,i,k)*Fkj(k,j) + X2ijkn(i,j,j,k)*Fkj(k,i))   ! (Eq. 169)
          enddo
          X2ijkl(i,j,j,i) =  tmp2  
          X2ijkl(i,j,i,j) =  tmp1  !> Important ordering when i=j
       enddo
    enddo

    X2energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E22(X2ijkl, noccEOS, Fragment1, Fragment2, MyFragment, X2energy, -1.0E0_realk)
    else    
       call get_mp2f12_sf_E22(X2ijkl, noccEOS, X2energy, -1.0E0_realk)
    endif

    Xenergy(2) = X2energy

    call mem_dealloc(X2ijkl)

    call mem_dealloc(X2ijkn)
    call mem_dealloc(X2ijnk)

    call mem_dealloc(Rijpq)
    call mem_dealloc(Rpqkn)
    call mem_dealloc(Rpqnk)

  end subroutine get_EX2

  !> Brief: Gives the single and pair fragment energy for X3 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EX3(Xenergy,Fkj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none

    real(realk), intent(inout) :: Xenergy(:)
    real(realk) :: X3energy, tmp1, tmp2

    real(realk), target, intent(in) :: Fkj(:,:)
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), pointer :: X3ijkl(:,:,:,:)  
    real(realk), pointer :: X3ijkn(:,:,:,:) ! sum_pq <ij|f12|mc> * <mc|f12|kn> 
    real(realk), pointer :: X3ijnk(:,:,:,:)  
 
    !real(realk), pointer :: Rijmc(:,:,:,:)  ! <ij|f12|ma'>
    !type(tensor) :: Rijmc  ! <ij|f12|ma'>
    real(realk), pointer :: Rijmc(:,:,:,:) 
    real(realk), pointer :: Rmckn(:,:,:,:)  
    real(realk), pointer :: Rmcnk(:,:,:,:)  
       
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS,noccAOStot
    integer :: ncabsMO
    integer :: m,k,n,i,j,c,l
!    integer :: dims(4)
  
    noccEOS = MyFragment%noccEOS
    noccAOS = MyFragment%noccAOS
    noccAOStot = MyFragment%nocctot
    ncabsMO = size(MyFragment%Ccabs,2) 

!    dims = [noccEOS, noccEOS, noccAOStot, ncabsMO]
!    call tensor_init(Rijmc,dims,4)
    call mem_alloc(Rijmc,  noccEOS, noccEOS, noccAOStot, ncabsMO)
    call mem_alloc(Rmckn,  noccAOStot, ncabsMO, noccEOS, noccAOS)

    !(Note INTSPEC is always stored as (2,4,1,3) )    
    ! For frozen core: index 4 (v) is only virtual
    ! Without frozen core: index 4 is core+virtual (effectively v=m) 
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iimc','RCRRG',Rijmc)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'mciv','CRRRG',Rmckn)
     
    !> Creating the X3ijkn
    m = noccEOS*noccEOS   ! <ij R mc> <mc R kn> = <ij X3 kn>
    k = noccAOStot*ncabsMO
    n = noccEOS*noccAOS
    
    call mem_alloc(X3ijkn, noccEOS, noccEOS, noccEOS, noccAOS)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc,m,Rmckn,k,0.0E0_realk,X3ijkn,m)
    !call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc%elm1,m,Rmckn,k,0.0E0_realk,X3ijkn,m)

    call mem_dealloc(Rmckn)
    call mem_alloc(Rmcnk,  noccAOStot, ncabsMO, noccAOS, noccEOS)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'mcvi','CRRRG',Rmcnk)
    ! FIXME: Generate Rmcnk by reordering Rmckn instead!!!

    call mem_alloc(X3ijnk, noccEOS, noccEOS, noccAOS, noccEOS)

    !> Creating the X3ijnk
    m = noccEOS*noccEOS   ! <ij R mc> <mc R nk> = <ij X3 nk>
    k = noccAOStot*ncabsMO
    n = noccEOS*noccAOS
    
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc,m,Rmcnk,k,0.0E0_realk,X3ijnk,m)
    !call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc%elm1,m,Rmcnk,k,0.0E0_realk,X3ijnk,m)

    !call tensor_free(Rijmc)
    call mem_dealloc(Rijmc)
    call mem_dealloc(Rmcnk)

    call mem_alloc(X3ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    !> Creating the X3ijkl
    do i=1,noccEOS
       do j=1,noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do k=1,noccAOS
             tmp1 = tmp1 - (X3ijkn(i,j,i,k)*Fkj(k,j) + X3ijnk(i,j,k,j)*Fkj(k,i))   ! (Eq. 159)
             tmp2 = tmp2 - (X3ijkn(j,i,i,k)*Fkj(k,j) + X3ijkn(i,j,j,k)*Fkj(k,i))   ! (Eq. 169)
          enddo
          X3ijkl(i,j,j,i) =  tmp2
          X3ijkl(i,j,i,j) =  tmp1
       enddo
    enddo
    call mem_dealloc(X3ijkn)
    call mem_dealloc(X3ijnk)

    X3energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E22(X3ijkl, noccEOS, Fragment1, Fragment2, MyFragment, X3energy, -1.0E0_realk)
    else    
       call get_mp2f12_sf_E22(X3ijkl, noccEOS, X3energy, -1.0E0_realk)
    endif

    Xenergy(3) = X3energy

    call mem_dealloc(X3ijkl)


  end subroutine get_EX3

  !> Brief: Gives the single and pair fragment energy for X4 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EX4(Xenergy,Fkj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none

    real(realk), intent(inout) :: Xenergy(:)
    real(realk) :: X4energy, tmp1, tmp2

    real(realk), target, intent(in) :: Fkj(:,:)
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOStot)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), pointer :: X4ijkl(:,:,:,:) 
    real(realk), pointer :: X4ijkn(:,:,:,:) 
    real(realk), pointer :: X4ijnk(:,:,:,:) ! sum_pq <ji|f12|cm> * <cm|f12|nk> = X3jink
    
    real(realk), pointer :: Rijcm(:,:,:,:)  ! <ij|f12|ma'>
    real(realk), pointer :: Rcmkn(:,:,:,:)   
    real(realk), pointer :: Rcmnk(:,:,:,:)   
  
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOStot,noccAOS
    integer :: ncabsMO
    integer :: m,k,n,i,j,c,l
  
    noccEOS = MyFragment%noccEOS
    noccAOS = MyFragment%noccAOS 
    noccAOStot = MyFragment%nocctot 
    ncabsMO = size(MyFragment%Ccabs,2) 
   
    call mem_alloc(X4ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(X4ijnk, noccEOS, noccEOS, noccAOS, noccEOS)
    call mem_alloc(X4ijkn, noccEOS, noccEOS, noccEOS, noccAOS)

    call mem_alloc(Rijcm,  noccEOS, noccEOS, ncabsMO, noccAOStot)
    call mem_alloc(Rcmkn,  ncabsMO, noccAOStot, noccEOS, noccAOS)
    call mem_alloc(Rcmnk,  ncabsMO, noccAOStot, noccAOS, noccEOS)
    
    !(Note INTSPEC is always stored as (2,4,1,3) )    
    ! For frozen core: index 4 (v) is only virtual
    ! Without frozen core: index 4 is core+virtual (effectively v=m) 
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iicm','RRRCG',Rijcm)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'cmiv','RRCRG',Rcmkn)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'cmvi','RRCRG',Rcmnk)

    !> Creating the X4ijkn
    m = noccEOS*noccEOS   ! <ij R mc> <mc R nk> = <m X4 n> instead of <ij R cm> <cm R kn> = <m X4 n>
    k = noccAOStot*ncabsMO
    n = noccEOS*noccAOS
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijcm,m,Rcmkn,k,0.0E0_realk,X4ijkn,m)
   
    !> Creating the X4ijnk
    m = noccEOS*noccEOS   ! <ij R mc> <mc R nk> = <ij X4 nk> 
    k = noccAOStot*ncabsMO
    n = noccEOS*noccAOS
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijcm,m,Rcmnk,k,0.0E0_realk,X4ijnk,m)

    !> Creating the X4ijkl
    do i=1,noccEOS
       do j=1,noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do k=1,noccAOS
             tmp1 = tmp1 - (X4ijkn(i,j,i,k)*Fkj(k,j) + X4ijnk(i,j,k,j)*Fkj(k,i))   ! (Eq. 159)
             tmp2 = tmp2 - (X4ijkn(j,i,i,k)*Fkj(k,j) + X4ijkn(i,j,j,k)*Fkj(k,i))   ! (Eq. 169)
          enddo
          X4ijkl(i,j,j,i) =  tmp2
          X4ijkl(i,j,i,j) =  tmp1
       enddo
    enddo 

    X4energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E22(X4ijkl, noccEOS, Fragment1, Fragment2, MyFragment, X4energy, -1.0E0_realk)
    else    
       call get_mp2f12_sf_E22(X4ijkl, noccEOS, X4energy, -1.0E0_realk)
    endif

    Xenergy(4) = X4energy

    call mem_dealloc(X4ijkl)
    call mem_dealloc(X4ijkn)
    call mem_dealloc(X4ijnk)

    call mem_dealloc(Rijcm)
    call mem_dealloc(Rcmkn)
    call mem_dealloc(Rcmnk)

  end subroutine get_EX4

  !> Brief: Gives the single and pair fragment energy for B1 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EB1(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none
    
    real(realk), intent(inout) :: Benergy(:)
    real(realk) :: B1energy, tmp1, tmp2

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    !> F12 integrals for the B1_term <ij|[[T,f12],f12]|kl> 
    real(realk), pointer :: B1ijkl(:,:,:,:)   

    integer :: noccEOS !number of occupied MO orbitals in EOS
   
    noccEOS = MyFragment%noccEOS
  
    call mem_alloc(B1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    
    !> B1-term
    !> B1ijkl <ij|[[T,f12],f12]|kl> stored as (i,j,k,l)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iiii','RRRRD',B1ijkl)

    B1energy = 0.0E0_realk
    
    if(dopair) then
       call get_mp2f12_pf_E23(B1ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B1energy,  1.0E0_realk)
    else
       call get_mp2f12_sf_E23(B1ijkl, noccEOS, B1energy,  1.0E0_realk)
    endif

    Benergy(1) = B1energy

    call mem_dealloc(B1ijkl)

  end subroutine get_EB1

  !> Brief: Gives the single and pair fragment energy for B2 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EB2(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none
    
    real(realk), intent(inout) :: Benergy(:)
    real(realk) :: B2energy, tmp1, tmp2

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    !> F12 integrals for the B2_term <ij|f12^2|rk>  r = RI MO   
    real(realk), pointer :: B2ijkl(:,:,:,:)     
    real(realk), pointer :: R2ijrk(:,:,:,:)
    
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: ncabsAO !number of CABS AO orbitals
    integer :: i,j,r
    
    noccEOS = MyFragment%noccEOS
    ncabsAO = size(MyFragment%Ccabs,1)    
      
    call mem_alloc(B2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(R2ijrk, noccEOS, noccEOS, ncabsAO, noccEOS)     

    !> B2-term
    !> R2ijrk <ij|f12^2|rk> stored as (i,j,r,k)    r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,&
         & CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iiri','RRRC2',R2ijrk)

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             tmp1 =  tmp1 + R2ijrk(i,j,r,j)*Myfragment%hJir(i,r)
             tmp2 =  tmp2 + R2ijrk(j,i,r,j)*Myfragment%hJir(i,r) 
          enddo
          B2ijkl(i,j,j,i) = tmp2
          B2ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    B2energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E23(B2ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B2energy,  1.0E0_realk)
    else
       call get_mp2f12_sf_E23(B2ijkl, noccEOS, B2energy,  1.0E0_realk)
    endif

    Benergy(2) = B2energy

    call mem_dealloc(B2ijkl)
    call mem_dealloc(R2ijrk)     

  end subroutine get_EB2

  !> Brief: Gives the single and pair fragment energy for B3 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EB3(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none
    
    real(realk), intent(inout) :: Benergy(:)
    real(realk) :: B3energy, tmp1, tmp2

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

   !> F12 integrals for the B3_term <ij|f12^2|kr>  r = RI MO         
    real(realk), pointer :: B3ijkl(:,:,:,:)
    real(realk), pointer :: R2ijkr(:,:,:,:)
 
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: ncabsAO !number of CABS AO orbitals
    integer :: i,j,r
    
    noccEOS = MyFragment%noccEOS
    ncabsAO = size(MyFragment%Ccabs,1)    
      
    call mem_alloc(B3ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(R2ijkr, noccEOS, noccEOS, noccEOS, ncabsAO)     

    !> B3-term
    !> R2ijkr <ij|f12^2|kr> stored as (i,j,k,r)    r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,&
         & CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iiir','RCRR2',R2ijkr)

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             tmp1 =  tmp1 + R2ijkr(i,j,i,r)*Myfragment%hJir(j,r)
             tmp2 =  tmp2 + R2ijkr(i,i,j,r)*Myfragment%hJir(j,r) 
          enddo
          B3ijkl(i,j,j,i) = tmp2
          B3ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    B3energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E23(B3ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B3energy,  1.0E0_realk)
    else
       call get_mp2f12_sf_E23(B3ijkl, noccEOS, B3energy,  1.0E0_realk)
    endif

    Benergy(3) = B3energy

    call mem_dealloc(B3ijkl)
    call mem_dealloc(R2ijkr)     

  end subroutine get_EB3

  !> Brief: Gives the single and pair fragment energy for B4 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EB4(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none
    
    real(realk), intent(inout) :: Benergy(:)
    real(realk) :: B4energy, tmp1, tmp2

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    !> F12 integrals for the B4_term
    real(realk), pointer :: B4ijkl(:,:,:,:)
    real(realk), pointer :: Rijrs(:,:,:,:)

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: ncabsAO !number of CABS AO orbitals
    integer :: i,j,r,s,t
    
    noccEOS = MyFragment%noccEOS
    ncabsAO = size(MyFragment%Ccabs,1)    
      
    call mem_alloc(B4ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijrs,  noccEOS, noccEOS, ncabsAO, ncabsAO)
    
    !> B4-term
    !> R2ijrs <ij|f12|rs> stored as (i,j,r,s)      r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iirr','RCRCG',Rijrs)

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             do s=1, ncabsAO
                do t=1, ncabsAO
                   tmp1 = tmp1 + Rijrs(i,j,r,s)*Myfragment%Krs(t,s)*Rijrs(i,j,r,t) + &
                        & Rijrs(i,j,s,r)*Myfragment%Krs(t,s)*Rijrs(i,j,t,r) 

                   tmp2 = tmp2 + Rijrs(j,i,r,s)*Myfragment%Krs(t,s)*Rijrs(i,j,r,t) + &
                        & Rijrs(j,i,s,r)*Myfragment%Krs(t,s)*Rijrs(i,j,t,r) 
                enddo
             enddo
          enddo
          B4ijkl(i,j,j,i) = tmp2
          B4ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    B4energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E23(B4ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B4energy,  1.0E0_realk)
    else
       call get_mp2f12_sf_E23(B4ijkl, noccEOS, B4energy,  1.0E0_realk)
    endif

    Benergy(4) = B4energy

    call mem_dealloc(B4ijkl)
    call mem_dealloc(Rijrs)     

  end subroutine get_EB4

  !> Brief: Gives the single and pair fragment energy for B5 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EB5(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none
    
    real(realk), intent(inout) :: Benergy(:)
    real(realk) :: B5energy, tmp1, tmp2

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    !> F12 integrals for the B5_term
    real(realk), pointer :: B5ijkl(:,:,:,:)
    real(realk), pointer :: Rijrm(:,:,:,:)

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: ncabsAO !number of CABS AO orbitals
    integer :: noccAOStot
    integer :: i,j,r,s,t,m
    
    noccEOS = MyFragment%noccEOS
    noccAOStot = MyFragment%nocctot
    ncabsAO = size(MyFragment%Ccabs,1)    

    call mem_alloc(B5ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijrm,  noccEOS, noccEOS, ncabsAO, noccAOStot)

    !> B5-term
    !> Rijrm <ij|f12|rm> stored as (i,j,r,m)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iirm','RRRCG',Rijrm)

    !> term5
    !> B5ijkl Brute force with memory savings
    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             do s=1, ncabsAO
                do m=1, noccAOStot
                   tmp1 = tmp1 + Rijrm(i,j,r,m)*Myfragment%Frs(s,r)*Rijrm(i,j,s,m) + &
                        & Rijrm(j,i,r,m)*Myfragment%Frs(s,r)*Rijrm(j,i,s,m) 
                   tmp2 = tmp2 + Rijrm(j,i,r,m)*Myfragment%Frs(s,r)*Rijrm(i,j,s,m) + &
                        & Rijrm(i,j,r,m)*Myfragment%Frs(s,r)*Rijrm(j,i,s,m)
                enddo
             enddo
          enddo
          B5ijkl(i,j,j,i) = tmp2
          B5ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    B5energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E23(B5ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B5energy, -1.0E0_realk)
    else
       call get_mp2f12_sf_E23(B5ijkl, noccEOS, B5energy, -1.0E0_realk)
    endif

    Benergy(5) = B5energy

    call mem_dealloc(B5ijkl)
    call mem_dealloc(Rijrm)

  end subroutine get_EB5

  !> Brief: Gives the single and pair fragment energy for B6 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EB6(Benergy,Fmn,Fab,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none
    
    real(realk), intent(inout) :: Benergy(:)
    real(realk) :: B6energy, tmp1, tmp2

    real(realk), target, intent(in) :: Fmn(:,:)
    real(realk), target, intent(in) :: Fab(:,:)
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    !> F12 integrals for the B6_term
    real(realk), pointer :: B6ijkl(:,:,:,:)
    real(realk), pointer :: Rijpa(:,:,:,:)

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nocvAOStot
    integer :: noccAOStot
    integer :: nvirtAOS 
    integer :: i,j,p,q,a
    integer :: offset

    if(Decinfo%Frozencore) then
       offset = MyFragment%ncore
    else
       offset = 0
    endif
    
    noccEOS = MyFragment%noccEOS
    noccAOS = MyFragment%noccAOS
    noccAOStot = MyFragment%nocctot
    nocvAOStot = MyFragment%nocctot + MyFragment%nvirtAOS
    nvirtAOS = MyFragment%nvirtAOS    

    call mem_alloc(B6ijkl, noccEOS, noccEOS,  noccEOS, noccEOS)
    call mem_alloc(Rijpa,  noccEOS, noccEOS,  nocvAOStot, nvirtAOS)
    
    !> B6-term
    !> Rijpa <ij|f12|pa> stored as (i,j,p,a)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iipa','RRRRG',Rijpa)

    !> term6
    !> Need to change this and separate this into two parts one for the Fij and one for the Fab
    !> B6ijkl Brute force with memory savings
    do i=1, noccEOS
       do j=1, noccEOS    
          tmp1  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do a=1, nvirtAOS

             do p=1, offset
                do q=1, offset
                   tmp1 = tmp1 + Rijpa(i,j,p,a)*MyFragment%ccfock(q,p)*Rijpa(i,j,q,a) + &
                      & Rijpa(j,i,p,a)*MyFragment%ccfock(q,p)*Rijpa(j,i,q,a) 

                   tmp2 = tmp2 + Rijpa(j,i,p,a)*MyFragment%ccfock(q,p)*Rijpa(i,j,q,a) + &
                      & Rijpa(i,j,p,a)*MyFragment%ccfock(q,p)*Rijpa(j,i,q,a) 
                enddo
             enddo

             do p=1+offset, noccAOStot
                do q=1+offset, noccAOStot
                   tmp1 = tmp1 + Rijpa(i,j,p,a)*MyFragment%ppfock(q-offset,p-offset)*Rijpa(i,j,q,a) + &
                      & Rijpa(j,i,p,a)*Myfragment%ppfock(q-offset,p-offset)*Rijpa(j,i,q,a) 

                   tmp2 = tmp2 + Rijpa(j,i,p,a)*MyFragment%ppfock(q-offset,p-offset)*Rijpa(i,j,q,a) + &
                      & Rijpa(i,j,p,a)*MyFragment%ppfock(q-offset,p-offset)*Rijpa(j,i,q,a) 
                enddo
             enddo

             do p=noccAOStot+1, nvirtAOS+noccAOStot
                do q=noccAOStot+1, nvirtAOS+noccAOStot
                   tmp1 = tmp1 + Rijpa(i,j,p,a)*MyFragment%qqfock(q-noccAOStot,p-noccAOStot)*Rijpa(i,j,q,a) + &
                        & Rijpa(j,i,p,a)*MyFragment%qqfock(q-noccAOStot,p-noccAOStot)*Rijpa(j,i,q,a) 

                   tmp2 = tmp2 + Rijpa(j,i,p,a)*MyFragment%qqfock(q-noccAOStot,p-noccAOStot)*Rijpa(i,j,q,a) + &
                        & Rijpa(i,j,p,a)*MyFragment%qqfock(q-noccAOStot,p-noccAOStot)*Rijpa(j,i,q,a) 
                enddo
             enddo
          
          enddo
          B6ijkl(i,j,i,j) = tmp1
          B6ijkl(i,j,j,i) = tmp2
       enddo
    enddo

    B6energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E23(B6ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B6energy, -1.0E0_realk)
    else
       call get_mp2f12_sf_E23(B6ijkl, noccEOS, B6energy, -1.0E0_realk)
    endif

    Benergy(6) = B6energy

    call mem_dealloc(B6ijkl)
    call mem_dealloc(Rijpa)

  end subroutine get_EB6

  !> Brief: Gives the single and pair fragment energy for B7 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EB7(Benergy,Fmn,Fab,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none
    
    real(realk), intent(inout) :: Benergy(:)
    real(realk) :: B7energy, tmp1, tmp2
    
    real(realk), target, intent(in) :: Fmn(:,:)
    real(realk), target, intent(in) :: Fab(:,:)
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    !> F12 integrals for the B7_term
    real(realk), pointer :: B7ijkl(:,:,:,:)
    real(realk), pointer :: Rijcm(:,:,:,:)

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: noccAOStot
    integer :: ncabsMO
    integer :: i,j,c,m,n
    integer :: offset

    if(Decinfo%Frozencore) then
       offset = MyFragment%ncore
    else
       offset = 0
    endif
    
    noccEOS = MyFragment%noccEOS
    noccAOS = MyFragment%noccAOS
    noccAOStot = MyFragment%nocctot
    ncabsMO = size(MyFragment%Ccabs,2)    

    call mem_alloc(B7ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijcm, noccEOS, noccEOS,  ncabsMO, noccAOStot)

    !> Rijcm <ij|f12|cm> stored as (i,j,c,m)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iicm','RRRCG',Rijcm)

    !> B7ijkl Brute force with memory savings
    do i=1, noccEOS
       do j=1, noccEOS    
          tmp1  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do c=1, ncabsMO
             
             do m=1, offset
                do n=1, offset
                   tmp1 = tmp1 + Rijcm(i,j,c,m)*MyFragment%ccfock(m,n)*Rijcm(i,j,c,n) + &
                        & Rijcm(j,i,c,m)*MyFragment%ccfock(m,n)*Rijcm(j,i,c,n) 

                   tmp2 = tmp2 + Rijcm(j,i,c,m)*MyFragment%ccfock(m,n)*Rijcm(i,j,c,n) + &
                        & Rijcm(i,j,c,m)*MyFragment%ccfock(m,n)*Rijcm(j,i,c,n) 
                enddo
             enddo
       
             do m=offset+1, noccAOStot
                do n=offset+1, noccAOStot
                   tmp1 = tmp1 + Rijcm(i,j,c,m)*MyFragment%ppfock(m-offset,n-offset)*Rijcm(i,j,c,n) + &
                        & Rijcm(j,i,c,m)*MyFragment%ppfock(m-offset,n-offset)*Rijcm(j,i,c,n) 

                   tmp2 = tmp2 + Rijcm(j,i,c,m)*MyFragment%ppfock(m-offset,n-offset)*Rijcm(i,j,c,n) + &
                        & Rijcm(i,j,c,m)*MyFragment%ppfock(m-offset,n-offset)*Rijcm(j,i,c,n) 
                enddo
             enddo
          
          
          enddo
          B7ijkl(i,j,j,i) = tmp2
          B7ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    B7energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E23(B7ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B7energy, 1.0E0_realk)
    else
       call get_mp2f12_sf_E23(B7ijkl, noccEOS, B7energy, 1.0E0_realk)
    endif

    Benergy(7) = B7energy

    call mem_dealloc(B7ijkl)
    call mem_dealloc(Rijcm)

  end subroutine get_EB7

  !> Brief: Gives the single and pair fragment energy for B7 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EB8(Benergy,Fmn,Fab,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none
    
    real(realk), intent(inout) :: Benergy(:)
    real(realk) :: B8energy, tmp1, tmp2

    real(realk), target, intent(in) :: Fmn(:,:)
    real(realk), target, intent(in) :: Fab(:,:)
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    !> F12 integrals for the B8_term
    real(realk), pointer :: B8ijkl(:,:,:,:)
    real(realk), pointer :: Rijcm(:,:,:,:)
    real(realk), pointer :: Rijcr(:,:,:,:) 
 
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: noccAOStot
    integer :: ncabsMO
    integer :: ncabsAO
    integer :: i,j,r,c,m
    integer :: offset

    noccEOS = MyFragment%noccEOS
    noccAOS = MyFragment%noccAOS
    noccAOStot = MyFragment%nocctot
    ncabsMO = size(MyFragment%Ccabs,2)    
    ncabsAO = size(MyFragment%Ccabs,1)    

    if(Decinfo%Frozencore) then
       offset = MyFragment%ncore
    else
       offset = 0
    endif

    call mem_alloc(B8ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijcr,  noccEOS, noccEOS, ncabsMO, ncabsAO) 
    call mem_alloc(Rijcm,  noccEOS, noccEOS, ncabsMO, noccAOStot) 
 
    !> Rijcm <ij|f12|cm> stored as (i,j,c,m)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iicm','RRRCG',Rijcm)
    !> Rijcr <ij|f12|cr> stored as (i,j,c,r)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iicr','RCRCG',Rijcr)

    !> B8ijkl Brute force with memory savings
    do i=1, noccEOS
       do j=1, noccEOS    
          tmp1 = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             do c=1, ncabsMO
                
                do m=1, noccAOStot
                   tmp1 = tmp1 + Rijcm(i,j,c,m)*MyFragment%Frm(r,m)*Rijcr(i,j,c,r) + &
                        & Rijcm(j,i,c,m)*MyFragment%Frm(r,m)*Rijcr(j,i,c,r) 

                   tmp2 = tmp2 + Rijcm(j,i,c,m)*MyFragment%Frm(r,m)*Rijcr(i,j,c,r) + &
                        & Rijcm(i,j,c,m)*MyFragment%Frm(r,m)*Rijcr(j,i,c,r) 
                enddo

             enddo
          enddo
          B8ijkl(i,j,j,i) = tmp2
          B8ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    !call ls_output(MyFragment%Frm,1,ncabsAO,1,noccAOStot,ncabsAO,noccAOStot,1,6)
   
    B8energy = 0.0E0_realk

    if(dopair) then
       call get_mp2f12_pf_E23(B8ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B8energy, -2.0E0_realk)
    else
       call get_mp2f12_sf_E23(B8ijkl, noccEOS, B8energy, -2.0E0_realk)
    endif

    Benergy(8) = B8energy

    call mem_dealloc(B8ijkl)
    call mem_dealloc(Rijcm)
    call mem_dealloc(Rijcr)
   
  end subroutine get_EB8

  !> Brief: Gives the single and pair fragment energy for B9 term in MP2F12
  !> Author: Yang M. Wang
  !> Date: August 2014
  subroutine get_EB9(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    implicit none
    
    real(realk), intent(inout) :: Benergy(:)
    real(realk) :: B9energy, tmp1, tmp2

    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    !> F12 integrals for the B9_term   
    real(realk), pointer :: B9ijkl(:,:,:,:)                                                                     
    real(realk), pointer :: Rijca(:,:,:,:)
    real(realk), pointer :: Rijpa(:,:,:,:)
 
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: ncabsMO
    integer :: nvirtAOS
    integer :: nocvAOStot
    integer :: i,j,r,c,a,p

    noccEOS  = MyFragment%noccEOS
    nocvAOStot = MyFragment%nocctot + MyFragment%nvirtAOS
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO = size(MyFragment%Ccabs,2)

    call mem_alloc(B9ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijpa,  noccEOS, noccEOS, nocvAOStot, nvirtAOS)
    call mem_alloc(Rijca,  noccEOS, noccEOS, ncabsMO, nvirtAOS)

    !> Rijpa <ij|f12|pa> stored as (i,j,p,a)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iipa','RRRRG',Rijpa)
    !> Rijca <ij|f12|ca> stored as (i,j,c,a)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,CocvAOStot,Ccabs,Cri,CvirtAOS,'iica','RRRCG',Rijca)

    do i=1, noccEOS
       do j=1, noccEOS    
          tmp1 = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do c=1, ncabsMO
             do a=1, nvirtAOS
                do p=1, nocvAOStot
                   tmp1 = tmp1 + Rijpa(i,j,p,a)*MyFragment%Fcp(c,p)*Rijca(i,j,c,a) + &
                        & Rijpa(j,i,p,a)*MyFragment%Fcp(c,p)*Rijca(j,i,c,a)

                   tmp2 = tmp2 + Rijpa(j,i,p,a)*MyFragment%Fcp(c,p)*Rijca(i,j,c,a) + &
                        & Rijpa(i,j,p,a)*MyFragment%Fcp(c,p)*Rijca(j,i,c,a)
                enddo
             enddo
          enddo
          B9ijkl(i,j,j,i) = tmp2
          B9ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    call mem_dealloc(Rijpa)
    call mem_dealloc(Rijca)

    B9energy = 0.0E0_realk
    if(dopair) then
       call get_mp2f12_pf_E23(B9ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B9energy, -2.0E0_realk)
    else
       call get_mp2f12_sf_E23(B9ijkl, noccEOS, B9energy, -2.0E0_realk)
    endif
    Benergy(9) = B9energy

    call mem_dealloc(B9ijkl)
 
  end subroutine get_EB9

  !> CCSD Routines

  !> Date: August 2014
  subroutine ccsdf12_Vijab_EV1(Venergy, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
       & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Taibj)

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V1energy, tmp1, tmp2
  
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Taibj(:,:,:,:)  
    real(realk), pointer :: Tabij(:,:,:,:) 
    
    real(realk), pointer :: Fijab(:,:,:,:)  
    real(realk), pointer :: V1ijkl(:,:,:,:)     !V1_ijkl <ij|f12*r^-1|kl>

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c

    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(Fijab, noccEOS,  noccEOS,  nvirtAOS, nvirtAOS) 
 
    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iiaa','RRRRF', Fijab)

    V1energy = 0.0E0_realk
    
    call mem_alloc(V1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Tabij, nvirtAOS, nvirtAOS,  noccEOS,  noccEOS)
    
    call array_reorder_4d(1.0E0_realk,Taibj,nvirtAOS,noccEOS,nvirtAOS,noccEOS,[1,3,2,4],0.0E0_realk,Tabij)

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do a=1, nvirtAOS
             do b=1, nvirtAOS
                tmp1 = tmp1 + Fijab(i,j,a,b)*Tabij(a,b,i,j)
                tmp2 = tmp2 + Fijab(i,j,a,b)*Tabij(a,b,j,i)
             enddo
          enddo
          V1ijkl(i,j,j,i) = tmp2
          V1ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V1ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V1energy, 1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V1ijkl, noccEOS, V1energy,  1.0E0_realk)
    endif

    Venergy(1) = V1energy

    call mem_dealloc(Fijab)
    call mem_dealloc(V1ijkl)
    call mem_dealloc(Tabij)

  end subroutine ccsdf12_Vijab_EV1

  !> Date: August 2014
  subroutine ccsdf12_Vijab_EV2(Venergy, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
       & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Taibj)

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V2energy, tmp1, tmp2
  
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Taibj(:,:,:,:)  
    real(realk), pointer :: Tabij(:,:,:,:) 
    
    real(realk), pointer :: Rijpq(:,:,:,:)  
    real(realk), pointer :: Gpqab(:,:,:,:)    
   
    real(realk), pointer :: V2ijab(:,:,:,:)   
    real(realk), pointer :: V2ijkl(:,:,:,:)   

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c, p, q

    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(Rijpq, noccEOS,  noccEOS,  nocvAOS,  nocvAOS) 
    call mem_alloc(Gpqab, nocvAOS,  nocvAOS,  nvirtAOS, nvirtAOS) 

    call mem_alloc(V2ijab, noccEOS,  noccEOS, nvirtAOS, nvirtAOS)
    call mem_alloc(V2ijkl, noccEOS,  noccEOS, noccEOS, noccEOS)
   
    call mem_alloc(Tabij, nvirtAOS, nvirtAOS, noccEOS, noccEOS)

    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iipp','RRRRG', Rijpq)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'ppaa','RRRRC', Gpqab)
    
    m = noccEOS*noccEOS  ! <ij R pq> <pq G ab> = <m V2 n> 
    k = nocvAOS*nocvAOS  
    n = nvirtAOS*nvirtAOS


    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijpq,m,Gpqab,k,0.0E0_realk,V2ijab,m)   
    
    V2energy = 0.0E0_realk
    
    call array_reorder_4d(1.0E0_realk,Taibj,nvirtAOS,noccEOS,nvirtAOS,noccEOS,[1,3,2,4],0.0E0_realk,Tabij)


    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do a=1, nvirtAOS
             do b=1, nvirtAOS
                tmp1 = tmp1 + V2ijab(i,j,a,b)*Tabij(a,b,i,j)
                tmp2 = tmp2 + V2ijab(i,j,a,b)*Tabij(a,b,j,i)
             enddo
          enddo
          V2ijkl(i,j,j,i) = tmp2
          V2ijkl(i,j,i,j) = tmp1
       enddo
    enddo


    if(dopair) then
       call get_mp2f12_pf_E21(V2ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V2energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V2ijkl, noccEOS, V2energy, -1.0E0_realk)
    endif

    Venergy(2) = V2energy

    call mem_dealloc(V2ijab)
    call mem_dealloc(V2ijkl)

    call mem_dealloc(Rijpq)
    call mem_dealloc(Gpqab)

    call mem_dealloc(Tabij)

  end subroutine ccsdf12_Vijab_EV2

  !> Date: August 2014
  subroutine ccsdf12_Vijab_EV3(Venergy, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
       & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Taibj)

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V3energy, V4energy, tmp1, tmp2
  
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Taibj(:,:,:,:)  
    real(realk), pointer :: Tabij(:,:,:,:) 
    
    real(realk), pointer :: Rijmc(:,:,:,:)  
    real(realk), pointer :: Gmcab(:,:,:,:)    

    real(realk), pointer :: V3ijkl(:,:,:,:)   
    real(realk), pointer :: V4ijkl(:,:,:,:)   


    real(realk), pointer :: V3ijab(:,:,:,:)   
    
    real(realk), pointer :: V4ijab(:,:,:,:)   

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c, p, q

    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(Rijmc, noccEOS, noccEOS, noccAOS, ncabsMO) 
    call mem_alloc(Gmcab, noccAOS, ncabsMO, nvirtAOS, nvirtAOS) 

    call mem_alloc(V3ijab, noccEOS,  noccEOS, nvirtAOS, nvirtAOS)
    call mem_alloc(V3ijkl, noccEOS,  noccEOS, noccEOS, noccEOS)

    call mem_alloc(V4ijab, noccEOS,  noccEOS, nvirtAOS, nvirtAOS)
    call mem_alloc(V4ijkl, noccEOS,  noccEOS, noccEOS, noccEOS)

    call mem_alloc(Tabij, nvirtAOS, nvirtAOS, noccEOS, noccEOS)

    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iimc','RCRRG', Rijmc)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'mcaa','CRRRC', Gmcab)
    
    m = noccEOS*noccEOS  ! <ij R pq> <pq G ab> = <m V2 n> 
    k = noccAOS*ncabsMO  
    n = nvirtAOS*nvirtAOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc,m,Gmcab,k,0.0E0_realk,V3ijab,m)   

    V3energy = 0.0E0_realk
    
    call array_reorder_4d(1.0E0_realk,Taibj,nvirtAOS,noccEOS,nvirtAOS,noccEOS,[1,3,2,4],0.0E0_realk,Tabij)
    call array_reorder_4d(1.0E0_realk,V3ijab,noccEOS,noccEOS,nvirtAOS,nvirtAOS,[2,1,4,3],0.0E0_realk,V4ijab)

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do a=1, nvirtAOS
             do b=1, nvirtAOS
                tmp1 = tmp1 + V3ijab(i,j,a,b)*Tabij(a,b,i,j)
                tmp2 = tmp2 + V3ijab(i,j,a,b)*Tabij(a,b,j,i)
             enddo
          enddo
          V3ijkl(i,j,j,i) = tmp2
          V3ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do a=1, nvirtAOS
             do b=1, nvirtAOS
                tmp1 = tmp1 + V4ijab(i,j,a,b)*Tabij(a,b,i,j)
                tmp2 = tmp2 + V4ijab(i,j,a,b)*Tabij(a,b,j,i)
             enddo
          enddo
          V4ijkl(i,j,j,i) = tmp2
          V4ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V3ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V3energy, -1.0E0_realk)
       call get_mp2f12_pf_E21(V4ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V4energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V3ijkl, noccEOS, V3energy, -1.0E0_realk)
       call get_mp2f12_sf_E21(V4ijkl, noccEOS, V4energy, -1.0E0_realk)
    endif

    Venergy(3) = V3energy
    Venergy(4) = V4energy

    call mem_dealloc(V3ijkl)
    call mem_dealloc(V4ijkl)

    call mem_dealloc(V3ijab)
    call mem_dealloc(V4ijab)

    call mem_dealloc(Rijmc)
    call mem_dealloc(Gmcab)

    call mem_dealloc(Tabij)

  end subroutine ccsdf12_Vijab_EV3
  
  !> Date: Okt 2014
  subroutine ccsdf12_Vijia_EV1(Venergy, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
       & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V1energy, tmp1, tmp2
  
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Tai(:,:)  
   
    real(realk), pointer :: Fijka(:,:,:,:)  
    real(realk), pointer :: Fijak(:,:,:,:)  
    real(realk), pointer :: V1ijkl(:,:,:,:)     !V1_ijkl <ij|f12*r^-1|kl>

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c

    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(Fijka, noccEOS, noccEOS, noccEOS, nvirtAOS) 
    call mem_alloc(Fijak, noccEOS, noccEOS, nvirtAOS, noccEOS) 
    call mem_alloc(V1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)

    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iiia','RRRRF', Fijka)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iiai','RRRRF', Fijak)
    
    V1energy = 0.0E0_realk
  
    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
       !   do k=1, noccEOS !HB Fijka using the fixed amplitude here and in addition in the mp2f12 energy
             do a=1, nvirtAOS
                tmp1 = tmp1 + Fijka(i,j,i,a)*Tai(a,j)
                tmp2 = tmp2 + Fijak(i,j,a,i)*Tai(a,j)
             enddo
      !    enddo
          V1ijkl(i,j,j,i) = tmp2
          V1ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V1ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V1energy, 1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V1ijkl, noccEOS, V1energy,  1.0E0_realk)
    endif

    Venergy(1) = V1energy

    call mem_dealloc(Fijka)
    call mem_dealloc(Fijak)
    call mem_dealloc(V1ijkl)

  end subroutine ccsdf12_Vijia_EV1

  subroutine ccsdf12_Vijia_EV2(Venergy, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,&
       & CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V2energy, tmp1, tmp2
  
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Tai(:,:)  
   
    real(realk), pointer :: Rijpq(:,:,:,:)  
    real(realk), pointer :: Gpqia(:,:,:,:)  
    real(realk), pointer :: Gpqai(:,:,:,:)  
    real(realk), pointer :: V2ijkl(:,:,:,:)     
    real(realk), pointer :: V2ijka(:,:,:,:)     
    real(realk), pointer :: V2ijak(:,:,:,:)
    
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c, p, q
    
    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(Rijpq,  noccEOS, noccEOS, nocvAOS,  nocvAOS)
    call mem_alloc(Gpqia,  nocvAOS, nocvAOS, noccEOS,  nvirtAOS)
    call mem_alloc(Gpqai,  nocvAOS, nocvAOS, nvirtAOS, noccEOS)
    call mem_alloc(V2ijkl, noccEOS, noccEOS, noccEOS,  noccEOS)
    call mem_alloc(V2ijka, noccEOS, noccEOS, noccEOS,  nvirtAOS)
    call mem_alloc(V2ijak, noccEOS, noccEOS, nvirtAOS, noccEOS)

    !> NB NB NB NB NB !> Thomas code G is F12 and R is Coulomb, in my code G is Coulomb and R is F12
    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         &CocvAOStot,Ccabs,Cri,CvirtAOS,'ppia','RRRRC', Gpqia)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iipp','RRRRG', Rijpq)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'ppai','RRRRC', Gpqai)
    
    V2energy = 0.0E0_realk

    m = noccEOS*noccEOS  ! <ij R pq> <pq G ia> = <m V2 n> 
    k = nocvAOS*nocvAOS  
    n = noccEOS*nvirtAOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijpq,m,Gpqia,k,0.0E0_realk,V2ijka,m)   

    m = noccEOS*noccEOS  ! <ij R pq> <pq G ai> = <m V2 n> 
    k = nocvAOS*nocvAOS  
    n = noccEOS*nvirtAOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijpq,m,Gpqai,k,0.0E0_realk,V2ijak,m)   
    
    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do a=1, nvirtAOS
             tmp1 = tmp1 + V2ijka(i,j,i,a)*Tai(a,j)
             tmp2 = tmp2 + V2ijak(i,j,a,i)*Tai(a,j)
          enddo
          V2ijkl(i,j,j,i) = tmp2
          V2ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V2ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V2energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V2ijkl, noccEOS, V2energy, -1.0E0_realk)
    endif

    Venergy(2) = V2energy

    call mem_dealloc(Gpqia)
    call mem_dealloc(Gpqai)
    call mem_dealloc(Rijpq)
    call mem_dealloc(V2ijkl)
    call mem_dealloc(V2ijka)
    call mem_dealloc(V2ijak)

  end subroutine ccsdf12_Vijia_EV2  

  subroutine ccsdf12_Vijia_EV3(Venergy, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
       & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V3energy, V4energy, tmp1, tmp2
  
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Tai(:,:)  
   
    real(realk), pointer :: Rijmc(:,:,:,:)  
    real(realk), pointer :: Gmcia(:,:,:,:)  
    real(realk), pointer :: Gmcai(:,:,:,:)  

    real(realk), pointer :: V3ijkl(:,:,:,:)     
    real(realk), pointer :: V3ijka(:,:,:,:)     
    real(realk), pointer :: V3ijak(:,:,:,:)

    real(realk), pointer :: V4ijkl(:,:,:,:)     
    real(realk), pointer :: V4jiak(:,:,:,:)
    real(realk), pointer :: V4jika(:,:,:,:)
    
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c

    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(Rijmc,  noccEOS, noccEOS, noccAOS,  ncabsMO)
    call mem_alloc(Gmcia,  noccAOS, ncabsMO, noccEOS,  nvirtAOS)
    call mem_alloc(Gmcai,  noccAOS, ncabsMO, nvirtAOS, noccEOS)
 
    call mem_alloc(V3ijka, noccEOS, noccEOS, noccEOS,  nvirtAOS)
    call mem_alloc(V3ijak, noccEOS, noccEOS, nvirtAOS, noccEOS)
    call mem_alloc(V3ijkl, noccEOS, noccEOS, noccEOS,  noccEOS)

    call mem_alloc(V4jika, noccEOS, noccEOS, noccEOS,  nvirtAOS)
    call mem_alloc(V4jiak, noccEOS, noccEOS, nvirtAOS, noccEOS)
    call mem_alloc(V4ijkl, noccEOS, noccEOS, noccEOS,  noccEOS)  

    !> NB NB NB NB NB !> Thomas code G is F12 and R is Coulomb, in my code G is Coulomb and R is F12
    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3) 
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iimc','RCRRG', Rijmc)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'mcia','CRRRC', Gmcia)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'mcai','CRRRC', Gmcai)
    
    V3energy = 0.0E0_realk
    V4energy = 0.0E0_realk

    m = noccEOS*noccEOS  ! <ij R mc> <mc G ia> = <m V3 n> 
    k = noccAOS*ncabsMO  
    n = noccEOS*nvirtAOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc,m,Gmcia,k,0.0E0_realk,V3ijka,m)   

    m = noccEOS*noccEOS  ! <ij R mc> <mc G ai> = <m V3 n> 
    k = noccAOS*ncabsMO  
    n = noccEOS*nvirtAOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc,m,Gmcai,k,0.0E0_realk,V3ijak,m)   

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
       !   do k=1, noccEOS !HB Fijka using the fixed amplitude here and in addition in the mp2f12 energy
             do a=1, nvirtAOS
                tmp1 = tmp1 + V3ijka(i,j,i,a)*Tai(a,j)
                tmp2 = tmp2 + V3ijak(i,j,a,i)*Tai(a,j)
             enddo
      !    enddo
          V3ijkl(i,j,j,i) = tmp2
          V3ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V3ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V3energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V3ijkl, noccEOS, V3energy, -1.0E0_realk)
    endif

    Venergy(3) = V3energy

    !> Venergy(4)
    !> Creating V4jika 
    call array_reorder_4d(1.0E0_realk,V3ijka,noccEOS,noccEOS,noccEOS,nvirtAOS,[2,1,3,4],0.0E0_realk,V4jika)
    
    !> Creating V4jiak 
    call array_reorder_4d(1.0E0_realk,V3ijak,noccEOS,noccEOS,nvirtAOS,noccEOS,[2,1,3,4],0.0E0_realk,V4jiak)

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
       !   do k=1, noccEOS !HB Fijka using the fixed amplitude here and in addition in the mp2f12 energy
             do a=1, nvirtAOS
               tmp1 = tmp1 + V4jiak(i,j,a,i)*Tai(a,j)
               tmp2 = tmp2 + V4jika(i,j,i,a)*Tai(a,j)
             enddo
      !    enddo
          V4ijkl(i,j,j,i) = tmp2
          V4ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V4ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V4energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V4ijkl, noccEOS, V4energy, -1.0E0_realk)
    endif

    Venergy(4) = V4energy
  
    call mem_dealloc(Gmcia)
    call mem_dealloc(Gmcai)
    call mem_dealloc(Rijmc)

    call mem_dealloc(V3ijkl)
    call mem_dealloc(V3ijka)
    call mem_dealloc(V3ijak)

    call mem_dealloc(V4jika)
    call mem_dealloc(V4jiak)
    call mem_dealloc(V4ijkl)

  end subroutine ccsdf12_Vijia_EV3


  !> Date: Okt 2014
  subroutine ccsdf12_Viajj_EV1(Venergy, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
       &CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V1energy, tmp1, tmp2
  
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Tai(:,:)  
   
    real(realk), pointer :: Fijka(:,:,:,:)  
    real(realk), pointer :: Fijak(:,:,:,:)  
    real(realk), pointer :: V1ijkl(:,:,:,:)     !V1_ijkl <ij|f12*r^-1|kl>

    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c

    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(V1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Fijka, noccEOS, noccEOS, noccEOS, nvirtAOS) 
    call mem_alloc(Fijak, noccEOS, noccEOS, nvirtAOS, noccEOS) 

    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iiia','RRRRF', Fijka)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iiai','RRRRF', Fijak)
    
    V1energy = 0.0E0_realk
   
    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
       !   do k=1, noccEOS !HB Fijka using the fixed amplitude here and in addition in the mp2f12 energy
             do a=1, nvirtAOS
                tmp1 = tmp1 + Fijak(i,j,a,j)*Tai(a,i)
                tmp2 = tmp2 + Fijka(i,j,j,a)*Tai(a,i)
             enddo
      !    enddo
          V1ijkl(i,j,j,i) = tmp2
          V1ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V1ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V1energy, 1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V1ijkl, noccEOS, V1energy,  1.0E0_realk)
    endif

    Venergy(1) = V1energy

    call mem_dealloc(Fijka)
    call mem_dealloc(Fijak)
    call mem_dealloc(V1ijkl)
    
  end subroutine ccsdf12_Viajj_EV1

  subroutine ccsdf12_Viajj_EV2(Venergy, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
       & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V2energy, tmp1, tmp2
  
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Tai(:,:)  
   
    real(realk), pointer :: Rijpq(:,:,:,:)  
    real(realk), pointer :: Gpqia(:,:,:,:)  
    real(realk), pointer :: Gpqai(:,:,:,:)  

    real(realk), pointer :: V2ijkl(:,:,:,:)     
    real(realk), pointer :: V2ijka(:,:,:,:)     
    real(realk), pointer :: V2ijak(:,:,:,:)
  
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c, p, q

    real(realk), pointer :: tmp(:)
    
    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(Rijpq,  noccEOS, noccEOS, nocvAOS,  nocvAOS)
    call mem_alloc(Gpqia,  nocvAOS, nocvAOS, noccEOS,  nvirtAOS)
    call mem_alloc(Gpqai,  nocvAOS, nocvAOS, nvirtAOS, noccEOS)
 
   call mem_alloc(V2ijkl, noccEOS, noccEOS, noccEOS,  noccEOS)
    call mem_alloc(V2ijka, noccEOS, noccEOS, noccEOS,  nvirtAOS)
    call mem_alloc(V2ijak, noccEOS, noccEOS, nvirtAOS, noccEOS)

    !> NB NB NB NB NB !> Thomas code G is F12 and R is Coulomb, in my code G is Coulomb and R is F12
    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iipp','RRRRG', Rijpq)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'ppia','RRRRC', Gpqia)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'ppai','RRRRC', Gpqai)
    
    V2energy = 0.0E0_realk

    !> Creating Gqpai 
    !call array_reorder_4d(1.0E0_realk,Gpqia,nocvAOS,nocvAOS,nocvAOS,nvirtAOS,[2,1,3,4],0.0E0_realk,Gqpia)
   
    m = noccEOS*noccEOS  ! <ij R pq> <pq G ia> = <m V2 n> 
    k = nocvAOS*nocvAOS  
    n = noccEOS*nvirtAOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijpq,m,Gpqai,k,0.0E0_realk,V2ijak,m)   

    m = noccEOS*noccEOS  ! <ij R pq> <pq G ai> = <m V2 n> 
    k = nocvAOS*nocvAOS  
    n = noccEOS*nvirtAOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijpq,m,Gpqia,k,0.0E0_realk,V2ijka,m)   
    
    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
          do a=1, nvirtAOS
             tmp1 = tmp1 + V2ijak(i,j,a,j)*Tai(a,i)
             tmp2 = tmp2 + V2ijka(i,j,j,a)*Tai(a,i)
          enddo
          V2ijkl(i,j,j,i) = tmp2
          V2ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V2ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V2energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V2ijkl, noccEOS, V2energy, -1.0E0_realk)
    endif

    Venergy(2) = V2energy
  
    call mem_dealloc(Gpqia)
    call mem_dealloc(Gpqai)
    call mem_dealloc(Rijpq)

    call mem_dealloc(V2ijkl)
    call mem_dealloc(V2ijka)
    call mem_dealloc(V2ijak)

  end subroutine ccsdf12_Viajj_EV2


  subroutine ccsdf12_Viajj_EV3(Venergy, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
       & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)

    real(realk), intent(inout) :: Venergy(:)
    real(realk) :: V3energy, V4energy, tmp1, tmp2
  
    type(decfrag),intent(inout) :: MyFragment
    type(decfrag),intent(in) :: Fragment1
    type(decfrag),intent(in) :: Fragment2

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    logical, intent(in) :: dopair 

    real(realk), intent(in), optional :: Tai(:,:)  
   
    real(realk), pointer :: Rijmc(:,:,:,:)  
    real(realk), pointer :: Gmcia(:,:,:,:)  
    real(realk), pointer :: Gmcai(:,:,:,:)  

    real(realk), pointer :: V3ijkl(:,:,:,:)     
    real(realk), pointer :: V3ijka(:,:,:,:)     
    real(realk), pointer :: V3ijak(:,:,:,:)

    real(realk), pointer :: V4ijkl(:,:,:,:)     
    real(realk), pointer :: V4jiak(:,:,:,:)
    real(realk), pointer :: V4jika(:,:,:,:)
    
    integer :: noccEOS !number of occupied MO orbitals in EOS
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c

    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(Rijmc,  noccEOS, noccEOS, noccAOS,  ncabsMO)
    call mem_alloc(Gmcia,  noccAOS, ncabsMO, noccEOS,  nvirtAOS)
    call mem_alloc(Gmcai,  noccAOS, ncabsMO, nvirtAOS, noccEOS)
 
    call mem_alloc(V3ijka, noccEOS, noccEOS, noccEOS,  nvirtAOS)
    call mem_alloc(V3ijak, noccEOS, noccEOS, nvirtAOS, noccEOS)
    call mem_alloc(V3ijkl, noccEOS, noccEOS, noccEOS,  noccEOS)

    call mem_alloc(V4jika, noccEOS, noccEOS, noccEOS,  nvirtAOS)
    call mem_alloc(V4jiak, noccEOS, noccEOS, nvirtAOS, noccEOS)
    call mem_alloc(V4ijkl, noccEOS, noccEOS, noccEOS,  noccEOS)  

    !> NB NB NB NB NB !> Thomas code G is F12 and R is Coulomb, in my code G is Coulomb and R is F12
    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3) 
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iimc','RCRRG', Rijmc)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'mcia','CRRRC', Gmcia)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'mcai','CRRRC', Gmcai)
    
    V3energy = 0.0E0_realk
    V4energy = 0.0E0_realk

    m = noccEOS*noccEOS  ! <ij R mc> <mc G ia> = <m V3 n> 
    k = noccAOS*ncabsMO  
    n = noccEOS*nvirtAOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc,m,Gmcai,k,0.0E0_realk,V3ijak,m)   

    m = noccEOS*noccEOS  ! <ij R mc> <mc G ai> = <m V3 n> 
    k = noccAOS*ncabsMO  
    n = noccEOS*nvirtAOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc,m,Gmcia,k,0.0E0_realk,V3ijka,m)   

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
       !   do k=1, noccEOS !HB Fijka using the fixed amplitude here and in addition in the mp2f12 energy
             do a=1, nvirtAOS
                tmp1 = tmp1 + V3ijak(i,j,a,j)*Tai(a,i)
                tmp2 = tmp2 + V3ijka(i,j,j,a)*Tai(a,i)
             enddo
      !    enddo
          V3ijkl(i,j,j,i) = tmp2
          V3ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V3ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V3energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V3ijkl, noccEOS, V3energy, -1.0E0_realk)
    endif

    Venergy(3) = V3energy

    !> Venergy(4)
    !> Creating V4jika 
    call array_reorder_4d(1.0E0_realk,V3ijka,noccEOS,noccEOS,noccEOS,nvirtAOS,[2,1,3,4],0.0E0_realk,V4jika)
    
    !> Creating V4jiak 
    call array_reorder_4d(1.0E0_realk,V3ijak,noccEOS,noccEOS,nvirtAOS,noccEOS,[2,1,3,4],0.0E0_realk,V4jiak)

    do i=1, noccEOS
       do j=1, noccEOS
          tmp1 = 0E0_realk
          tmp2 = 0E0_realk
       !   do k=1, noccEOS !HB Fijka using the fixed amplitude here and in addition in the mp2f12 energy
             do a=1, nvirtAOS
               tmp1 = tmp1 + V4jika(i,j,j,a)*Tai(a,i)
               tmp2 = tmp2 + V4jiak(i,j,a,j)*Tai(a,i)
             enddo
      !    enddo
          V4ijkl(i,j,j,i) = tmp2
          V4ijkl(i,j,i,j) = tmp1
       enddo
    enddo

    if(dopair) then
       call get_mp2f12_pf_E21(V4ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V4energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V4ijkl, noccEOS, V4energy, -1.0E0_realk)
    endif

    Venergy(4) = V4energy
  
    call mem_dealloc(Gmcia)
    call mem_dealloc(Gmcai)
    call mem_dealloc(Rijmc)

    call mem_dealloc(V3ijkl)
    call mem_dealloc(V3ijka)
    call mem_dealloc(V3ijak)

    call mem_dealloc(V4jika)
    call mem_dealloc(V4jiak)
    call mem_dealloc(V4ijkl)

  end subroutine ccsdf12_Viajj_EV3

  subroutine get_ccsd_energy(CCSDenergy,MyFragment,CoccEOS,CoccAOStot,&
       & CvirtAOS,CocvAOStot,Ccabs,Cri,Tai,Taibj)

    real(realk), intent(inout) :: CCSDenergy
    type(decfrag),intent(inout) :: MyFragment

    real(realk), target, intent(in) :: CoccEOS(:,:)  !CoccEOS(nbasis,noccEOS)
    real(realk), target, intent(in) :: CoccAOStot(:,:)  !CoccEOS(nbasis,noccAOS)
    real(realk), target, intent(in) :: CocvAOStot(:,:)  !CocvAOStot(nbasis, nocvAOS)
    real(realk), target, intent(in) :: Ccabs(:,:)    !Ccabs(ncabsAO, ncabsMO)
    real(realk), target, intent(in) :: Cri(:,:)      !Cri(ncabsAO,ncabsAO)
    real(realk), target, intent(in) :: CvirtAOS(:,:) !CvritAOS(nbasis,nvirtAOS)

    real(realk), intent(in) :: Tai(:,:)  
    real(realk), intent(in) :: Taibj(:,:,:,:)     
   
    real(realk), pointer :: Gijab(:,:,:,:)  

    real(realk) :: Ttmp, Gtmp
    
    integer :: noccEOS 
    integer :: noccAOS
    integer :: nvirtAOS
    integer :: ncabsMO
    integer :: nocvAOS
    integer :: i, j, m, k, n, a, b, c

    noccEOS  = MyFragment%noccEOS
    noccAOS  = MyFragment%noccAOS 
    nvirtAOS = MyFragment%nvirtAOS
    ncabsMO  = size(MyFragment%Ccabs,2)
    nocvAOS  = MyFragment%noccAOS + MyFragment%nvirtAOS

    call mem_alloc(Gijab,  noccEOS, noccEOS, nvirtAOS, nvirtAOS)
 
    !> NB NB NB NB NB !> Thomas code G is F12 and R is Coulomb, in my code G is Coulomb and R is F12
    !> Get integrals <ij|f12*r^-1|ab> stored as (i,j,a,b)  (Note INTSPEC is always stored as (2,4,1,3))      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3) 
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOStot,&
         & CocvAOStot,Ccabs,Cri,CvirtAOS,'iiaa','RRRRC', Gijab)
    
    ! Calculate standard CCSD energy (brainless summation in this test code)
    CCSDenergy=0.0E0_realk
    do j=1,noccEOS
       do b=1,nvirtAOS
          do i=1,noccEOS
             do a=1,nvirtAOS
                ! Energy = sum_{ijab} ( Tai*Tbj + Taibj) * (ai | bj)
                Ttmp = Tai(a,i)*Tai(b,j) + Taibj(a,i,b,j)
                Gtmp = 2.0E0_realk * Gijab(i,j,a,b) - Gijab(i,j,b,a)
                CCSDenergy = CCSDenergy + Ttmp * Gtmp
             end do
          end do
       end do
    end do

    call mem_dealloc(Gijab)

  end subroutine get_ccsd_energy

  !> Brief: Gives the single and pair fragment energy for MP2F12
  !> Author: Yang M. Wang
  !> Date: April 2013
  subroutine get_f12_fragment_energy(MyFragment, Taibj, Tai, case, Fragment1, Fragment2)
    implicit none

    !> Atomic fragment to be determined (Single or Pair fragment)
    type(decfrag), intent(inout) :: MyFragment
    !> t2EOS amplitudes stored in the order T(a,i,b,j)
    real(realk), intent(in), pointer :: Taibj(:,:,:,:) 
    !> t1EOS amplitudes stored in the order T(a,i)
    real(realk), intent(in), pointer :: Tai(:,:) 
    !> Case MODEL
    integer, intent(in) :: case
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(in), optional :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(in), optional :: Fragment2
    !> Logical variable to check if this is a pair fragment
    logical :: dopair

    ! ***********************************************************
    !   Allocating for Coefficient matrix
    ! ***********************************************************
    !> MO coefficient matrix for the occupied EOS
    real(realk), pointer :: CoccEOS(:,:)
    !> MO coefficient matrix for the occupied AOS (core + valence, both with and without frozen core)
    real(realk), pointer :: CoccAOStot(:,:)
    !> MO coefficient matrix for the virtual AOS
    real(realk), pointer :: CvirtAOS(:,:)
    !> MO coefficient matrix for the occupied + virtual AOS
    real(realk), pointer :: CocvAOStot(:,:)
    !> MO coefficient matrix for the CABS MOs
    real(realk), pointer :: Ccabs(:,:)
    !> MO coefficient matrix for the RI MOs
    real(realk), pointer :: Cri(:,:)

    ! ***********************************************************
    !   Allocating for C matrix
    ! ***********************************************************
    !> Fock Fkj 
    real(realk), pointer :: Fkj(:,:)
    !> Fock Fmn 
    real(realk), pointer :: Fmn(:,:)
    !> Fock Fab
    real(realk), pointer :: Fab(:,:)
    !> Fock Fpq
    real(realk), pointer :: Fpq(:,:)

    ! ***********************************************************
    !   Allocating integer space sizes
    ! ***********************************************************
    !> number of AO orbitals
    integer :: nbasis
    !> number of occupied MO orbitals in EOS 
    integer :: noccEOS, nvirtEOS, noccfull
    !> number of occupied MO orbitals in AOS 
    integer :: noccAOS
    !> number of virtual MO orbitals in AOS 
    integer :: nvirtAOS
    !> number of occupied + virtual MO orbitals in EOS 
    integer :: nocvAOStot  

    !> number of CABS AO orbitals
    integer :: ncabsAO
    !> number of CABS MO orbitals
    integer :: ncabsMO

    integer :: ix, iy, i, j, m, n, k, l, p, q, c, r, s, t, a, b

    real(realk) :: V1energy, V2energy, V3energy, V4energy, V5energy
    real(realk) :: X1energy, X2energy, X3energy, X4energy 
    real(realk) :: B1energy, B2energy, B3energy, B4energy
    real(realk) :: B5energy, B6energy, B7energy, B8energy, B9energy  
    real(realk) :: E_21, E_21C, E_22, E_23, E_F12
    real(realk) :: tmp, energy, tmp2
    real(realk) :: temp
    real(realk) :: MP2energy, CCSDenergy
    
    real(realk) :: ECCSD_E21
    real(realk), pointer :: ECCSD_Vijab(:)
    real(realk), pointer :: ECCSD_Vijia(:)
    real(realk), pointer :: ECCSD_Vijaj(:)

    real(realk), pointer :: Venergy(:)
    real(realk), pointer :: Xenergy(:)
    real(realk), pointer :: Benergy(:)

    !> Timings
    real(realk) :: tcpu,twall
    logical :: Master,Collaborate,DoBasis
    integer :: n1,n2,n3,n4,Tain1,Tain2,noccAOStot,offset
#ifdef VAR_MPI
    Master = infpar%lg_mynum .EQ. infpar%master
    Collaborate = infpar%lg_nodtot .GT. 1
#else
    Master = .TRUE.
    Collaborate = .FALSE.
#endif

    ! ***********************************************************
    !   Sanity Check if we do the pair calculation
    ! ***********************************************************
    if((present(Fragment1) .AND. (.NOT. present(Fragment2))) .OR. &
         & (present(Fragment2) .AND. (.NOT. present(Fragment1)))) then
       call lsquit("get_f12_fragment_energy: Missing optional arguments Fragment1 and Fragment2")
    endif

    dopair = .FALSE.
    if(present(Fragment1) .AND. present(Fragment2)) then
       dopair = .TRUE.
    endif

    IF(Collaborate.And.Master)THEN
#ifdef VAR_MPI
       !wake up slaves
       call ls_mpibcast(F12_INTEGRAL_CALCULATION,infpar%master,infpar%lg_comm)
       !Broadcast info to slaves
       DoBasis = .TRUE.
       IF(dopair)THEN
          call decmpi_bcast_f12_info(MyFragment, Taibj,Tai,case,dopair,Dobasis,&
               & Fragment1,Fragment2)
       ELSE
          call decmpi_bcast_f12_info(MyFragment, Taibj,Tai,case,dopair,Dobasis)
       ENDIF
#endif
    ENDIF

    IF(Master)THEN

    call LSTIMER('START',tcpu,twall,DECinfo%output)
    
    
    nbasis   = MyFragment%nbasis
    noccEOS  = MyFragment%noccEOS
    nvirtEOS = MyFragment%nvirtEOS
    noccfull = noccEOS
       
    noccAOS = MyFragment%noccAOS
    nvirtAOS = MyFragment%nvirtAOS
    nvirtAOS = MyFragment%nvirtAOS

    ! For frozen core: noccAOS in only valence. But
    ! noccAOStot is core+valence, and nocvAOStot is core+valence+virtual,
    ! both with and without frozen core.
    noccAOStot = MyFragment%nocctot 
    nocvAOStot = noccAOStot + nvirtAOS

    ! Offset: Used for frozen core
    if(DECinfo%frozencore) then
       offset = MyFragment%ncore
    else
       offset=0
    end if
    
    ncabsAO = size(MyFragment%Ccabs,1)    
    ncabsMO = size(MyFragment%Ccabs,2)
    
    ! ***********************************************************
    !   Printing Input variables 
    ! ***********************************************************
    if(DECinfo%F12debug) then
       print *, "-------------------------------------------------"
       print *, "     F12-integrals.F90                           "
       print *, "-------------------------------------------------"
       print *, "nbasis:    ", nbasis
       print *, "noccEOS:   ", noccEOS
       print *, "nvirtEOS: ", nvirtEOS
       print *, "-------------------------------------------------"
       print *, "noccAOS    ", noccAOS
       print *, "noccAOStot ", noccAOStot
       print *, "nocvAOStot    ", nocvAOStot
       print *, "nvirtAOS   ", nvirtAOS
       print *, "ncabsAO    ", ncabsAO
       print *, "ncabsMO    ", ncabsMO
    end if
    
    ! Creating a CoccEOS matrix 
    call mem_alloc(CoccEOS, MyFragment%nbasis, noccEOS)
    do i=1, MyFragment%noccEOS
       ix = MyFragment%idxo(i)
       CoccEOS(:,i) = MyFragment%Co(:,ix)
    end do

    ! Creating a CoccAOStot matrix (always core+valence, also for frozen core)
    call mem_alloc(CoccAOStot, MyFragment%nbasis, noccAOStot)
    ! Only for frozen core: Include core MOs explicitly
    ! (without frozen core: core MOs are included in MyFragment%Co already)
    do i=1,offset
       CoccAOStot(:,i) = MyFragment%CoreMO(:,i)
    end do
    do i=1, MyFragment%noccAOS
       CoccAOStot(:,i+offset) = MyFragment%Co(:,i)
    end do

    ! Creating a CvirtAOS matrix 
    call mem_alloc(CvirtAOS, MyFragment%nbasis, nvirtAOS)
    do i=1, nvirtAOS
       CvirtAOS(:,i) = MyFragment%Cv(:,i)
    end do

    ! Creating a CocvAOStot matrix 
    call mem_alloc(CocvAOStot, MyFragment%nbasis, nocvAOStot)
    do i=1,noccAOStot
       CocvAOStot(:,i) = CoccAOStot(:,i)
    end do
    do i=1,nvirtAOS
       CocvAOStot(:,i+noccAOStot) = MyFragment%Cv(:,i)
    end do

    ! Creating a Ccabs matrix 
    call mem_alloc(Ccabs, ncabsAO, ncabsMO)
    do i=1, ncabsMO
       Ccabs(:,i) = MyFragment%Ccabs(:,i)
    end do

    ! Creating a Cri matrix 
    call mem_alloc(Cri, ncabsAO, ncabsAO)
    do i=1, ncabsAO
       Cri(:,i) = MyFragment%Cri(:,i)
    end do

    call mem_alloc(Venergy,5)

    WRITE(DECinfo%output,*) "Memory statistics after allocation of Venergy:"  
    call stats_globalmem(DECinfo%output)

    call get_EV1(Venergy, Fragment1, Fragment2, MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri) 
    call LSTIMER('get_EV1_timing: ',tcpu,twall,DECinfo%output)
   
    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EV1:"  
    call stats_globalmem(DECinfo%output)

    call get_EV2(Venergy, Fragment1, Fragment2, MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri) 
    call LSTIMER('get_EV2_timing: ',tcpu,twall,DECinfo%output)

    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EV2:"  
    call stats_globalmem(DECinfo%output)
    
    call get_EV3(Venergy, Fragment1, Fragment2, MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri) 
    call LSTIMER('get_EV3_timing: ',tcpu,twall,DECinfo%output)
    
    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EV3:"  
    call stats_globalmem(DECinfo%output)

    call get_EV4(Venergy, Fragment1, Fragment2, MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Taibj) 
    call LSTIMER('get_EV4_timing: ',tcpu,twall,DECinfo%output)

    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EV4:"  
    call stats_globalmem(DECinfo%output)
        
    E_21 = 0.0E0_realk
    E_21 = Venergy(1) + Venergy(2) + Venergy(3) + Venergy(4) + Venergy(5)

    ! MP2F12 CCoupling
    if(DECinfo%F12Ccoupling) then
        call MP2F12_Ccoupling_energy(MyFragment,E_21C)
        E_21 = E_21 + E_21C 
    endif

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, ' E21 V term                             '
       print *, '----------------------------------------'
       print *, " E21_CC_term:  ", E_21C
       print *, " E21_V_term1:  ", Venergy(1)
       print *, " E21_V_term2:  ", Venergy(2)
       print *, " E21_V_term3:  ", Venergy(3)
       print *, " E21_V_term4:  ", Venergy(4)
       print *, " E21_V_term5:  ", Venergy(5)
       print *, '----------------------------------------'
       print *, " E21_Vsum:     ", E_21
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E21 V term                             '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " E21_CC_term:  ", E_21C
       write(DECinfo%output,*) " E21_V_term1:  ", Venergy(1)
       write(DECinfo%output,*) " E21_V_term2:  ", Venergy(2)
       write(DECinfo%output,*) " E21_V_term3:  ", Venergy(3)
       write(DECinfo%output,*) " E21_V_term4:  ", Venergy(4)
       write(DECinfo%output,*) " E21_V_term5:  ", Venergy(5)
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " E21_Vsum:     ", E_21
    end if

    call mem_dealloc(Venergy)

    ! ***********************************************************
    ! Creating the F matrix 
    ! ***********************************************************
    ! Creating a Fkj MO matrix occ EOS

    call mem_alloc(Fkj, noccAOS, noccEOS)
    Fkj = 0E0_realk 
    do j=1, noccEOS      
       iy = MyFragment%idxo(j)
       Fkj(:,j) = MyFragment%ppfock(:,iy)
    end do
    ! Note that Fkj contains only valence orbitals since F(core,valence)=0.

    WRITE(DECinfo%output,*) "Memory statistics after allocation of Fkj:"  
    call stats_globalmem(DECinfo%output)

    ! Creating a Fmn MO matrix occ AOS
    call mem_alloc(Fmn, noccAOS, noccAOS)
    Fmn = 0E0_realk 

    WRITE(DECinfo%output,*) "Memory statistics after allocation of Fmn:"  
    call stats_globalmem(DECinfo%output)

    !> Double Storage! This need to be changed, just for conceptual reasons
    Fmn = MyFragment%ppfock

    ! Creating a Fab MO matrix virt AOS
    call mem_alloc(Fab, nvirtAOS, nvirtAOS)
    Fab = 0E0_realk 

    WRITE(DECinfo%output,*) "Memory statistics after allocation of Fab:"  
    call stats_globalmem(DECinfo%output)

    !> Double storage! This need to be changed, just for conceptual reasons
    Fab = MyFragment%qqfock
    
    call mem_alloc(Xenergy,4)

    WRITE(DECinfo%output,*) "Memory statistics after allocation of Xenergy:"  
    call stats_globalmem(DECinfo%output)

    call get_EX1(Xenergy,Fkj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    call LSTIMER('get_EX1_timing: ',tcpu,twall,DECinfo%output)
    
    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EX1:"  
    call stats_globalmem(DECinfo%output)

    call get_EX2(Xenergy,Fkj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    call LSTIMER('get_EX2_timing: ',tcpu,twall,DECinfo%output)
    
    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EX2:"  
    call stats_globalmem(DECinfo%output)

    call get_EX3(Xenergy,Fkj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    call LSTIMER('get_EX3_timing: ',tcpu,twall,DECinfo%output)
    
    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EX3:"  
    call stats_globalmem(DECinfo%output)

    call get_EX4(Xenergy,Fkj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    call LSTIMER('get_EX4_timing: ',tcpu,twall,DECinfo%output)

    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EX4:"  
    call stats_globalmem(DECinfo%output)

    E_22 = 0.0E0_realk
    E_22 = Xenergy(1) + Xenergy(2) + Xenergy(3) + Xenergy(4) 

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, ' E_22 X term                            '
       print *, '----------------------------------------'
       print *, " E22_X_term1: ", Xenergy(1)
       print *, " E22_X_term2: ", Xenergy(2)
       print *, " E22_X_term3: ", Xenergy(3)
       print *, " E22_X_term4: ", Xenergy(4)
       print *, '----------------------------------------'
       print *, " E22_Xsum: ", E_22

       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E_22 X term                            '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " E22_X_term1: ", Xenergy(1)
       write(DECinfo%output,*) " E22_X_term2: ", Xenergy(2)
       write(DECinfo%output,*) " E22_X_term3: ", Xenergy(3)
       write(DECinfo%output,*) " E22_X_term4: ", Xenergy(4)
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " E22_Xsum: ", E_22
    end if

    call mem_dealloc(Xenergy)

    ! ***********************************************************
    !   Creating the B matrix 
    ! ***********************************************************

    call mem_alloc(Benergy,9)
    
    WRITE(DECinfo%output,*) "Memory statistics after allocation of Benergy:"  
    call stats_globalmem(DECinfo%output)

    call get_EB1(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)
    call LSTIMER('get_EB1_timing: ',tcpu,twall,DECinfo%output)
    
    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EB1:"  
    call stats_globalmem(DECinfo%output)

    call get_EB2(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)   
    call LSTIMER('get_EB2_timing: ',tcpu,twall,DECinfo%output)
    
    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EB2:"  
    call stats_globalmem(DECinfo%output)

    call get_EB3(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)   
    call LSTIMER('get_EB3_timing: ',tcpu,twall,DECinfo%output)

    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EB3:"  
    call stats_globalmem(DECinfo%output)
   
    call get_EB4(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)   
    call LSTIMER('get_EB4_timing: ',tcpu,twall,DECinfo%output)

    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EB4:"  
    call stats_globalmem(DECinfo%output)

    call get_EB5(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)   
    call LSTIMER('get_EB5_timing: ',tcpu,twall,DECinfo%output)

    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EB5:"  
    call stats_globalmem(DECinfo%output)

    call get_EB6(Benergy,Fmn,Fab,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
         & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri)   
    call LSTIMER('get_EB6_timing: ',tcpu,twall,DECinfo%output)

    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EB6:"  
    call stats_globalmem(DECinfo%output)

    call get_EB7(Benergy,Fmn,Fab,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,&
         & CvirtAOS,CocvAOStot,Ccabs,Cri)   
    call LSTIMER('get_EB7_timing: ',tcpu,twall,DECinfo%output)

    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EB7:"  
    call stats_globalmem(DECinfo%output)

    call get_EB8(Benergy,Fmn,Fab,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,&
         & CvirtAOS,CocvAOStot,Ccabs,Cri)   
    call LSTIMER('get_EB8_timing: ',tcpu,twall,DECinfo%output)
    
    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EB8:"  
    call stats_globalmem(DECinfo%output)

    call get_EB9(Benergy,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,CoccAOStot,&
         & CvirtAOS,CocvAOStot,Ccabs,Cri)   
    call LSTIMER('get_EB9_timing: ',tcpu,twall,DECinfo%output)

    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EB9:"  
    call stats_globalmem(DECinfo%output)

    B1energy = Benergy(1)
    B2energy = Benergy(2)  
    B3energy = Benergy(3)  
    B4energy = Benergy(4)  
    B5energy = Benergy(5)  
    B6energy = Benergy(6)  
    B7energy = Benergy(7)  
    B8energy = Benergy(8)  
    B9energy = Benergy(9)  

    call mem_dealloc(Benergy)

    E_23 = 0.0E0_realk
    E_23 = B1energy + B2energy + B3energy + B4energy + B5energy + B6energy + B7energy &
         & + B8energy + B9energy

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, ' E_22 B term                            '
       print *, '----------------------------------------'
       print *, " E23_B_term1: ", B1energy
       print *, " E23_B_term2: ", B2energy   
       print *, " E23_B_term3: ", B3energy   
       print *, " E23_B_term4: ", B4energy   
       print *, " E23_B_term5: ", B5energy   
       print *, " E23_B_term6: ", B6energy   
       print *, " E23_B_term7: ", B7energy   
       print *, " E23_B_term8: ", B8energy   
       print *, " E23_B_term9: ", B9energy   
       print *, '----------------------------------------'
       print *, " E23_B_sum: ", E_23

       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E_22 B term                            '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " E23_B_term1: ", B1energy
       write(DECinfo%output,*) " E23_B_term2: ", B2energy   
       write(DECinfo%output,*) " E23_B_term3: ", B3energy   
       write(DECinfo%output,*) " E23_B_term4: ", B4energy   
       write(DECinfo%output,*) " E23_B_term5: ", B5energy   
       write(DECinfo%output,*) " E23_B_term6: ", B6energy   
       write(DECinfo%output,*) " E23_B_term7: ", B7energy   
       write(DECinfo%output,*) " E23_B_term8: ", B8energy   
       write(DECinfo%output,*) " E23_B_term9: ", B9energy   
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " E23_B_sum: ", E_23
    end if

    E_F12 = 0.0E0_realk
    E_F12 = E_21 + E_22 + E_23

    !> MP2-energy from an MP2-calculation
    MP2energy = Myfragment%energies(FRAGMODEL_OCCMP2)

   ! print *, "MP2energy: ", MP2energy

    if(DECinfo%F12debug) then
       print *,   '----------------------------------------------------------------'
       print *,   '                   DEC-MP2-F12 CALCULATION                      '
       print *,   '----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: MP2 CORRELATION ENERGY (For CC) =  ', MP2energy
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E_21
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E_22
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E_23
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 CORRECTION TO ENERGY =         ', E_F12
       print *, '-------------------------------------------------------'
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: TOTAL CORRELATION ENERGY (For CC) =', MP2energy+E_F12
    end if

    write(DECinfo%output,'(1X,a,f20.10)') '----------------------------------------------------------------'
    write(DECinfo%output,'(1X,a,f20.10)') '                  WANGY DEC-MP2-F12 CALCULATION                 '
    write(DECinfo%output,'(1X,a,f20.10)') '----------------------------------------------------------------'
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: MP2 CORRELATION ENERGY (For CC) =  ', MP2energy
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E_21
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E_22
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E_23
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 CORRECTION TO ENERGY =         ', E_F12
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: MP2-F12 CORRELATION ENERGY (CC) =  ', MP2energy+E_F12

    !> Setting the MP2-F12 correction
    Myfragment%energies(FRAGMODEL_MP2f12) = E_F12
    
    ! Which model? CCSD 
    WhichCCmodel: select case(case) 

    case(MODEL_NONE) ! SKip calculation

       return

    case(MODEL_CCSD)

!!$    print *,  '----------------------------------------'
!!$    print *, '  Tai                                   '
!!$    print *, '----------------------------------------'
!!$    DO i=1, noccEOS
!!$       DO a=1, nvirtAOS
!!$          print *, "a i value: ", a,i,Tai(a,i)
!!$       ENDDO
!!$    ENDDO

       ! **********************************************
       !   CCSD Vijab
       ! **********************************************
       call mem_alloc(ECCSD_Vijab,5)

       WRITE(DECinfo%output,*) "Memory statistics after allocation of ECCSD_Vijab energy:"  
       call stats_globalmem(DECinfo%output)

       call ccsdf12_Vijab_EV1(ECCSD_Vijab, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Taibj)
       call LSTIMER('ccsdf12_Vijab_EV1_timing: ',tcpu,twall,DECinfo%output)

       WRITE(DECinfo%output,*) "Memory statistics after subroutine ccsdf12_Vijab_EV1:"  
       call stats_globalmem(DECinfo%output)

       call ccsdf12_Vijab_EV2(ECCSD_Vijab, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Taibj)
       call LSTIMER('ccsdf12_Vijab_EV2_timing: ',tcpu,twall,DECinfo%output)

       WRITE(DECinfo%output,*) "Memory statistics after subroutine ccsdf12_Vijab_EV2:"  
       call stats_globalmem(DECinfo%output)

       call ccsdf12_Vijab_EV3(ECCSD_Vijab, Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Taibj)
       call LSTIMER('ccsdf12_Vijab_EV3_timing: ',tcpu,twall,DECinfo%output)

       WRITE(DECinfo%output,*) "Memory statistics after subroutine ccsdf12_Vijab_EV3:"  
       call stats_globalmem(DECinfo%output)

       call get_EV4(ECCSD_Vijab, Fragment1, Fragment2, MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Taibj) 
       call LSTIMER('get_EV4_timing: ',tcpu,twall,DECinfo%output)

       WRITE(DECinfo%output,*) "Memory statistics after subroutine get_EV4:"  
       call stats_globalmem(DECinfo%output)

       print *, '----------------------------------------'
       print *, ' E_CCSD_Vijab energies                  '
       print *, '----------------------------------------'
       print *, " ECCSD_Vijab_term1: ", ECCSD_Vijab(1)
       print *, " ECCSD_Vijab_term2: ", ECCSD_Vijab(2)
       print *, " ECCSD_Vijab_term3: ", ECCSD_Vijab(3)
       print *, " ECCSD_Vijab_term4: ", ECCSD_Vijab(4)
       print *, " ECCSD_Vijab_term5: ", ECCSD_Vijab(5)
       print *, '----------------------------------------'
       print *, " sum: ", ECCSD_Vijab(1) + ECCSD_Vijab(2) + ECCSD_Vijab(3) + ECCSD_Vijab(4) + ECCSD_Vijab(5) 

       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E_CCSD_Vijab energies                  '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " ECCSD_Vijab_term1: ", ECCSD_Vijab(1)
       write(DECinfo%output,*) " ECCSD_Vijab_term2: ", ECCSD_Vijab(2)
       write(DECinfo%output,*) " ECCSD_Vijab_term3: ", ECCSD_Vijab(3)
       write(DECinfo%output,*) " ECCSD_Vijab_term4: ", ECCSD_Vijab(4)
       write(DECinfo%output,*) " ECCSD_Vijab_term5: ", ECCSD_Vijab(5)
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " sum: ", ECCSD_Vijab(1) + ECCSD_Vijab(2) + ECCSD_Vijab(3) + ECCSD_Vijab(4) + ECCSD_Vijab(5) 
       ECCSD_E21 = ECCSD_Vijab(1) + ECCSD_Vijab(2) + ECCSD_Vijab(3) + ECCSD_Vijab(4) + ECCSD_Vijab(5) 

       call mem_dealloc(ECCSD_Vijab)

       ! ***********************************************************
       !    Next term V_ij^ia (V_ijia) or (Viija)
       ! **********************************************************
       call mem_alloc(ECCSD_Vijia,4)
       WRITE(DECinfo%output,*) "Memory statistics after allocation of ECCSD_Vijia energy:"  
       call stats_globalmem(DECinfo%output)
       
       call ccsdf12_Vijia_EV1(ECCSD_Vijia,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)
       call LSTIMER('ccsdf12_Vijia_EV1_timing: ',tcpu,twall,DECinfo%output)
       
       WRITE(DECinfo%output,*) "Memory statistics after subroutine ccsdf12_Vijia_EV1:"  
       call stats_globalmem(DECinfo%output)

       call ccsdf12_Vijia_EV2(ECCSD_Vijia,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)
       call LSTIMER('ccsdf12_Vijia_EV2_timing: ',tcpu,twall,DECinfo%output)
   
       WRITE(DECinfo%output,*) "Memory statistics after subroutine ccsdf12_Vijia_EV2:"  
       call stats_globalmem(DECinfo%output)

       call ccsdf12_Vijia_EV3(ECCSD_Vijia,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)
       call LSTIMER('ccsdf12_Vijia_EV3_timing: ',tcpu,twall,DECinfo%output)
   
       WRITE(DECinfo%output,*) "Memory statistics after subroutine ccsdf12_Vijia_EV3:"  
       call stats_globalmem(DECinfo%output)

       print *, '----------------------------------------'
       print *, ' E_CCSD_Vijia energies                  '
       print *, '----------------------------------------'
       print *, " ECCSD_Vijia_term1: ", ECCSD_Vijia(1)
       print *, " ECCSD_Vijia_term2: ", ECCSD_Vijia(2)
       print *, " ECCSD_Vijia_term3: ", ECCSD_Vijia(3)
       print *, " ECCSD_Vijia_term4: ", ECCSD_Vijia(4)
       print *, '----------------------------------------'
       print *, " sum: ", ECCSD_Vijia(1) + ECCSD_Vijia(2) + ECCSD_Vijia(3) + ECCSD_Vijia(4) 

       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E_CCSD_Vijia energies                  '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " ECCSD_Vijia_term1: ", ECCSD_Vijia(1)
       write(DECinfo%output,*) " ECCSD_Vijia_term2: ", ECCSD_Vijia(2)
       write(DECinfo%output,*) " ECCSD_Vijia_term3: ", ECCSD_Vijia(3)
       write(DECinfo%output,*) " ECCSD_Vijia_term4: ", ECCSD_Vijia(4)
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " sum: ", ECCSD_Vijia(1) + ECCSD_Vijia(2) + ECCSD_Vijia(3) + ECCSD_Vijia(4) 
       ECCSD_E21 = ECCSD_E21 + ECCSD_Vijia(1) + ECCSD_Vijia(2) + ECCSD_Vijia(3) + ECCSD_Vijia(4) 

       call mem_dealloc(ECCSD_Vijia)
       
       ! ***********************************************************
       !    Last term V_ij^aj 
       ! **********************************************************
       call mem_alloc(ECCSD_Vijaj,4)
       WRITE(DECinfo%output,*) "Memory statistics after allocation of ECCSD_Vijaj energy:"  
       call stats_globalmem(DECinfo%output)
       
       call ccsdf12_Viajj_EV1(ECCSD_Vijaj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)
       call LSTIMER('ccsdf12_Viajj_EV1_timing: ',tcpu,twall,DECinfo%output)
  
       WRITE(DECinfo%output,*) "Memory statistics after subroutine ccsdf12_Viajj_EV1:"  
       call stats_globalmem(DECinfo%output)
       
       call ccsdf12_Viajj_EV2(ECCSD_Vijaj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)
       call LSTIMER('ccsdf12_Viajj_EV2_timing: ',tcpu,twall,DECinfo%output)
     
       WRITE(DECinfo%output,*) "Memory statistics after subroutine ccsdf12_Viajj_EV2:"  
       call stats_globalmem(DECinfo%output)
       
       call ccsdf12_Viajj_EV3(ECCSD_Vijaj,Fragment1,Fragment2,MyFragment,dopair,CoccEOS,&
            & CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai)
       call LSTIMER('ccsdf12_Viajj_EV3_timing: ',tcpu,twall,DECinfo%output)
  
       WRITE(DECinfo%output,*) "Memory statistics after subroutine ccsdf12_Viajj_EV3:"  
       call stats_globalmem(DECinfo%output)

       print *, '----------------------------------------'
       print *, ' E_CCSD_Vijaj energies                  '
       print *, '----------------------------------------'
       print *, " ECCSD_Vijaj_term1: ", ECCSD_Vijaj(1)
       print *, " ECCSD_Vijaj_term2: ", ECCSD_Vijaj(2)
       print *, " ECCSD_Vijaj_term3: ", ECCSD_Vijaj(3)
       print *, " ECCSD_Vijaj_term4: ", ECCSD_Vijaj(4)
       print *, '----------------------------------------'
       print *, " sum: ", ECCSD_Vijaj(1) + ECCSD_Vijaj(2) + ECCSD_Vijaj(3) + ECCSD_Vijaj(4) 

       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E_CCSD_Vijaj energies                  '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " ECCSD_Vijaj_term1: ", ECCSD_Vijaj(1)
       write(DECinfo%output,*) " ECCSD_Vijaj_term2: ", ECCSD_Vijaj(2)
       write(DECinfo%output,*) " ECCSD_Vijaj_term3: ", ECCSD_Vijaj(3)
       write(DECinfo%output,*) " ECCSD_Vijaj_term4: ", ECCSD_Vijaj(4)
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) " sum: ", ECCSD_Vijaj(1) + ECCSD_Vijaj(2) + ECCSD_Vijaj(3) + ECCSD_Vijaj(4) 
       ECCSD_E21 = ECCSD_E21 + ECCSD_Vijaj(1) + ECCSD_Vijaj(2) + ECCSD_Vijaj(3) + ECCSD_Vijaj(4) 

    call mem_dealloc(ECCSD_Vijaj)

    E_F12 = ECCSD_E21 + E_21+E_22+E_23
    
    call get_ccsd_energy(CCSDenergy,MyFragment,CoccEOS,CoccAOStot,CvirtAOS,CocvAOStot,Ccabs,Cri,Tai,Taibj)
    call LSTIMER('get_ccsd_energy_timings: ',tcpu,twall,DECinfo%output)
  
    WRITE(DECinfo%output,*) "Memory statistics after subroutine get_ccsd_energy:"  
    call stats_globalmem(DECinfo%output)
    
    if(DECinfo%F12debug) then
       print *,   '----------------------------------------------------------------'
       print *,   '                   DEC-CCSD-F12 CALCULATION                     '
       print *,   '----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: E21 MP2 CORRECTION TO ENERGY = ', E_21
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: E21 CCSD    CORRECTION TO ENERGY = ', ECCSD_E21
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: E22+E23 MP2 CORRECTION TO ENERGY = ', E_22+E_23
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: CCSD-F12 CORRECTION TO ENERGY =    ', E_F12
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: CCSD ENERGY =                      ', CCSDenergy
       print *, '----------------------------------------'
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: TOTAL CORRELATION ENERGY =         ', CCSDenergy+E_F12
    end if

    write(DECinfo%output,'(1X,a,f20.10)') '----------------------------------------------------------------'
    write(DECinfo%output,'(1X,a,f20.10)') '                 WANGY DEC-CCSD-F12 CALCULATION                 '
    write(DECinfo%output,'(1X,a,f20.10)') '----------------------------------------------------------------'
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: E21 MP2 CORRECTION TO ENERGY = ', E_21
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: E21 CCSD    CORRECTION TO ENERGY = ', ECCSD_E21
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: E22+E23 MP2 CORRECTION TO ENERGY = ', E_22+E_23
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: CCSD-F12 CORRECTION TO ENERGY =    ', E_F12
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: CCSD ENERGY =                      ', CCSDenergy
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: TOTAL CORRELATION ENERGY =         ', CCSDenergy+E_F12

   !> Setting the CCSD-F12 correction (Needs to be changed)
    Myfragment%energies(FRAGMODEL_CCSDf12) = E_F12

    end select WhichCCmodel

    ! ***********************************************************
    !    Free Memory
    ! ***********************************************************

    !> Need to be free to avoid memory leak for the type(matrix) CMO_RI in CABS.F90
    ! call free_cabs()

    !> F-term
    call mem_dealloc(Fkj)
    call mem_dealloc(Fmn)
    call mem_dealloc(Fab)

    !> Coeff
    call mem_dealloc(CoccEOS)
    call mem_dealloc(CoccAOStot)
    call mem_dealloc(CocvAOStot)
    call mem_dealloc(Ccabs)
    call mem_dealloc(Cri)
    call mem_dealloc(CvirtAOS)
    ENDIF !if master

  end subroutine get_f12_fragment_energy

  !> Brief: Gives the single and pair fragment energy for MP2F12
  !> Author: Yang M. Wang
  !> Date: April 2013
  subroutine get_f12_fragment_energy_slave()
    implicit none
    !> Atomic fragment to be determined (Single or Pair fragment)
    type(decfrag) :: MyFragment
    !> t2EOS amplitudes stored in the order T(a,i,b,j)
    real(realk), pointer :: Taibj(:,:,:,:) 
    !> t1EOS amplitudes stored in the order T(a,i)
    real(realk), pointer :: Tai(:,:) 
    !> Case MODEL
    integer :: case
    !> Fragment 1 in the pair fragment
    type(decfrag) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag) :: Fragment2
    !
    logical :: DoBasis,PairFrag
    integer :: n1,n2,n3,n4,Tain1,Tain2
#ifdef VAR_MPI
    !Recieve info from master and allocates MyFragment and amplitudes    
    call decmpi_bcast_f12_info(MyFragment, Taibj, Tai, case,PairFrag,DoBasis,&
         & Fragment1, Fragment2)

    IF(PairFrag)THEN
       call get_f12_fragment_energy(MyFragment, Taibj, Tai, case, Fragment1, Fragment2)
       call atomic_fragment_free(Fragment1)
       call atomic_fragment_free(Fragment2)
    ELSE
       call get_f12_fragment_energy(MyFragment, Taibj, Tai, case)
    ENDIF

    call atomic_fragment_free(MyFragment)
    call mem_dealloc(Taibj)
    IF(case.EQ.MODEL_CCSD)THEN
       call mem_dealloc(Tai)
    ENDIF
#endif
  end subroutine get_f12_fragment_energy_slave

  !> Brief: MP2-F12 correction for the single fragment of term for the energies related to E21
  !> Author: Yang M. Wang
  !> Data: April 2013
  subroutine get_mp2f12_sf_E21(ijkl, nocc, energy, scalar)
    implicit none

    real(realk),intent(in)  :: scalar
    real(realk),intent(out) :: energy
    real(realk),intent(in)  :: ijkl(nocc,nocc,nocc,nocc)
    integer,intent(in)      :: nocc
    !
    integer     :: i,j
    real(realk) :: tmp,tmp2

    tmp = 0.0E0_realk
    do i=1, nocc
       tmp = tmp + ijkl(i,i,i,i)
    enddo

    energy = -1.0E0_realk*tmp ! The valeev factor
    tmp = 0E0_realk           ! NB Important reset

    do j=1, nocc
       do i=j+1, nocc 
          tmp = tmp  + 5.0E0_realk*ijkl(i,j,i,j) - ijkl(i,j,j,i)
       enddo
    enddo

    energy = energy - 0.5E0_realk*tmp ! The valeev factor
    energy = energy*scalar
  end subroutine get_mp2f12_sf_E21

  !> Brief: MP2-F12 correction for the single fragment of term for the energies related to E22
  !> Author: Yang M. Wang
  !> Data: August 2013
  subroutine get_mp2f12_sf_E22(Xijkl, n1, energy, scalar)
    implicit none

    real(realk),intent(in)  :: scalar
    real(realk),intent(out) :: energy
    integer,intent(in)      :: n1
    real(realk),intent(in)  :: Xijkl(n1,n1,n1,n1)
    !
    integer     :: i,j,k
    real(realk) :: tmp

    tmp = 0.0E0_realk
    energy = 0.0E0_realk

    do i=1, n1
       tmp = tmp + Xijkl(i,i,i,i)
    enddo

    energy = 0.25E0_realk*tmp

    tmp  = 0E0_realk         ! NB Important reset

    do j=1, n1
       do i=j+1, n1
          tmp = tmp +  7.0E0_realk * Xijkl(i,j,i,j) + Xijkl(i,j,j,i)
       enddo
    enddo

    energy = energy + 0.0625E0_realk*tmp
    energy = energy*scalar

  end subroutine get_mp2f12_sf_E22

  subroutine get_mp2f12_sf_E23(ijkl, nocc, energy, scalar)
    implicit none
    real(realk),intent(in)  :: scalar
    real(realk),intent(out) :: energy
    real(realk),intent(in)  :: ijkl(nocc,nocc,nocc,nocc)
    integer,intent(in)      :: nocc
    !
    integer     :: i,j
    real(realk) :: tmp
    
    tmp = 0.0E0_realk
    energy = 0.0E0_realk

    do i=1, nocc
       tmp = tmp + ijkl(i,i,i,i)
    enddo

    energy = 0.25E0_realk*tmp
    tmp = 0E0_realk         ! NB Important reset

    do j=1, nocc
       do i=j+1, nocc 
          tmp = tmp + 7.0E0_realk * ijkl(i,j,i,j) + ijkl(i,j,j,i)
       enddo
    enddo

    energy = energy + 0.0625E0_realk*tmp
    energy = energy*scalar

  end subroutine get_mp2f12_sf_E23

  !> Brief: MP2-F12 correction for the pair fragment of terms E21: V-terms 
  !> Author: Yang M. Wang
  !> Data: April 2013
  subroutine get_mp2f12_pf_E21(ijkl, Fragment1, Fragment2, PairFragment, nocc, energy, scalar)
    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(in) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(in) :: Fragment2
    !> Pair fragment formed from fragment 1 and 2
    type(decfrag), intent(inout) :: PairFragment
    !> The pair fragment energy
    real(realk),intent(out) :: energy
    !> Scalar to be multiplied with the energy
    real(realk),intent(in) :: scalar
    !> The integrals
    real(realk),intent(in) :: ijkl(nocc,nocc,nocc,nocc)
    !> Number of occupied MO orbitals in EOS space
    integer,intent(in) :: nocc
    !
    integer     :: i,j
    real(realk) :: tmp
    logical,pointer :: dopair_occ(:,:)

    call mem_alloc(dopair_occ,nocc,nocc)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
    
    tmp = 0.0E0_realk
    energy = 0.0E0_realk

    do j=1, nocc
       do i=1, nocc
          if(dopair_occ(i,j)) then !Do Pair 1 and 2 
             tmp = tmp + 5E0_realk * ijkl(i,j,i,j) - ijkl(i,j,j,i)
          endif
       enddo
    enddo
    energy = -0.25E0_realk*tmp
    energy = scalar*energy

    !Free memory
    call mem_dealloc(dopair_occ)

  end subroutine get_mp2f12_pf_E21

  !> Brief: MP2-F12 correction for the pair fragment of term for the energies related to E22
  !> Author: Yang M. Wang
  !> Data: August 2013
  subroutine get_mp2f12_pf_E22(Xijkl, n1, Fragment1, Fragment2, PairFragment, energy, scalar)
    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(in) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(in) :: Fragment2
    !> Pair fragment formed from fragment 1 and 2
    type(decfrag), intent(inout) :: PairFragment

    real(realk),intent(in)  :: scalar
    real(realk),intent(out) :: energy
    integer,intent(in)  :: n1
    real(realk),intent(in)  :: Xijkl(n1,n1,n1,n1)
    !
    integer     :: i,j,k
    real(realk) :: tmp,tmp2,tmp3
    logical,pointer :: dopair_occ(:,:)
   
    call mem_alloc(dopair_occ,n1,n1)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
    
    tmp = 0.0E0_realk
    energy = 0.0E0_realk

    do j=1, n1
       do i=j+1, n1 
          if(dopair_occ(i,j)) then !Do Pair 1 and 2   
             tmp = tmp +  7.0E0_realk * Xijkl(i,j,i,j) + Xijkl(i,j,j,i)
          endif
       enddo
    enddo

    energy = energy + 0.0625E0_realk*tmp
    energy = energy*scalar

    call mem_dealloc(dopair_occ)

  end subroutine get_mp2f12_pf_E22

  !> Brief: MP2-F12 correction for the pair fragment of term for the energies related to E23
  !> Author: Yang M. Wang
  !> Data: August 2013
  subroutine get_mp2f12_pf_E23(ijkl, Fragment1, Fragment2, PairFragment, nocc, energy, scalar)
    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(in) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(in) :: Fragment2
    !> Pair fragment formed from fragment 1 and 2
    type(decfrag), intent(inout) :: PairFragment
    
    real(realk),intent(in)  :: scalar
    real(realk),intent(out) :: energy
    real(realk),intent(in)  :: ijkl(nocc,nocc,nocc,nocc)
    integer,intent(in)      :: nocc
    !
    integer     :: i,j
    real(realk) :: tmp
    logical,pointer :: dopair_occ(:,:)

    call mem_alloc(dopair_occ,nocc,nocc)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)

    tmp = 0.0E0_realk
    energy = 0.0E0_realk

    do j=1, nocc
       do i=j+1, nocc 
          if(dopair_occ(i,j)) then !Do Pair 1 and 2   
             tmp = tmp + 7.0E0_realk * ijkl(i,j,i,j) + ijkl(i,j,j,i)
          endif
       enddo
    enddo
    energy = energy + 0.0625E0_realk*tmp
    energy = energy*scalar

    call mem_dealloc(dopair_occ)

  end subroutine get_mp2f12_pf_E23

#else

  subroutine wangy_dummy_sub12()
    implicit none
  end subroutine wangy_dummy_sub12

#endif

end module f12_integrals_module

