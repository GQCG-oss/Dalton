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
  use Integralparameters!, only: MP2INAMP
  use IntegralInterfaceMod!, only: II_getBatchOrbitalInfo
  use IntegralInterfaceDEC!, only: II_precalc_DECScreenMat,&
  !       & II_getBatchOrbitalScreen, II_GET_DECPACKED4CENTER_J_ERI

  use ccintegrals!, only: get_full_AO_integrals,get_AO_hJ,get_AO_K,get_AO_Fock
  use f12_routines_module!, only: MO_transform_AOMatrix

  ! Patricks mat_transpose routine 
  use reorder_frontend_module!, only: mat_transpose(rows,column,pref1,A,pref2,AT)
       
  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
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

  use wangy_playground_module

  public :: f12_single_fragment_energy, f12_pair_fragment_energy!, energy

  private

contains
  !> Brief: Gives the single fragment energy for MP2F12
  !> Author: Yang M. Wang
  !> Date: April 2013
  subroutine f12_single_fragment_energy(MyFragment)
    implicit none

    !> Atomic fragment to be determined  (NOT pair fragment)
    type(ccatom), intent(inout) :: MyFragment
    !> MO coefficient matrix for the occupied EOS
    real(realk), pointer :: CoccEOS(:,:)
    !> MO coefficient matrix for the occupied + virtual EOS
    real(realk), pointer :: CocvAOS(:,:)
    !> MO coefficient matrix for the CABS
    real(realk), pointer :: Ccabs(:,:)

    !> F12 integrals for the V1_term <ij|gr|kl>
    real(realk), pointer :: V1ijkl(:,:,:,:) 
    !> F12 integrals for the V2_term sum_pq <ij|g|pq> * <pq|r|kl>
    real(realk), pointer :: V2ijkl(:,:,:,:)   
    !> F12 integrals for the V2_term <ij|g|pq>
    real(realk), pointer :: Gijpq(:,:,:,:)
    real(realk), pointer :: Gpqij(:,:,:,:)
    real(realk), pointer :: Gpjqi(:,:,:,:)
    !> F12 integrals for the V2_term <ij|r|pq>
    real(realk), pointer :: Rpqij(:,:,:,:)     
    real(realk), pointer :: Rpjqi(:,:,:,:)
    !> F12 integrals for the V3_term <ij|g|ma'>
    real(realk), pointer :: Gmjci(:,:,:,:)
    !> F12 integrals for the V3_term <ij|r|ma'>
    real(realk), pointer :: Rmjci(:,:,:,:)

    !> F12 Fock Fij 
    real(realk), pointer :: Fij(:,:)
    !> F12 integral X1ijkl
    real(realk), pointer :: X1ijkl(:,:,:,:)
    !> F12 integrals for the X2_term sum_pq <ij|g|pq> <pq|g|kl>
    real(realk), pointer :: X2ijkl(:,:,:,:) 

    !> number of AO orbitals
    integer :: number_basis, nbasis
    !> number of occupied MO orbitals in EOS 
    integer :: noccEOS, nunoccEOS, noccfull
    !> number of occupied + virtual MO orbitals in EOS 
    integer :: nocvAOS  
  
    !> number of CABS AO orbitals
    integer :: ncabsAO
    !> number of CABS MO orbitals
    integer :: ncabsMO

    integer :: ix, i, j, m, n, k, l, p, q
    real(realk) :: V1energy, V2energy, V3energy, V4energy
    real(realk) :: X1energy, X2energy, X3energy, X4energy
    real(realk) :: tmp

    nbasis   = MyFragment%number_basis
    noccEOS  = MyFragment%noccEOS
    nunoccEOS = MyFragment%nunoccEOS
    noccfull = noccEOS

    nocvAOS = MyFragment%noccAOS + MyFragment%nunoccAOS
    
    ncabsAO = size(MyFragment%cabsMOs,1)    
    ncabsMO = size(MyFragment%cabsMOs,2)
   
    ! ***********************************************************
    ! Allocating memory for V matrix
    ! ***********************************************************
    call mem_alloc(V1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    V1ijkl = 0.0E0_realk
    call mem_alloc(V2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)  
    V2ijkl = 0.0E0_realk
    
    call mem_alloc(Gijpq, noccEOS, noccEOS, nocvAOS, nocvAOS)
    call mem_alloc(Gpqij, nocvAOS, nocvAOS, noccEOS, noccEOS)
    call mem_alloc(Gpjqi, nocvAOS, noccEOS, nocvAOS, noccEOS)
    Gijpq = 0.0E0_realk

    call mem_alloc(Rpqij, nocvAOS, nocvAOS, noccEOS, noccEOS) 
    call mem_alloc(Rpjqi, nocvAOS, noccEOS, nocvAOS, noccEOS) 
    Rpqij = 0.0E0_realk

    call mem_alloc(Gmjci, noccEOS, noccEOS, ncabsMO, noccEOS)
    call mem_alloc(Rmjci, noccEOS, noccEOS, ncabsMO, noccEOS)
    Gmjci = 0.0E0_realk
    Rmjci = 0.0E0_realk

 
    ! ***********************************************************
    ! Allocating memory for X matrix
    ! ***********************************************************
    call mem_alloc(X1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(X2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)

    ! ***********************************************************
    ! Allocating memory for F matrix
    ! ***********************************************************
    call mem_alloc(Fij, noccEOS, noccEOS)
!!$    
    ! ***********************************************************
    ! Creating a CoccEOS matrix 
    ! ***********************************************************
    call mem_alloc(CoccEOS, MyFragment%number_basis, noccEOS)
    do i=1, MyFragment%noccEOS
       ix = MyFragment%idxo(i)
       CoccEOS(:,i) = MyFragment%ypo(:,ix)
    end do

    ! ***********************************************************
    ! Creating a CocvAOS matrix 
    ! ***********************************************************
    call mem_alloc(CocvAOS, MyFragment%number_basis, nocvAOS)
     do i=1, MyFragment%noccAOS
       CocvAOS(:,i) = MyFragment%ypo(:,i)
    end do

    do i=1, MyFragment%nunoccAOS
       CocvAOS(:,i+MyFragment%noccAOS) = MyFragment%ypv(:,i)
    end do

    ! ***********************************************************
    ! Creating a Ccabs matrix 
    ! ***********************************************************
    call mem_alloc(Ccabs, ncabsAO, ncabsMO)
     do i=1, ncabsMO
       Ccabs(:,i) = MyFragment%cabsMOs(:,i)
    end do
    
    print *, '(Norm of CABS):'
    print *, '----------------------------------------'
    print *, 'norm2(Ccabs):', norm2(Ccabs)

    ! ***********************************************************
    ! Creating the V matrix 
    ! ***********************************************************
    ! Get integrals <ij|gr|kl> stored as (i,j,k,l) 
    call get_mp2f12_MOmatrix_ijkl(MyFragment%MyLsitem%Setting,nbasis,noccEOS,nocvAOS,CoccEOS,CocvAOS,V1ijkl,'RRRRF')  
    ! Get integrals <ij|g|pq> stored as (p,j,q,i)
    call get_mp2f12_MOmatrix_ijpq(MyFragment%MyLsitem%Setting,nbasis,noccEOS,nocvAOS,CoccEOS,CocvAOS,Gpjqi,'RRRRG')
    ! Get integrals <ij|r|pq> stored as (p,j,q,i)
    call get_mp2f12_MOmatrix_ijpq(MyFragment%MyLsitem%Setting,nbasis,noccEOS,nocvAOS,CoccEOS,CocvAOS,Rpjqi,'RRRRC')
    ! Get integrals <ij|g|mc> stored as (m,j,c,i) where c = cabs (a')
    !call get_mp2f12_MOmatrix_ijmc(MyFragment%MyLsitem%Setting,nbasis, &
    !     & noccEOS,ncabsAO,ncabsMO,CoccEOS,Ccabs,Gmjci,'RCRRG')

    !call get_mp2f12_MOmatrix_ijmc(MyFragment%MyLsitem%Setting,nbasis, &
    !     & noccEOS,nocvAOS,ncabsAO,ncabsMO,CoccEOS,CocvAOS,Ccabs,Rmjci,'RCRRC')
    
    ! Reorder the matrix Gpjqi to Gijpq with (nocvAOS,noccEOS,nocvAOS,noccEOS) as the dimensions of Gpjqi  
    call array_reorder_4d(1.0E0_realk,Gpjqi,nocvAOS,noccEOS,nocvAOS,noccEOS,[4,2,1,3],0.0E0_realk,Gijpq)
    ! Reorder the matrix Rpjqi to Rpqij with (nocvAOS,noccEOS,nocvAOS,noccEOS) as the dimensions of Rpjqi  
    call array_reorder_4d(1.0E0_realk,Rpjqi,nocvAOS,noccEOS,nocvAOS,noccEOS,[1,3,4,2],0.0E0_realk,Rpqij)
    ! Matrix multiplication of sum_pq <ij|g|pq><pq|r|kl>
    m = noccEOS*noccEOS
    k = nocvAOS*nocvAOS
    n = noccEOS*noccEOS
    call dec_simple_dgemm(m,k,n,Gijpq,Rpqij,V2ijkl, 'n', 'n')

    ! ***********************************************************
    ! Creating the F matrix 
    ! ***********************************************************
    ! Get integrals <i|F|j> stored as (i,j) !Fij
    Fij = Myfragment%ppfock
    
    ! ***********************************************************
    ! Creating the X matrix 
    ! ***********************************************************
    call get_mp2f12_MOmatrix_ijkl(MyFragment%MyLsitem%Setting,nbasis,noccEOS,nocvAOS,CoccEOS,CocvAOS,X1ijkl,'RRRR2')  

    m = noccEOS*noccEOS
    k = nocvAOS*nocvAOS
    n = noccEOS*noccEOS

    call mat_transpose(m,k,1.0E0_realk,Gijpq,0.0E0_realk,Gpqij)
    call dec_simple_dgemm(m,k,n,Gijpq,Gpqij,X2ijkl, 'n', 'n')

    ! ***********************************************************
    ! Printout statements 
    ! ***********************************************************
    print *, '----------------------------------------'
    print *, '(E21 Terms):'
    print *, '----------------------------------------'
    print *, 'norm2(V1ijkl):', norm2(V1ijkl)
    print *, '----------------------------------------'   
    print *, 'norm2(Rpqij):', norm2(Rpqij)
    print *, 'norm2(Gijpq):', norm2(Gijpq)
    print *, '----------------------------------------'
    print *, 'norm2(Gmjci):', norm2(Gmjci)
    print *, 'norm2(Rmjci):', norm2(Rmjci)
    print *, '----------------------------------------'
    print *, 'norm2(V1ijkl):', norm2(V1ijkl)
    print *, 'norm2(V2ijkl):', norm2(V2ijkl)
    print *, '----------------------------------------'

    V1energy = 0.0E0_realk
    V2energy = 0.0E0_realk
    X1energy = 0.0E0_realk
    X2energy = 0.0E0_realk

    print *, '(Single Fragment Energies):'
    print *, '----------------------------------------'
    call get_mp2f12_sf_E21(V1ijkl, noccEOS, V1energy, 1.0E0_realk)
    call get_mp2f12_sf_E21(V2ijkl, noccEOS, V2energy, -1.0E0_realk)
    print *, "E_21_V_term1:", V1energy
    print *, "E_21_V_term2:", V2energy
    print *, '----------------------------------------'
    print *, "E_21_Vsum:", V1energy + V2energy
    print *, '----------------------------------------' 

    print *, '(E22 Terms):'
    print *, '----------------------------------------'
    print *, 'norm2(Fij):'  , norm2(Fij)
    print *, 'norm2(Xijkl):', norm2(X1ijkl)
    print *, 'norm2(X2ijkl):', norm2(X2ijkl)
    print *, 'norm2(Gijpq):', norm2(Gijpq)
    print *, 'norm2(Gpqij):', norm2(Gpqij)
    print *, '----------------------------------------'
    print *, '(Single Fragment Energies):'
    print *, '----------------------------------------'
    call get_mp2f12_sf_E22(Fij, X1ijkl, noccEOS, X1energy, 1.0E0_realk)
    call get_mp2f12_sf_E22(Fij, X2ijkl, noccEOS, X2energy, -1.0E0_realk)
    print *, "E_22_X_term1:", X1energy
    print *, "E_22_X_term2:", X2energy
    print *, '----------------------------------------'
    print *, "E_22_Xsum:", X1energy + X2energy 
    print *, '----------------------------------------'
   
    myfragment%energies(14) = V1energy + V2energy + X1energy + X2energy 
    
    ! Free memory
    call mem_dealloc(CoccEOS)
    call mem_dealloc(CocvAOS)
    call mem_dealloc(Ccabs)

    call mem_dealloc(V1ijkl)

    call mem_dealloc(Gijpq)
    call mem_dealloc(Gpqij)
    call mem_dealloc(Gpjqi)

    call mem_dealloc(Rpqij)
    call mem_dealloc(Rpjqi)

    call mem_dealloc(Gmjci)
    call mem_dealloc(Rmjci)
    
    call mem_dealloc(V2ijkl)

    call mem_dealloc(Fij)
    call mem_dealloc(X1ijkl)
    call mem_dealloc(X2ijkl)

  end subroutine f12_single_fragment_energy

  !> Brief: Gives the pair fragment energy for MP2F12
  !> Author: Yang M. Wang
  !> Date: April 2013
  subroutine f12_pair_fragment_energy(Fragment1, Fragment2, PairFragment, natoms)

    implicit none
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Fragment 1 in the pair fragment
    type(ccatom),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(ccatom),intent(inout) :: Fragment2
    !> Pair fragment formed from fragment 1 and 2
    type(ccatom), intent(inout) :: PairFragment
    !> MO coefficient matrix for the occupied EOS
    real(realk), pointer :: CoccEOS(:,:)
    !> MO coefficient matrix for the occupied + virtual EOS
    real(realk), pointer :: CocvAOS(:,:)
    !> MO coefficient matrix for the CABS
    real(realk), pointer :: Ccabs(:,:)

    !> F12 integrals for the V1_term <ij|gr|kl>
    real(realk), pointer :: Fijkl(:,:,:,:) 
    !> F12 integrals for the V2_term <ij|g|pq>
    real(realk), pointer :: Gijpq(:,:,:,:)
    real(realk), pointer :: Gpjqi(:,:,:,:)
    !> F12 integrals for the V2_term <ij|r|pq>
    real(realk), pointer :: Rpqij(:,:,:,:)     
    real(realk), pointer :: Rpjqi(:,:,:,:)
    
    !> F12 integrals for the V2_term sum_pq <ij|g|pq> <pq|r|kl>
    real(realk), pointer :: V2ijkl(:,:,:,:) 

    !> number of AO orbitals
    integer :: number_basis, nbasis
    !> number of occupied MO orbitals in EOS 
    integer :: noccEOS 
    !> number of occupied + virtual MO orbitals in EOS 
    integer :: nocvAOS  

    !> number of CABS AO orbitals
    integer :: ncabsAO
    !> number of CABS MO orbitals
    integer :: ncabsMO

    integer :: ix, i, j, m, n, k, l, p, q
    real(realk) :: V1energy, V2energy
    real(realk) :: tmp
    logical,pointer :: dopair_occ(:,:)

    nbasis  = PairFragment%number_basis
    noccEOS = PairFragment%noccEOS
    nocvAOS = PairFragment%noccAOS + PairFragment%nunoccAOS

    call mem_alloc(CoccEOS, Pairfragment%number_basis, noccEOS)
    ! ***********************************************************
    ! Creating a CoccEOS matrix 
    ! ***********************************************************
    do i=1, PairFragment%noccEOS
       ix = PairFragment%idxo(i)
       CoccEOS(:,i) = PairFragment%ypo(:,ix)
    end do

    call mem_alloc(CocvAOS, PairFragment%number_basis, nocvAOS)
    ! ***********************************************************
    ! Creating a CocvAOS matrix 
    ! ***********************************************************
    do i=1, PairFragment%noccAOS
       CocvAOS(:,i) = PairFragment%ypo(:,i)
    end do

    do i=1, PairFragment%nunoccAOS
       CocvAOS(:,i+PairFragment%noccAOS) = PairFragment%ypv(:,i)
    end do
    
    ! ***********************************************************
    ! Creating a Ccabs matrix 
    ! ***********************************************************
    ncabsMO = size(PairFragment%cabsMOs,2)
    ncabsAO = size(PairFragment%cabsMOs,1)
    call mem_alloc(Ccabs, ncabsAO, ncabsMO)
    do i=1, ncabsMO
       Ccabs(:,i) = PairFragment%cabsMOs(:,i)
    end do
    
    call mem_alloc(dopair_occ,noccEOS,noccEOS)
    call mem_alloc(Fijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    Fijkl = 0.0E0_realk    
    call mem_alloc(V2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)  
   
    call mem_alloc(Gijpq, noccEOS, noccEOS, nocvAOS, nocvAOS)
    call mem_alloc(Gpjqi, nocvAOS, noccEOS, nocvAOS, noccEOS)
    Gijpq = 0.0E0_realk

    call mem_alloc(Rpqij, nocvAOS, nocvAOS, noccEOS, noccEOS) 
    call mem_alloc(Rpjqi, nocvAOS, noccEOS, nocvAOS, noccEOS) 
    Rpqij = 0.0E0_realk

    ! Get integrals <ij|gr|kl> stored as  (i,j,k,l)
    call get_mp2f12_MOmatrix_ijkl(PairFragment%MyLsitem%Setting,nbasis,noccEOS,nocvAOS,CoccEOS,CocvAOS,Fijkl,'RRRRF')  
    ! Get integrals <ij|g|pq> stored as  (p,j,p,i)  
    call get_mp2f12_MOmatrix_ijpq(PairFragment%MyLsitem%Setting,nbasis,noccEOS,nocvAOS,CoccEOS,CocvAOS,Gpjqi,'RRRRG')          
    ! Get integrals <ij|r|pq> stored as  (p,j,p,i)  
    call get_mp2f12_MOmatrix_ijpq(PairFragment%MyLsitem%Setting,nbasis,noccEOS,nocvAOS,CoccEOS,CocvAOS,Rpjqi,'RRRRC')   
    ! Reorder the matrix Gpjqi to Gijpq with (nocvAOS,noccEOS,nocvAOS,noccEOS) as the dimensions of Gpjqi  
    call array_reorder_4d(1.0E0_realk,Gpjqi,nocvAOS,noccEOS,nocvAOS,noccEOS,[4,2,1,3],0.0E0_realk,Gijpq)
    ! Reorder the matrix Rpjqi to Rpqij with (nocvAOS,noccEOS,nocvAOS,noccEOS) as the dimensions of Rpjqi  
    call array_reorder_4d(1.0E0_realk,Rpjqi,nocvAOS,noccEOS,nocvAOS,noccEOS,[1,3,4,2],0.0E0_realk,Rpqij)
    ! Matrix multiplication of sum_pq <ij|g|pq><pq|r|kl>
    m = noccEOS*noccEOS
    k = nocvAOS*nocvAOS
    n = noccEOS*noccEOS
    call dec_simple_dgemm(m,k,n,Gijpq,Rpqij,V2ijkl, 'n', 'n')

    V1energy = 0.0E0_realk
    V2energy = 0.0E0_realk
    
    print *, '(Inside Pair Fragment):'
    print *, '----------------------------------------'
    print *, "norm2(Fijkl):", norm2(Fijkl) 
    print *, '----------------------------------------'
    
    print *, '(Pair Fragment Energies):'
    call get_mp2f12_pf_E21(Fijkl,  Fragment1, Fragment2, PairFragment, noccEOS, V1energy, 1.0E0_realk )
    print *, "E_21_V_term1:", V1energy
    call get_mp2f12_pf_E21(V2ijkl, Fragment1, Fragment2, PairFragment, noccEOS, V2energy, -1.0E0_realk )
    print *, "E_21_V_term2:", V2energy
    print *, '----------------------------------------'
    print *, "E_21_Vsum:", V1energy + V2energy
    print *, '----------------------------------------'
    
    ! Input for Dec Driver
    pairfragment%energies(14) = V1energy + V2energy

    ! Free memory
    call mem_dealloc(dopair_occ)
  
    ! Free memory
    call mem_dealloc(CoccEOS)
    call mem_dealloc(CocvAOS)
    call mem_dealloc(Ccabs)

    call mem_dealloc(Fijkl)

    call mem_dealloc(Gijpq)
    call mem_dealloc(Gpjqi)

    call mem_dealloc(Rpqij)
    call mem_dealloc(Rpjqi)

    call mem_dealloc(V2ijkl)

  end subroutine f12_pair_fragment_energy

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

    tmp = 0E0_realk
    do i=1, nocc
       tmp = tmp + ijkl(i,i,i,i)
    enddo

    energy = -1E0_realk*tmp ! The valeev factor
    tmp = 0E0_realk         ! NB Important reset

    do i=1, nocc
       do j=i+1, nocc 
          tmp = tmp  + 5E0_realk*ijkl(i,j,i,j) - 1E0_realk*ijkl(i,j,j,i)
       enddo
    enddo
    energy = energy - 0.5E0_realk*tmp ! The valeev factor
    energy = energy*scalar

  end subroutine get_mp2f12_sf_E21

!> Brief: MP2-F12 correction for the single fragment of term for the energies related to E22
  !> Author: Yang M. Wang
  !> Data: August 2013
  subroutine get_mp2f12_sf_E22(Fij, Xijkl, nocc, energy, scalar)
    implicit none

    real(realk),intent(in)  :: scalar
    real(realk),intent(out) :: energy
    real(realk),intent(in)  :: Fij(nocc,nocc)
    real(realk),intent(in)  :: Xijkl(nocc,nocc,nocc,nocc)
    integer,intent(in)      :: nocc
    !
    integer     :: i,j
    real(realk) :: tmp,tmp2

    real(realk), pointer :: Bijkl(:,:,:,:)
    tmp = 0E0_realk

    call mem_alloc(Bijkl,nocc,nocc,nocc,nocc)

    do j=1,nocc
       do i=1,nocc
          tmp2 = Fij(i,i) + Fij(j,j)
          Bijkl(i,j,i,j) = -1.0E0_realk*tmp2*Xijkl(i,j,i,j)
          Bijkl(i,j,j,i) = -1.0E0_realk*tmp2*Xijkl(i,j,j,i)
       enddo
    enddo
    
    do i=1, nocc
       tmp = tmp + Bijkl(i,i,i,i)
    enddo

    energy = 0.25E0_realk*tmp
    tmp = 0E0_realk         ! NB Important reset

    do i=1, nocc
       do j=i+1, nocc 
          tmp = tmp +  7.0E0_realk * Bijkl(i,j,i,j) + Bijkl(i,j,j,i)
       enddo
    enddo
    energy = energy + 0.0625E0_realk*tmp
    energy = energy*scalar

    call mem_dealloc(Bijkl)

  end subroutine get_mp2f12_sf_E22

  !> Brief: MP2-F12 correction for the pair fragment of term V1: E_PQ(V1) 
  !> Author: Yang M. Wang
  !> Data: April 2013
  subroutine get_mp2f12_pf_E21(ijkl, Fragment1, Fragment2, PairFragment, nocc, energy, scalar)
    implicit none
    !> Fragment 1 in the pair fragment
    type(ccatom),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(ccatom),intent(inout) :: Fragment2
    !> Pair fragment formed from fragment 1 and 2
    type(ccatom), intent(inout) :: PairFragment
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
    tmp = 0.0_realk

    call mem_alloc(dopair_occ,nocc,nocc)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)

    do i=1,nocc
       do j=1,nocc
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

  !> Brief: Integral print
  !> Author: Yang M. Wang
  !> Data: August 2013
  subroutine matrix_print_4d(A, p, q, r, s)
    implicit none

    real(realk),intent(in)  :: A(p,q,r,s)
    integer,intent(in)      :: p,q,r,s
    !
    integer :: i,j,k,l
 
    do i=1, p
       do j=1, q
          do k=1, r
             do l=1, s 
                if(abs(A(i,j,k,l)) > 1E-10_realk) then
                   print *, i,j,k,l, A(i,j,k,l)
                else
                   print *, i,j,k,l, 0E0_realk
                endif
             enddo
          enddo
       enddo
    enddo
    
  end subroutine matrix_print_4d
  
  
  !> Brief: Get <ij|OPERATOR|kl> integrals stored in the order (k,l,m,n).
  !> Author: Yang M. Wang
  !> Data: June 2013
  subroutine get_mp2f12_MOmatrix_ijkl(MySetting,nbasis,noccEOS,nocvAOS,CoccEOS,CocvAOS,ijkl,INTSPEC) 
    implicit none
    
    !> Integrals settings
    type(lssetting), intent(inout) :: mysetting
    !> Number of basis functions AO
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals MO in EOS space
    integer,intent(in) :: noccEOS
    !> Number of occupied and virtual MO in AOS space 
    integer,intent(in) :: nocvAOS
    !> Occupied MO coefficients
    real(realk),intent(in),dimension(nbasis,noccEOS) :: CoccEOS
    !> Occupied and virtual MO coefficients in AOS
    real(realk),intent(in),dimension(nbasis,nocvAOS) :: CocvAOS
    !>  <ij |gr| ij> integrals stored in the order (i,j,j,i)
    real(realk), intent(inout) :: ijkl(noccEOS,noccEOS,noccEOS,noccEOS)   
    real(realk), pointer :: kjli(:,:,:,:)

    integer :: alphaB,gammaB,dimAlpha,dimGamma,GammaStart, GammaEnd, AlphaStart, AlphaEnd
    real(realk),pointer :: tmp1(:),tmp2(:),CoccEOST(:,:),CocvAOST(:,:)
    integer(kind=long) :: dim1,dim2
    integer :: m,k,n,idx
    logical :: FullRHS,doscreen
    integer :: MaxActualDimAlpha,nbatchesAlpha,nbatches,MaxActualDimGamma,nbatchesGamma,iorb
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
    type(batchtoorb), pointer :: batch2orbAlpha(:),batch2orbGamma(:)
#ifdef VAR_OMP
    integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
    TYPE(DECscreenITEM)   :: DecScreen
    Character            :: intSpec(5)
    integer :: MinAObatchSize, MaxAObatchSize, GammaBatchSize, AlphaBatchSize

    ! ***********************************************************
    ! For efficiency when calling dgemm, save transposed matrices
    ! ***********************************************************
    call mem_alloc(CoccEOST,noccEOS,nbasis)
    call mat_transpose(nbasis,noccEOS, 1.0E0_realk, CoccEOS, 0.0E0_realk, CoccEOST)

    ! ************************************************************
    ! Allocate mem space for a temporary V1 that will be reordered
    ! ************************************************************
    call mem_alloc(kjli,noccEOS,noccEOS,noccEOS,noccEOS)

    ! ************************
    ! Determine AO batch sizes
    ! ************************
    ! NOTE: Ideally the batch sizes should be optimized according to the available memory
    ! (as is done e.g. in get_optimal_batch_sizes_for_mp2_integrals).
    ! For simplicity we simply choose the gamma batch to contain all basis functions,
    ! while we make the alpha batch as small as possible

    ! Minimum AO batch size
    call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,MinAObatchSize,'R')

    ! Maximum AO batch size (all basis functions)
    MaxAObatchSize = nbasis
    ! Setting MinAO to AO batch size for debug purposes
    MinAObatchSize = nbasis

    ! Set alpha and gamma batch size as written above
    GammaBatchSize = MaxAObatchSize
    AlphaBatchSize = MinAObatchSize
    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
    ! Orbital to batch information
    ! ----------------------------

    call mem_alloc(orb2batchGamma,nbasis)
    call build_batchesofAOS(DECinfo%output,mysetting,GammaBatchSize,nbasis,MaxActualDimGamma,&
         & batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,orb2BatchGamma,'R')
    ! Batch to orbital information
    ! ----------------------------
    call mem_alloc(batch2orbGamma,nbatchesGamma)
    do idx=1,nbatchesGamma
       call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
       batch2orbGamma(idx)%orbindex = 0
       batch2orbGamma(idx)%norbindex = 0
    end do
    do iorb=1,nbasis
       idx = orb2batchGamma(iorb)
       batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
       K = batch2orbGamma(idx)%norbindex
       batch2orbGamma(idx)%orbindex(K) = iorb
    end do

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************
    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchAlpha,nbasis)
    call build_batchesofAOS(DECinfo%output,mysetting,AlphaBatchSize,nbasis,&
         & MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,'R')

    ! Batch to orbital information
    ! ----------------------------
    call mem_alloc(batch2orbAlpha,nbatchesAlpha)
    do idx=1,nbatchesAlpha
       call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
       batch2orbAlpha(idx)%orbindex = 0
       batch2orbAlpha(idx)%norbindex = 0
    end do
    do iorb=1,nbasis
       idx = orb2batchAlpha(iorb)
       batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
       K = batch2orbAlpha(idx)%norbindex
       batch2orbAlpha(idx)%orbindex(K) = iorb
    end do

    ! Integral screening stuff
    doscreen = Mysetting%scheme%cs_screen .or. Mysetting%scheme%ps_screen
    call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mysetting,&
         & nbatches,nbatchesAlpha,nbatchesGamma,INTSPEC)
    IF(doscreen)then
       call II_getBatchOrbitalScreen(DecScreen,mysetting,&
            & nbasis,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,DECinfo%output,DECinfo%output)
    endif
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
       dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
       GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
       GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
 
       BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
          dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
          AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
          AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch

          ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
          ! ************************************************************************************
          dim1 = i8*nbasis*nbasis*dimAlpha*dimGamma   ! dimension for integral array

          call mem_alloc(tmp1,dim1)
          ! Store integral in tmp1(1:dim1) array in (beta,delta,alphaB,gammaB) order
          IF(doscreen) mysetting%LST_GAB_RHS => DECSCREEN%masterGabRHS
          IF(doscreen) mysetting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p
          call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
               & mysetting, tmp1, batchindexAlpha(alphaB), batchindexGamma(gammaB), &
               & batchsizeAlpha(alphaB), batchsizeGamma(gammaB), nbasis, nbasis, dimAlpha, dimGamma,FullRHS,&
               & nbatches,INTSPEC)

          !print *, 'norm2(tmp1)', norm2(tmp1)

          ! Transform beta to occupied index "j".
          ! *************************************
          ! Note: ";" indicates the place where the array is transposed:
          ! tmp2(delta,alphaB,gammaB,j) = sum_{beta} tmp1^T(beta;delta,alphaB,gammaB) Cocc_{beta j}
          m = nbasis*dimGamma*dimAlpha   ! # elements in "delta alphaB gammaB" dimension of tmp1^T
          k = nbasis                     ! # elements in "beta" dimension of tmp1^T
          n = noccEOS                    ! # elements in second dimension of Cocc
          dim2 = i8*noccEOS*nbasis*dimAlpha*dimGamma  ! dimension of tmp2 array

          call mem_alloc(tmp2,dim2)
          call dec_simple_dgemm(m,k,n,tmp1,CoccEOS,tmp2, 't', 'n')
          call mem_dealloc(tmp1)

          !print *, 'norm2(tmp2)', norm2(tmp2)
          
          ! Transform beta to unoccupied index "b".
          ! ***************************************
          ! tmp3(b,alphaB,gammaB,j) = sum_{delta} CunoccT(b,delta) tmp2(delta,alphaB,gammaB,j)
          ! Note: We have stored the transposed Cunocc matrix, so no need to transpose in
          ! the call to dgemm.
          
          m = noccEOS
          k = nbasis
          n = dimAlpha*dimGamma*noccEOS
          dim1 = i8*noccEOS*noccEOS*dimAlpha*dimGamma  ! dimension of tmp2 array
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,CoccEOST,tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2) 

          !print *, 'norm2(CoccEOST)', norm2(CoccEOST)
          !print *, 'norm2(tmp3)', norm2(tmp1)

          ! Transpose to make alphaB and gammaB indices available
          ! *****************************************************
          dim2=dim1
          call mem_alloc(tmp2,dim2)
          ! tmp2(gammaB, j, b, alphaB) = tmp1^T(b, alphaB; gammaB, j)
          m = noccEOS*dimAlpha      ! dimension of "row" in tmp1 array (to be "column" in tmp2
          n = noccEOS*dimGamma      ! dimension of "column" in tmp1 array (to be "row" in tmp2)

          call mat_transpose(m, n, 1.0E0_realk, tmp1, 0.0E0_realk,tmp2)
          call mem_dealloc(tmp1)

          !print *, 'norm2(tmp4) (transposed)', norm2(tmp2)

          ! Transform gamma batch index to occupied index
          ! *********************************************
          ! tmp4(i,j,b,alphaB) = sum_{gamma in gammaBatch} CoccT(i,gamma) tmp2(gamma,j,b,alphaB)
          m = noccEOS
          k = dimGamma
          n = noccEOS*noccEOS*dimAlpha
          dim1 = i8*noccEOS*noccEOS*noccEOS*dimAlpha
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,CoccEOST(:,GammaStart:GammaEnd),tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2)
        
          !print *, 'norm2(tmp5)', norm2(tmp1)

          ! Transform alpha batch index to unoccupied index and update output integral
          ! **************************************************************************
          ! ijba(i,j,b,a) =+ sum_{alpha in alphaBatch} tmp1(i,j,b,alpha)  CoccEOS(alpha,a)
          m = noccEOS*noccEOS*noccEOS
          k = dimAlpha
          n = noccEOS
          dim2 = i8*m*n
          call dec_simple_dgemm_update(m,k,n,tmp1,CoccEOST(:,AlphaStart:AlphaEnd),kjli, 'n', 't')
          call mem_dealloc(tmp1)

          !print *, 'norm2(kjli)', norm2(kjli)

          ! Note: To have things consecutive in memory it is better to pass CunoccT to the dgemm
          ! routine and then transpose (insted of passing CoccEOS and not transpose).

       end do BatchAlpha
    end do BatchGamma

    !> Reorder tmp(k,j,l,i) -> tmp(i,j,k,l)
    call array_reorder_4d(1.0E0_realk,kjli,noccEOS,noccEOS,noccEOS,noccEOS,[4,2,1,3],0.0E0_realk,ijkl)

    ! Free and nullify stuff
    ! **********************
    nullify(mysetting%LST_GAB_LHS)
    nullify(mysetting%LST_GAB_RHS)
    call free_decscreen(DECSCREEN)

    ! Free gamma batch stuff
    call mem_dealloc(orb2batchGamma)
    call mem_dealloc(batchdimGamma)
    call mem_dealloc(batchsizeGamma)
    call mem_dealloc(batchindexGamma)
    orb2batchGamma => null()
    batchdimGamma => null()
    batchsizeGamma => null()
    batchindexGamma => null()
    do idx=1,nbatchesGamma
       call mem_dealloc(batch2orbGamma(idx)%orbindex)
       batch2orbGamma(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbGamma)
    batch2orbGamma => null()

    ! Free alpha batch stuff
    call mem_dealloc(orb2batchAlpha)
    call mem_dealloc(batchdimAlpha)
    call mem_dealloc(batchsizeAlpha)
    call mem_dealloc(batchindexAlpha)
    orb2batchAlpha => null()
    batchdimAlpha => null()
    batchsizeAlpha => null()
    batchindexAlpha => null()
    do idx=1,nbatchesAlpha
       call mem_dealloc(batch2orbAlpha(idx)%orbindex)
       batch2orbAlpha(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbAlpha)
    batch2orbAlpha => null()

    ! Free F12 related pointers
    call mem_dealloc(CoccEOST)
    !call mem_dealloc(CocvAOST)
    call mem_dealloc(kjli)
    
  end subroutine get_mp2f12_MOmatrix_ijkl


  !> Brief: Get <ij|OPERATOR|pq> integrals stored in the order (i,j,p,q).
  !> Author: Yang M. Wang
  !> Data: June 2013
  subroutine get_mp2f12_MOmatrix_ijpq(MySetting,nbasis,noccEOS,nocvAOS,CoccEOS,CocvAOS,pjqi,INTSPEC) 
    implicit none

    !> Integrals settings
    type(lssetting), intent(inout) :: mysetting
    !> Number of basis functions AO
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals MO in EOS
    integer,intent(in) :: noccEOS
    !> Number of occupied and virtual orbitals MO in AOS
    integer,intent(in) :: nocvAOS
    !> Occupied MO coefficients
    real(realk),intent(in),dimension(nbasis,noccEOS) :: CoccEOS
    !> Occupied and virtual MO coefficients in AOS
    real(realk),intent(in),dimension(nbasis,nocvAOS) :: CocvAOS
    !>  <ij |OPERATOR| pq> integrals stored in the order (i,j,p,q)
    real(realk), intent(inout) :: pjqi(nocvAOS,noccEOS,nocvAOS,noccEOS)   

    integer :: alphaB,gammaB,dimAlpha,dimGamma,GammaStart, GammaEnd, AlphaStart, AlphaEnd
    real(realk),pointer :: tmp1(:),tmp2(:),CoccEOST(:,:),CocvAOST(:,:)
    integer(kind=long) :: dim1,dim2
    integer :: m,k,n,idx
    logical :: FullRHS,doscreen
    integer :: MaxActualDimAlpha,nbatchesAlpha,nbatches,MaxActualDimGamma,nbatchesGamma,iorb
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
    type(batchtoorb), pointer :: batch2orbAlpha(:),batch2orbGamma(:)
#ifdef VAR_OMP
    integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
    TYPE(DECscreenITEM)   :: DecScreen
    Character            :: intSpec(5)
    integer :: MinAObatchSize, MaxAObatchSize, GammaBatchSize, AlphaBatchSize

    ! ***********************************************************
    ! For efficiency when calling dgemm, save transposed matrices
    ! ***********************************************************
    call mem_alloc(CoccEOST,noccEOS,nbasis)
    call mem_alloc(CocvAOST,nocvAOS,nbasis)

    call mat_transpose(nbasis,noccEOS,1E0_realk,CoccEOS,0E0_realk,CoccEOST)
    call mat_transpose(nbasis,nocvAOS,1E0_realk,CocvAOS,0E0_realk,CocvAOST)

    ! ************************************************************
    ! Allocate mem space for a temporary V1 that will be reordered
    ! ************************************************************
    ! call mem_alloc(pjqi,nocvAOS,noccEOS,nocvAOS,noccEOS)

    ! ************************
    ! Determine AO batch sizes
    ! ************************
    ! NOTE: Ideally the batch sizes should be optimized according to the available memory
    ! (as is done e.g. in get_optimal_batch_sizes_for_mp2_integrals).
    ! For simplicity we simply choose the gamma batch to contain all basis functions,
    ! while we make the alpha batch as small as possible

    ! Minimum AO batch size
    call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,MinAObatchSize,'R')

    ! Maximum AO batch size (all basis functions)
    MaxAObatchSize = nbasis
    ! Setting AO batch size (all basis functions) For testing!
    ! MinAObatchSize = nbasis

    ! Set alpha and gamma batch size as written above
    GammaBatchSize = MaxAObatchSize
    AlphaBatchSize = MinAObatchSize
    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
    ! Orbital to batch information
    ! ----------------------------

    call mem_alloc(orb2batchGamma,nbasis)
    call build_batchesofAOS(DECinfo%output,mysetting,GammaBatchSize,nbasis,MaxActualDimGamma,&
         & batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,orb2BatchGamma,'R')
    ! Batch to orbital information
    ! ----------------------------
    call mem_alloc(batch2orbGamma,nbatchesGamma)
    do idx=1,nbatchesGamma
       call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
       batch2orbGamma(idx)%orbindex = 0
       batch2orbGamma(idx)%norbindex = 0
    end do
    do iorb=1,nbasis
       idx = orb2batchGamma(iorb)
       batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
       K = batch2orbGamma(idx)%norbindex
       batch2orbGamma(idx)%orbindex(K) = iorb
    end do

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************
    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchAlpha,nbasis)
    call build_batchesofAOS(DECinfo%output,mysetting,AlphaBatchSize,nbasis,&
         & MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,'R')

    ! Batch to orbital information
    ! ----------------------------
    call mem_alloc(batch2orbAlpha,nbatchesAlpha)
    do idx=1,nbatchesAlpha
       call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
       batch2orbAlpha(idx)%orbindex = 0
       batch2orbAlpha(idx)%norbindex = 0
    end do
    do iorb=1,nbasis
       idx = orb2batchAlpha(iorb)
       batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
       K = batch2orbAlpha(idx)%norbindex
       batch2orbAlpha(idx)%orbindex(K) = iorb
    end do

    ! Integral screening stuff
    doscreen = Mysetting%scheme%cs_screen .or. Mysetting%scheme%ps_screen
    !doscreen = .FALSE.
    call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mysetting,&
         & nbatches,nbatchesAlpha,nbatchesGamma,INTSPEC)
    IF(doscreen)then
       call II_getBatchOrbitalScreen(DecScreen,mysetting,&
            & nbasis,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,DECinfo%output,DECinfo%output)
    endif
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
    pjqi = 0.0_realk

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
       dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
       GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
       GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
 
       BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
          dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
          AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
          AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch

          ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
          ! ************************************************************************************
          dim1 = i8*nbasis*nbasis*dimAlpha*dimGamma   ! dimension for integral array

          call mem_alloc(tmp1,dim1)
          ! Store integral in tmp1(1:dim1) array in (beta,delta,alphaB,gammaB) order
          IF(doscreen) mysetting%LST_GAB_RHS => DECSCREEN%masterGabRHS
          IF(doscreen) mysetting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p
          call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
               & mysetting, tmp1, batchindexAlpha(alphaB), batchindexGamma(gammaB), &
               & batchsizeAlpha(alphaB), batchsizeGamma(gammaB), nbasis, nbasis, dimAlpha, dimGamma,FullRHS,&
               & nbatches,INTSPEC)

          !print *,"norm2(tmp1):", norm2(tmp1)

          ! Transform beta to occupied index "j".
          ! *************************************
          ! Note: ";" indicates the place where the array is transposed:
          ! tmp2(delta,alphaB,gammaB,j) = sum_{beta} tmp1^T(beta;delta,alphaB,gammaB) Cocc_{beta j}
          m = nbasis*dimGamma*dimAlpha   ! # elements in "delta alphaB gammaB" dimension of tmp1^T
          k = nbasis                     ! # elements in "beta" dimension of tmp1^T
          n = noccEOS                    ! # elements in second dimension of Cocc
          dim2 = i8*m*n ! dimension of tmp2 array

          call mem_alloc(tmp2,dim2)
          call dec_simple_dgemm(m,k,n,tmp1,CoccEOS,tmp2, 't', 'n')
          call mem_dealloc(tmp1)

          !print *,"norm2(tmp2):", norm2(tmp2)

          ! Transform beta to unoccupied index "b".
          ! ***************************************
          ! tmp3(b,alphaB,gammaB,j) = sum_{delta} CunoccT(b,delta) tmp2(delta,alphaB,gammaB,j)
          ! Note: We have stored the transposed Cunocc matrix, so no need to transpose in
          ! the call to dgemm.
          
          m = nocvAOS
          k = nbasis
          n = dimAlpha*dimGamma*noccEOS
          dim1 = i8*m*n  ! dimension of tmp2 array
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,CocvAOST,tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2) 

          !print *,"norm2(tmp3):", norm2(tmp1)


          ! Transpose to make alphaB and gammaB indices available
          ! *****************************************************
          dim2=dim1
          call mem_alloc(tmp2,dim2)
          ! tmp2(gammaB, j, b, alphaB) = tmp1^T(b, alphaB; gammaB, j)
          m = nocvAOS*dimAlpha      ! dimension of "row" in tmp1 array (to be "column" in tmp2
          n = noccEOS*dimGamma      ! dimension of "column" in tmp1 array (to be "row" in tmp2)

          call mat_transpose( m, n, 1.0E0_realk, tmp1, 0.0E0_realk, tmp2)
          call mem_dealloc(tmp1)

          !print *,"norm2(tmp4):", norm2(tmp2)

          ! Transform gamma batch index to occupied index
          ! *********************************************
          ! tmp4(i,j,b,alphaB) = sum_{gamma in gammaBatch} CoccT(i,gamma) tmp2(gamma,j,b,alphaB)
          m = nocvAOS
          k = dimGamma
          n = nocvAOS*noccEOS*dimAlpha
          dim1 = i8*m*n
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,CocvAOST(:,GammaStart:GammaEnd),tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2)

          !print *,"norm2(tmp5):", norm2(tmp1)
        
          ! Transform alpha batch index to unoccupied index and update output integral
          ! **************************************************************************
          ! ijba(i,j,b,a) =+ sum_{alpha in alphaBatch} tmp1(i,j,b,alpha)  CoccEOS(alpha,a)
          m = nocvAOS*noccEOS*nocvAOS
          k = dimAlpha
          n = noccEOS
          call dec_simple_dgemm_update(m,k,n,tmp1,CoccEOST(:,AlphaStart:AlphaEnd),pjqi, 'n', 't')
          call mem_dealloc(tmp1)

          ! Note: To have things consecutive in memory it is better to pass CunoccT to the dgemm
          ! routine and then transpose (insted of passing CoccEOS and not transpose).

       end do BatchAlpha
    end do BatchGamma

    !> Reorder tmp(p,j,q,i) -> tmp(i,j,p,q)
    !call array_reorder_4d(1.0E0_realk,pjqi,noccEOS,noccEOS,nocvAOS,nocvAOS,[4,2,1,3],0.0E0_realk,ijpq)
    
    !print *,"norm2(pjqi)", norm2(pjqi)
    !print *,"norm2(ijpq)", norm2(ijpq)

    !print *, ijpq

    ! Free and nullify stuff
    ! **********************
    nullify(mysetting%LST_GAB_LHS)
    nullify(mysetting%LST_GAB_RHS)
    call free_decscreen(DECSCREEN)

    ! Free gamma batch stuff
    call mem_dealloc(orb2batchGamma)
    call mem_dealloc(batchdimGamma)
    call mem_dealloc(batchsizeGamma)
    call mem_dealloc(batchindexGamma)
    orb2batchGamma => null()
    batchdimGamma => null()
    batchsizeGamma => null()
    batchindexGamma => null()
    do idx=1,nbatchesGamma
       call mem_dealloc(batch2orbGamma(idx)%orbindex)
       batch2orbGamma(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbGamma)
    batch2orbGamma => null()

    ! Free alpha batch stuff
    call mem_dealloc(orb2batchAlpha)
    call mem_dealloc(batchdimAlpha)
    call mem_dealloc(batchsizeAlpha)
    call mem_dealloc(batchindexAlpha)
    orb2batchAlpha => null()
    batchdimAlpha => null()
    batchsizeAlpha => null()
    batchindexAlpha => null()
    do idx=1,nbatchesAlpha
       call mem_dealloc(batch2orbAlpha(idx)%orbindex)
       batch2orbAlpha(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbAlpha)
    batch2orbAlpha => null()

    ! Free F12 related pointers
    call mem_dealloc(CoccEOST)
    call mem_dealloc(CocvAOST)
    
  end subroutine get_mp2f12_MOmatrix_ijpq

  !> Brief: Get <ij|OPERATOR|pq> integrals stored in the order (i,j,p,q).
  !> Author: Yang M. Wang
  !> Data: June 2013
  subroutine get_mp2f12_MOmatrix_ijmc(MySetting,nbasis,noccEOS,ncabsAO,ncabsMO,CoccEOS,Ccabs,mjci,INTSPEC) 
    implicit none

    !> Integrals settings
    type(lssetting), intent(inout) :: mysetting
    !> Number of basis functions AO
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals MO in EOS
    integer,intent(in) :: noccEOS
    !> Number of Cabs AO
    integer,intent(in) :: ncabsAO
    !> Number of Cabs MO
    integer,intent(in) :: ncabsMO

    !> Occupied MO coefficients
    real(realk),intent(in),dimension(nbasis,noccEOS) :: CoccEOS
    !> Cabs MO Coefficients
    real(realk),intent(in),dimension(ncabsAO,ncabsMO) :: Ccabs
        
    !>  <ij |OPERATOR| mc> integrals stored in the order (m,j,c,i)
    real(realk), intent(inout) :: mjci(noccEOS,noccEOS,ncabsMO,noccEOS)   

    integer :: alphaB,gammaB,dimAlpha,dimGamma,GammaStart, GammaEnd, AlphaStart, AlphaEnd
    real(realk),pointer :: tmp1(:),tmp2(:),CoccEOST(:,:),CcabsT(:,:)
    integer(kind=long) :: dim1,dim2
    integer :: m,k,n,idx
    logical :: FullRHS,doscreen
    integer :: MaxActualDimAlpha,nbatchesAlpha,nbatches,MaxActualDimGamma,nbatchesGamma,iorb
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
    type(batchtoorb), pointer :: batch2orbAlpha(:),batch2orbGamma(:)
#ifdef VAR_OMP
    integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
    TYPE(DECscreenITEM)   :: DecScreen
    Character            :: intSpec(5)
    integer :: MinAObatchSize, MaxAObatchSize, GammaBatchSize, AlphaBatchSize

    ! ***********************************************************
    ! For efficiency when calling dgemm, save transposed matrices
    ! ***********************************************************
    call mem_alloc(CoccEOST,noccEOS,nbasis)
    call mat_transpose(nbasis,noccEOS,1.0E0_realk,CoccEOS,0.0E0_realk,CoccEOST)

    call mem_alloc(CcabsT,ncabsMO,ncabsAO)
    call mat_transpose(ncabsAO,ncabsMO,1.0E0_realk,Ccabs,0.0E0_realk,CcabsT)
  
    ! ************************************************************
    ! Allocate mem space for a temporary V1 that will be reordered
    ! ************************************************************
    ! call mem_alloc(pjqi,nocvAOS,noccEOS,nocvAOS,noccEOS)

    ! ************************
    ! Determine AO batch sizes
    ! ************************
    ! NOTE: Ideally the batch sizes should be optimized according to the available memory
    ! (as is done e.g. in get_optimal_batch_sizes_for_mp2_integrals).
    ! For simplicity we simply choose the gamma batch to contain all basis functions,
    ! while we make the alpha batch as small as possible

    ! Minimum AO batch size
    call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,MinAObatchSize,'R')

    ! Maximum AO batch size (all basis functions)
    MaxAObatchSize = nbasis
    ! Setting AO batch size (all basis functions) For testing!
    MinAObatchSize = nbasis

    ! Set alpha and gamma batch size as written above
    ! GammaBatchSize = MaxAObatchSize
    GammaBatchSize = MaxAObatchSize
    AlphaBatchSize = MinAObatchSize
    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
    ! Orbital to batch information
    ! ----------------------------

    call mem_alloc(orb2batchGamma,nbasis)
    call build_batchesofAOS(DECinfo%output,mysetting,GammaBatchSize,nbasis,MaxActualDimGamma,&
         & batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,orb2BatchGamma,'R')
    ! Batch to orbital information
    ! ----------------------------
    call mem_alloc(batch2orbGamma,nbatchesGamma)
    do idx=1,nbatchesGamma
       call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
       batch2orbGamma(idx)%orbindex = 0
       batch2orbGamma(idx)%norbindex = 0
    end do
    do iorb=1,nbasis
       idx = orb2batchGamma(iorb)
       batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
       K = batch2orbGamma(idx)%norbindex
       batch2orbGamma(idx)%orbindex(K) = iorb
    end do

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************
    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchAlpha,nbasis)
    call build_batchesofAOS(DECinfo%output,mysetting,AlphaBatchSize,nbasis,&
         & MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,'R')

    ! Batch to orbital information
    ! ----------------------------
    call mem_alloc(batch2orbAlpha,nbatchesAlpha)
    do idx=1,nbatchesAlpha
       call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
       batch2orbAlpha(idx)%orbindex = 0
       batch2orbAlpha(idx)%norbindex = 0
    end do
    do iorb=1,nbasis
       idx = orb2batchAlpha(iorb)
       batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
       K = batch2orbAlpha(idx)%norbindex
       batch2orbAlpha(idx)%orbindex(K) = iorb
    end do

    ! Integral screening stuff
    doscreen = Mysetting%scheme%cs_screen .or. Mysetting%scheme%ps_screen
    !doscreen = .FALSE.
    call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mysetting,&
         & nbatches,nbatchesAlpha,nbatchesGamma,INTSPEC)
    IF(doscreen)then
       call II_getBatchOrbitalScreen(DecScreen,mysetting,&
            & nbasis,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,DECinfo%output,DECinfo%output)
    endif
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
    mjci = 0.0_realk

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
       dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
       GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
       GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
 
       BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
          dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
          AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
          AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch

          ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
          ! ************************************************************************************
          dim1 = i8*nbasis*ncabsAO*dimAlpha*dimGamma   ! dimension for integral array

          call mem_alloc(tmp1,dim1)
          tmp1 = 0.0_realk

          ! Store integral in tmp1(1:dim1) array in (beta,delta,alphaB,gammaB) order
          IF(doscreen) mysetting%LST_GAB_RHS => DECSCREEN%masterGabRHS
          IF(doscreen) mysetting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p
          call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
               & mysetting, tmp1, batchindexAlpha(alphaB), batchindexGamma(gammaB), &
               & batchsizeAlpha(alphaB), batchsizeGamma(gammaB), nbasis, ncabsAO, dimAlpha, dimGamma,FullRHS,&
               & nbatches,'RCRRG')

          print *,"FullRHS:", FullRHS
          print *,"norm2(tmp1):", norm2(tmp1)
          call mem_dealloc(tmp1)

          STOP "wangy hack"
          
!!$          ! Transform beta to occupied index "j".
!!$          ! *************************************
!!$          ! Note: ";" indicates the place where the array is transposed:
!!$          ! tmp2(delta,alphaB,gammaB,j) = sum_{beta} tmp1^T(beta;delta,alphaB,gammaB) Cocc_{beta j}
!!$          m = ncabsAO*dimGamma*dimAlpha   ! # elements in "delta alphaB gammaB" dimension of tmp1^T
!!$          k = nbasis                      ! # elements in "beta" dimension of tmp1^T
!!$          n = noccEOS                     ! # elements in second dimension of Cocc
!!$          dim2 = i8*m*n ! dimension of tmp2 array
!!$
!!$          call mem_alloc(tmp2,dim2)
!!$          tmp2 = 0.0_realk
!!$
!!$          call dec_simple_dgemm(m,k,n,tmp1,CoccEOS,tmp2, 't', 'n')
!!$          call mem_dealloc(tmp1)
!!$
!!$          print *,"norm2(tmp2):", norm2(tmp2)
!!$          !call mem_dealloc(tmp2)
!!$
!!$          ! Transform beta to unoccupied index "b".
!!$          ! ***************************************
!!$          ! tmp3(b,alphaB,gammaB,j) = sum_{delta} CunoccT(b,delta) tmp2(delta,alphaB,gammaB,j)
!!$          ! Note: We have stored the transposed Cunocc matrix, so no need to transpose in
!!$          ! the call to dgemm.
!!$          m = ncabsMO
!!$          k = ncabsAO
!!$          n = dimAlpha*dimGamma*noccEOS
!!$          dim1 = i8*m*n  ! dimension of tmp2 array
!!$          call mem_alloc(tmp1,dim1)
!!$          tmp1 = 0.0_realk
!!$
!!$          print *,"m, k, n:", m, k, n 
!!$          print *, "shape(tmp1)", shape(tmp1)
!!$          print *, "shape(Ccabs)", shape(Ccabs)
!!$          print *, "shape(tmp2)", shape(tmp2)
!!$
!!$          print *,"norm2(Ccabs) before dgemm:", norm2(Ccabs)
!!$          print *,"norm2(tmp2) before dgemm:", norm2(tmp2)
!!$
!!$          call dec_simple_dgemm(m,k,n,Ccabs,tmp2,tmp1,'t','n')
!!$
!!$          print *,"norm2(tmp2) after dgemm:", norm2(tmp2)
!!$          print *,"norm2(Ccabs) after dgemm:", norm2(Ccabs)
!!$          print *,"-----------------------------------"
!!$          print *,"norm2(tmp3):", norm2(tmp1)
!!$ 
!!$
!!$          print *,"ncabsMO, ncabsAO:", ncabsMO, ncabsAO
!!$          print *,"noccEOS nbasis", noccEOS, nbasis
!!$         
!!$          !print *,"-----------------------------------"
!!$          !print *,"norm2(Ccabs):", norm2(Ccabs)
!!$          
!!$          call mem_dealloc(tmp2) 
!!$          call mem_dealloc(tmp1)


!!$          ! Transpose to make alphaB and gammaB indices available
!!$          ! *****************************************************
!!$          dim2=dim1
!!$          call mem_alloc(tmp2,dim2)
!!$          ! tmp2(gammaB, j, b, alphaB) = tmp1^T(b, alphaB; gammaB, j)
!!$          m = ncabsMO*dimAlpha      ! dimension of "row" in tmp1 array (to be "column" in tmp2
!!$          n = noccEOS*dimGamma      ! dimension of "column" in tmp1 array (to be "row" in tmp2)
!!$
!!$          call mat_transpose(tmp1, m, n, tmp2)
!!$          call mem_dealloc(tmp1)
!!$
!!$          print *,"norm2(tmp4):", norm2(tmp2)
!!$
!!$          ! Transform gamma batch index to occupied index
!!$          ! *********************************************
!!$          ! tmp4(i,j,b,alphaB) = sum_{gamma in gammaBatch} CoccT(i,gamma) tmp2(gamma,j,b,alphaB)
!!$          m = noccEOS
!!$          k = dimGamma
!!$          n = noccEOS*ncabsMO*dimAlpha
!!$          dim1 = i8*m*n
!!$          call mem_alloc(tmp1,dim1)
!!$          call dec_simple_dgemm(m,k,n,CoccEOST(:,GammaStart:GammaEnd),tmp2,tmp1, 'n', 'n')
!!$          call mem_dealloc(tmp2)
!!$
!!$          print *,"norm2(tmp5):", norm2(tmp1)
!!$        
!!$          ! Transform alpha batch index to unoccupied index and update output integral
!!$          ! **************************************************************************
!!$          ! ijba(i,j,b,a) =+ sum_{alpha in alphaBatch} tmp1(i,j,b,alpha)  CoccEOS(alpha,a)
!!$          m = noccEOS*noccEOS*ncabsMO
!!$          k = dimAlpha
!!$          n = noccEOS
!!$          call dec_simple_dgemm_update(m,k,n,tmp1,CoccEOST(:,AlphaStart:AlphaEnd),mjci,'n','t')
          
          ! call mem_dealloc(tmp1)
          ! Note: To have things consecutive in memory it is better to pass CunoccT to the dgemm
          ! routine and then transpose (insted of passing CoccEOS and not transpose).

       end do BatchAlpha
    end do BatchGamma

    print *,

    ! Free and nullify stuff
    ! **********************
    nullify(mysetting%LST_GAB_LHS)
    nullify(mysetting%LST_GAB_RHS)
    call free_decscreen(DECSCREEN)

    ! Free gamma batch stuff
    call mem_dealloc(orb2batchGamma)
    call mem_dealloc(batchdimGamma)
    call mem_dealloc(batchsizeGamma)
    call mem_dealloc(batchindexGamma)
    orb2batchGamma => null()
    batchdimGamma => null()
    batchsizeGamma => null()
    batchindexGamma => null()
    do idx=1,nbatchesGamma
       call mem_dealloc(batch2orbGamma(idx)%orbindex)
       batch2orbGamma(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbGamma)
    batch2orbGamma => null()

    ! Free alpha batch stuff
    call mem_dealloc(orb2batchAlpha)
    call mem_dealloc(batchdimAlpha)
    call mem_dealloc(batchsizeAlpha)
    call mem_dealloc(batchindexAlpha)
    orb2batchAlpha => null()
    batchdimAlpha => null()
    batchsizeAlpha => null()
    batchindexAlpha => null()
    do idx=1,nbatchesAlpha
       call mem_dealloc(batch2orbAlpha(idx)%orbindex)
       batch2orbAlpha(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbAlpha)
    batch2orbAlpha => null()

    ! Free F12 related pointers
    call mem_dealloc(CoccEOST)
    call mem_dealloc(CcabsT)
    
  end subroutine get_mp2f12_MOmatrix_ijmc

  !> Brief: Fock matrix elements
  !> Author: Yang M. Wang
  !> Data: August 2013
  subroutine get_mp2f12_Fij(Fij, MyLsitem, MyFragment, nocc, noccfull, nvirt, nbasis, ncabsAO)
    implicit none

    !> Full molecule info
    type(ccatom), intent(in) :: MyFragment
    type(lsitem), intent(inout) :: mylsitem
    integer :: nocc, noccfull, nvirt, nbasis, ncabsAO
    type(matrix) :: Fcc
    type(matrix) :: Fii
    type(matrix) :: Dmat

    

    !> Fock Occupied MO coefficients
    real(realk), intent(inout) :: Fij(nocc,nocc)   
    
    

    ! Mixed AO/AO full MO Fock matrix
    print *, "nbasis, ncabsAO", nbasis, ncabsAO
    call mat_init(Fcc,nbasis,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'RRRRC')    
    
    !Fii
    call mat_init(Fii, nocc, nocc)
    !call MO_transform_AOMatrix(mylsitem, nbasis, nocc, noccfull, nvirt,&
    !     & MyFragment%ypo, MyFragment%ypv,'ii',Fcc,Fii)

    !call mat_to_full(Fii,1.0E0_realk,Fij)
 
    !print *, norm2(Fii%elms)
    
    call mat_free(Fcc)
    call mat_free(Fii)
  !  call mat_free(Dmat)

  end subroutine get_mp2f12_Fij

  

end module f12_integrals_module
 
