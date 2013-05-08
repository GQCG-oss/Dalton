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

#ifdef VAR_LSMPI
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

  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
#ifdef VAR_LSMPI
  use decmpi_module !, only: mpi_communicate_mp2_int_and_amp
#endif

  use dec_fragment_utils!,only: calculate_fragment_memory, &
  !       & dec_simple_dgemm_update,start_flop_counter,&
  !       & end_flop_counter, dec_simple_dgemm, mypointer_init, &
  !       & get_currently_available_memory, atomic_fragment_free
  use array2_simple_operations!, only: array2_free, array2_extract_EOS, &
  !       & get_mp2_integral_transformation_matrices, get_mp2_integral_transformation_matrices_fc, &
  !      & extract_occupied_eos_mo_indices, extract_virtual_EOS_MO_indices,array2_init,array2_print
  use array4_simple_operations!, only: array4_delete_file, array4_init_file, &
  !       & array4_init_standard, array4_free, array4_reorder, array4_init, &
  !       & array4_contract1, array4_open_file, array4_write_file_type2, &
  !       & array4_close_file, array4_write_file_type1, mat_transpose, &
  !     & array4_read_file_type2

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
#ifdef MOD_UNRELEASED
    !> MO coefficient matrix for the occupied EOS
    real(realk), pointer :: CoccEOS(:,:)
    
    !> F12 integrals for the V_term: (gr)_ij^ij integrals
    real(realk), pointer :: gr_ijji(:,:,:,:)

    !> number of AO orbitals
    integer :: number_basis
    
    !> number of occupied MO orbitals in EOS 
    integer :: noccEOS

    integer :: ix, i, j
    real(realk) :: energy

    noccEOS = MyFragment%noccEOS
    call mem_alloc(CoccEOS, MyFragment%number_basis, noccEOS)
    call mem_alloc(gr_ijji, noccEOS, noccEOS, noccEOS, noccEOS)

    ! ***********************************************************
    ! Creating a CoccEOS matrix 
    ! ***********************************************************
    do i=1, MyFragment%noccEOS
       ix = MyFragment%idxo(i)
       CoccEOS(:,i) = MyFragment%ypo(:,ix)
    end do

    print *, "---------------get_f12_V_gr_ijij------------- "

    ! Get integrals <ij|gr|ij> stored as gr(i,j,j,i)
    call get_f12_V_gr_ijij(MyFragment%MyLsitem%Setting,MyFragment%number_basis,noccEOS,CoccEOS,gr_ijji)

    print *, "---------------E_P(V1)------------- "

    call mp2f12_E21_single_fragment(gr_ijji, noccEOS, energy)

    print *, energy
    myfragment%energies(14) = energy

    ! Free memory
    call mem_dealloc(CoccEOS)
    call mem_dealloc(gr_ijji)
#endif
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
#ifdef MOD_UNRELEASED
    real(realk), pointer :: CoccEOS(:,:)
    real(realk), pointer :: gr_ijji(:,:,:,:)

    integer :: number_basis
    integer :: noccEOS
    integer :: ix, i, j
    real(realk) :: energy

    logical,pointer :: dopair_occ(:,:)
    call mem_alloc(dopair_occ,noccEOS,noccEOS)

    noccEOS = Pairfragment%noccEOS
    call mem_alloc(CoccEOS, Pairfragment%number_basis, noccEOS)
    call mem_alloc(gr_ijji, noccEOS, noccEOS, noccEOS, noccEOS)

    ! ***********************************************************
    ! Creating a CoccEOS matrix 
    ! ***********************************************************
    do i=1, PairFragment%noccEOS
       ix = PairFragment%idxo(i)
       CoccEOS(:,i) = PairFragment%ypo(:,ix)
    end do

    ! Get integrals <ij|gr|ij> stored as gr(i,j,j,i)
    call get_f12_V_gr_ijij(PairFragment%MyLsitem%Setting,PairFragment%number_basis,noccEOS,CoccEOS,gr_ijji)

    print *, "---------------E_PQ(V2)------------- "

    call mp2f12_E21_pair_fragment(gr_ijji,  Fragment1, Fragment2, PairFragment, noccEOS, energy)

    print *, energy
    pairfragment%energies(14) = energy

    ! Input for Dec Driver
    ! myfragment%energies(14) = energy

    ! Free memory
    call mem_dealloc(dopair_occ)
    call mem_dealloc(CoccEOS)
    call mem_dealloc(gr_ijji)
#endif
  end subroutine f12_pair_fragment_energy

#ifdef MOD_UNRELEASED
  !> Brief: MP2-F12 correction for the single fragment of term V1: E_P(V1)
  !> Author: Yang M. Wang
  !> Data: April 2013
  subroutine mp2f12_E21_single_fragment(gr_ijij, nocc, energy)
    implicit none

    real(realk),intent(out) :: energy
    real(realk),intent(in)  :: gr_ijij(nocc,nocc,nocc,nocc)
    integer,intent(in)      :: nocc
    !
    integer     :: i,j
    real(realk) :: tmp

    tmp = 0E0_realk
    do i=1,nocc
       tmp = tmp + gr_ijij(i,i,i,i)
    enddo

    energy = -0.5E0_realk*tmp

    tmp = 0E0_realk
    do j=1,nocc
       do i=j+1,nocc
          tmp = tmp + 5E0_realk * gr_ijij(i,j,j,i) - gr_ijij(i,j,i,j)
       enddo
    enddo

    energy = energy - 0.25E0_realk*tmp
    energy = 2E0_realk*energy

  end subroutine mp2f12_E21_single_fragment


  !> Brief: MP2-F12 correction for the pair fragment of term V1: E_PQ(V1) 
  !> Author: Yang M. Wang
  !> Data: April 2013
  subroutine mp2f12_E21_pair_fragment(gr_ijij, Fragment1, Fragment2, PairFragment, noccEOS, energy)
    implicit none

    !> Fragment 1 in the pair fragment
    type(ccatom),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(ccatom),intent(inout) :: Fragment2
    !> Pair fragment formed from fragment 1 and 2
    type(ccatom), intent(inout) :: PairFragment
    !> The pair fragment energy
    real(realk),intent(out) :: energy
    !> The integrals
    real(realk),intent(in)  :: gr_ijij(noccEOS,noccEOS,noccEOS,noccEOS)
    !> Number of occupied MO orbitals in EOS space
    integer,intent(in)      :: noccEOS
    !
    integer     :: i,j
    real(realk) :: tmp
    logical,pointer :: dopair_occ(:,:)

    tmp = 0.0_realk

    call mem_alloc(dopair_occ,noccEOS,noccEOS)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)

    do i=1,noccEOS
       do j=1,noccEOS

          if( dopair_occ(i,j) ) then  !DoPair1and2   
             tmp = tmp + 5E0_realk * gr_ijij(i,j,j,i) - gr_ijij(i,j,i,j)
          endif

       enddo
    enddo

    energy = -0.25E0_realk*tmp

    !Free memory
    call mem_dealloc(dopair_occ)

  end subroutine mp2f12_E21_pair_fragment

  !> Brief: Get <ij|gr|ij> integrals stored in the order (i,j,j,i).
  !> Author: Yang M. Wang
  !> Data: April 2013
  subroutine get_f12_V_gr_ijij(MySetting,nbasis,noccEOS,CoccEOS,ijji)
    implicit none

    !> Integrals settings
    type(lssetting), intent(inout) :: mysetting
    !> Number of basis functions AO
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals MO
    integer,intent(in) :: noccEOS
    !> Occupied MO coefficients
    real(realk),intent(in),dimension(nbasis,noccEOS) :: CoccEOS
    !>  <ij |gr| ij) integrals stored in the order (i,j,j,i)
    real(realk), intent(inout) :: ijji(noccEOS,noccEOS,noccEOS,noccEOS)

    integer :: alphaB,gammaB,dimAlpha,dimGamma,GammaStart, GammaEnd, AlphaStart, AlphaEnd
    real(realk),pointer :: tmp1(:),tmp2(:),CoccEOST(:,:)
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
    call mat_transpose(CoccEOS,nbasis,noccEOS,CoccEOST)

    ! ************************
    ! Determine AO batch sizes
    ! ************************
    ! NOTE: Ideally the batch sizes should be optimized according to the available memory
    ! (as is done e.g. in get_optimal_batch_sizes_for_mp2_integrals).
    ! For simplicity we simply choose the gamma batch to contain all basis functions,
    ! while we make the alpha batch as small as possible

    ! Minimum AO batch size
    call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,MinAObatchSize)

    ! Maximum AO batch size (all basis functions)
    MaxAObatchSize = nbasis

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
         & batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,orb2BatchGamma)

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
         & MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha)

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

    ! *****************
    ! Set integral info
    ! *****************
    INTSPEC(1)='R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)='R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)='R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)='R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)='F' !F = f_12 r_12^-1 operator

    !call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRF')
    !call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
    !   &                          MyMolecule%ypo, MyMolecule%ypv,'iiii',gAO,Fijkl)

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
    ijji = 0.0_realk

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
               & mysetting, tmp1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
               & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),nbasis,nbasis,dimAlpha,dimGamma,FullRHS,&
               & nbatches,INTSPEC)

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

          ! Transform beta to unoccupied index "b".
          ! ***************************************
          ! tmp1(b,alphaB,gammaB,j) = sum_{delta} CunoccT(b,delta) tmp2(delta,alphaB,gammaB,j)
          ! Note: We have stored the transposed Cunocc matrix, so no need to transpose in
          ! the call to dgemm.
          m = noccEOS
          k = nbasis
          n = dimAlpha*dimGamma*noccEOS
          dim1 = i8*noccEOS*noccEOS*dimAlpha*dimGamma  ! dimension of tmp2 array
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,CoccEOST,tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2)

          ! Transpose to make alphaB and gammaB indices available
          ! *****************************************************
          dim2=dim1
          call mem_alloc(tmp2,dim2)
          ! tmp2(gammaB,j,b,alphaB) = tmp1^T(b,alphaB;gammaB,j)
          m = noccEOS*dimAlpha    ! dimension of "row" in tmp1 array (to be "column" in tmp2)
          n = noccEOS*dimGamma      ! dimension of "column" in tmp1 array (to be "row" in tmp2)
          call mat_transpose(tmp1,m,n,tmp2)
          call mem_dealloc(tmp1)

          ! Transform gamma batch index to occupied index
          ! *********************************************
          ! tmp1(i,j,b,alphaB) = sum_{gamma in gammaBatch} CoccT(i,gamma) tmp2(gamma,j,b,alphaB)
          m = noccEOS
          k = dimGamma
          n = noccEOS*noccEOS*dimAlpha
          dim1 = i8*noccEOS*noccEOS*noccEOS*dimAlpha
          call mem_alloc(tmp1,dim1)
          call dec_simple_dgemm(m,k,n,CoccEOST(:,GammaStart:GammaEnd),tmp2,tmp1, 'n', 'n')
          call mem_dealloc(tmp2)

          ! Transform alpha batch index to unoccupied index and update output integral
          ! **************************************************************************
          ! ijba(i,j,b,a) =+ sum_{alpha in alphaBatch} tmp1(i,j,b,alpha)  CoccEOS(alpha,a)
          m = noccEOS*noccEOS*noccEOS
          k = dimAlpha
          n = noccEOS
          call dec_simple_dgemm_update(m,k,n,tmp1,CoccEOST(:,AlphaStart:AlphaEnd),ijji, 'n', 't')
          call mem_dealloc(tmp1)
          ! Note: To have things consecutive in memory it is better to pass CunoccT to the dgemm
          ! routine and then transpose (insted of passing Cunocc and not transpose).

       end do BatchAlpha
    end do BatchGamma

    print *, "---------------Inside get_f12_V_gr_ijij-------------- "

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

    call mem_dealloc(CoccEOST)

  end subroutine get_f12_V_gr_ijij
#endif

end module f12_integrals_module
