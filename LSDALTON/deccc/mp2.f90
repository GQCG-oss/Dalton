!> @file
!> Two-electron integrals and amplitudes for MP2.
!> \ author Kasper Kristensen

module mp2_module

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


  !> Calculate integrals and amplitudes for MP2 energy (and possibly first order properties)
  interface MP2_integrals_and_amplitudes
     !> Integrals for energy and MP2 amplitudes
     module procedure MP2_integrals_and_amplitudes_energy
     !> Integrals for energy, integrals for first properties, and MP2 amplitudes
     module procedure MP2_integrals_and_amplitudes_energy_and_first_order
  end interface


contains



  !> \brief Get maximum batch dimension encountered in integral program.
  !> \author Kasper Kristensen
  !> \date February 2011
  function max_batch_dimension(mylsitem,nbasis) result(maxdim)

    implicit none

    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of basis function
    integer, intent(in) :: nbasis
    integer :: maxdim
    integer, pointer :: orb2batch(:), batchdim(:)
    integer :: i, nbatches

    ! Initialize stuff
    nullify(orb2batch)
    nullify(batchdim)
    call mem_alloc(orb2batch,nbasis)

    ! Get batch info
    call II_getBatchOrbitalInfo(mylsitem%setting,nbasis,&
         & orb2Batch,nbatches,DECinfo%output,DECinfo%output)

    ! Vector containing dimensions for each batch
    call mem_alloc(batchdim,nbatches)
    batchdim = 0
    do i=1,nbasis
       batchdim(orb2batch(i)) = batchdim(orb2batch(i))+1
    end do

    ! Find maximum batch dimension
    maxdim=0
    do i=1,nbatches
       if( batchdim(i) > maxdim ) maxdim=batchdim(i)
    end do

    ! Clean up
    call mem_dealloc(batchdim)
    call mem_dealloc(orb2batch)

  end function max_batch_dimension



  !> \brief Workhorse for calculating EOS integrals and EOS amplitudes for MP2 calculation -
  !> both for occupied and virtual partitioning schemes.
  !> If requested, also EOS integrals for first order MP2 properties are calculated.
  !> See index convention inside subroutine!
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine MP2_integrals_and_amplitudes_workhorse(MyFragment,goccEOS, toccEOS, &
       & gvirtEOS, tvirtEOS, djik,blad,bat,first_order_integrals)

    implicit none

    ! For frozen core all occupied indices refer to only the valence space EXCEPT
    ! for first_order_integrals where the "k" index for blak and djik is BOTH core+valence
    ! (See below)

    !> Atomic fragment (or pair fragment)
    type(ccatom), intent(inout) :: MyFragment
    !> Integrals for occ EOS: (d j|c i) in the order (d,j,c,i) [see notation inside]
    type(array4),intent(inout) :: goccEOS
    !> Amplitudes for occ EOS in the order (d,j,c,i) [see notation inside]
    type(array4),intent(inout) :: toccEOS
    !> Integrals for virt EOS: (b l|a k) in the order (b,l,a,k) [see notation inside]
    !> NOTE: If first_order_integrals AND frozen core, then l is valence and k is both core+valence
    type(array4),intent(inout) :: gvirtEOS
    !> Amplitudes for virt EOS in the order (b,l,a,k) [see notation inside]
    type(array4),intent(inout) :: tvirtEOS
    !> Occ EOS integrals (d j | i k) in the order (d,j,i,k)  [only if first_order_integrals is true]
    !> NOTE: If first_order_integrals AND frozen core, then j is valence and k is both core+valence
    type(array4),intent(inout) :: djik
    !> Virt EOS integrals (b l | a d) in the order (b,l,a,d)  [only if first_order_integrals is true]
    type(array4),intent(inout) :: blad
    !>Batch sizes used for MP2 integral/amplitude calculation
    !> For MPI:  master rank - this is determined inside subroutine
    !> For MPI:  slave rank - this is determined based on input (effectively intent(in))
    type(mp2_batch_construction),intent(inout) :: bat
    !> Determines whether djik and blad are calculated (true)
    !> or not (false --> djik and blad not initialized)
    logical,intent(inout) :: first_order_integrals
    type(array2) :: CDIAGocc, CDIAGvirt, Uocc, Uvirt
    real(realk), pointer :: EVocc(:), EVvirt(:)
    integer :: nbasis,nocc,nvirt, noccEOS, nvirtEOS
    integer :: alpha,gamma,alphaB,gammaB,dimAlpha,dimGamma
    real(realk), pointer :: UvirtEOST(:,:)
    type(mypointer) :: tmp1,tmp2,tmp3,tmp4
    type(mypointer),pointer :: b1(:),b2(:),b3(:)
    real(realk),pointer :: UoccEOST(:,:),  UoccEOS(:,:), UvirtEOS(:,:),CvirtT(:,:),gvirt2(:,:,:,:)
    real(realk),pointer :: gocc(:,:,:,:),tocc(:,:,:,:),gvirt(:,:,:,:),tvirt(:,:,:,:),CoccT(:,:)
    integer,pointer :: V(:,:)
    integer(kind=long) :: dim1,dim2,dim3,dim4,idx,idx2,max1,max2,max3,maxdim,start
    integer:: Astart, Aend,dimA, A,B,I,J,counter,siz,arrsize
    real(realk) :: flops
    integer,dimension(4) :: dimocc, dimvirt
    integer :: m,k,n, nvbatches, Abat, GammaStart, GammaEnd, AlphaStart, AlphaEnd,c,d,l
    real(realk) :: deltaEPS
    integer :: iorb,thread_idx,nthreads
    integer,dimension(4) :: dims
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    Character(80)        :: FilenameCS,FilenamePS
    logical :: FoundInMem,FullRHS,doscreen
    real(realk),pointer :: VVVO(:,:,:,:), OOOV(:,:,:,:)
    real(realk),pointer :: LvirtEOST(:,:), LoccEOST(:,:), LvirtT(:,:)
    type(array2) :: LoccEOS,LvirtEOS, tmparray2, LoccTALL,CDIAGoccTALL,UoccALL
    real(realk) :: tcpu, twall,tcpuTOT,twallTOT,tcpu1,twall1,tcpu2,twall2
    real(realk) :: tcpu_start,twall_start, tcpu_end,twall_end
    integer :: MaxActualDimAlpha,nbatchesAlpha,nbatches
    integer :: MaxActualDimGamma,nbatchesGamma
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
#ifdef VAR_OMP
    integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
    TYPE(DECscreenITEM)   :: DecScreen
    logical :: master,wakeslave
    real(realk) :: twmpi1,twmpi2, tcmpi1, tcmpi2, tmpidiff
#ifdef VAR_LSMPI
    INTEGER(kind=ls_mpik) :: HSTATUS
    CHARACTER*(MPI_MAX_PROCESSOR_NAME) ::  HNAME
!    this really should be
!    character*(MPI_MAX_PROCESSOR_NAME) :: HNAME
    integer,pointer :: decmpitasks(:)
    integer(kind=ls_mpik) :: masterrank = 0
#endif
    integer(kind=ls_mpik) :: ierr
    integer :: myload,ncore
    real(realk),pointer :: arr(:)
    integer :: num,extra,narrays,nocctot
    type(mypointer),pointer :: CvirtTspecial(:,:)
    real(realk),pointer :: mini1(:),mini2(:),mini3(:),mini4(:)
    logical :: ts,fc
    Character            :: intSpec(5)
    myload = 0


! If MPI is not used, consider the single node to be "master"
master=.true.
#ifdef VAR_LSMPI
if(infpar%lg_mynum /= 0) then  ! this is a local slave
master=.false.
end if
#endif


    ! Lots of timings
    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpuTOT,twallTOT,DECinfo%output)
    call LSTIMER('START',tcpu_start,twall_start,DECinfo%output)

    doscreen = MyFragment%mylsitem%setting%scheme%cs_screen.OR.&
         & MyFragment%mylsitem%setting%scheme%ps_screen


    ! Note on the indices
    ! *******************
    !
    ! Occupied AOS, diagonal basis: I,J
    ! Virtual  AOS, diagonal basis: A,B
    ! Occupied EOS, local basis: i,j
    ! Occupied AOS, local basis: k,l
    ! Virtual  EOS, local basis: a,b
    ! Virtual  AOS, local basis: c,d
    ! AO indices: alpha,beta,gamma,delta
    ! Batch (AO) indices : alphaB, gammaB
    !
    ! The indices are connected to the two-electron distribution (Mulliken notation) as:
    ! ( beta delta | alpha gamma) = g(beta,delta,alpha,gamma)
    !
    ! beta  <--> B
    ! delta <--> J
    ! alpha <--> A
    ! gamma <--> I
    !
    ! Occupied partitioning (two occ EOS, two virt AOS):
    ! B <--> d
    ! J <--> j
    ! A <--> c
    ! I <--> i
    !
    ! Virtual partitioning (two occ EOS, two virt AOS):
    ! B <--> b
    ! J <--> l
    ! A <--> a
    ! I <--> k
    !



    ! Initialize stuff
    nullify(orb2batchAlpha)
    nullify(batchdimAlpha)
    nullify(batchsizeAlpha)
    nullify(batch2orbAlpha)
    nullify(batchindexAlpha)
    nullify(orb2batchGamma)
    nullify(batchdimGamma)
    nullify(batchsizeGamma)
    nullify(batch2orbGamma)
    nullify(batchindexGamma)
    nullify(mini1,mini2,mini3,mini4)
    nbasis = MyFragment%number_basis
    nocc = MyFragment%noccAOS   ! occupied AOS (only valence for frozen core)
    nvirt = MyFragment%nunoccAOS   ! virtual AOS
    noccEOS = MyFragment%noccEOS  ! occupied EOS
    nvirtEOS = MyFragment%nunoccEOS  ! virtual EOS
    nocctot = MyFragment%nocctot     ! total occ: core+valence (identical to nocc without frozen core)
    ncore = MyFragment%ncore   ! number of core orbitals
    ! For frozen core energy calculation, we never need core orbitals
    ! (but we do if first order integrals are required)
    if(DECinfo%frozencore .and. (.not. first_order_integrals)) nocctot = nocc

    ! In general, for frozen core AND first order integrals, special care must be taken
    ! No frozen core OR frozen core calculation for just energy uses the same
    ! code from now on because the frozen core approximation is "built into" the fragment,
    if(DECinfo%frozencore .and. first_order_integrals) then
       fc = .true.
    else
       fc =.false.
    end if

    if(first_order_integrals) then
       if(DECinfo%PL>0) write(DECinfo%output,*) 'Calculating MP2 integrals (both energy and density) and MP2 amplitudes...'
    else
       if(DECinfo%PL>0) write(DECinfo%output,*) 'Calculating MP2 integrals (only energy) and MP2 amplitudes...'
    end if
    if(MyFragment%nEOSatoms==2) then ! pair fragment
       if(master) write(DECinfo%output,'(a,3i8)') '#PAIRDIMS# basis,occ,virt ', nbasis,nocc,nvirt
    else ! single fragment
       if(master) write(DECinfo%output,'(a,3i8)') '#SINGLEDIMS# basis,occ,virt ', nbasis,nocc,nvirt
    end if


    ! Size of EOS arrays used for updating inside integral loop
    ! *********************************************************

    ! occupied EOS dimension during integral loop (different from output dimensions!)
    dimocc=[nvirt,noccEOS,noccEOS,nvirt]
    call mem_alloc(gocc,dimocc(1),dimocc(2),dimocc(3),dimocc(4) )  ! occ EOS integrals
    call mem_alloc(tocc,dimocc(1),dimocc(2),dimocc(3),dimocc(4) )  ! occ EOS amplitudes
    gocc=0E0_realk
    tocc=0E0_realk

    ! virtual EOS dimension during integral loop (different from output dimensions)
    dimvirt=[nvirtEOS,nvirtEOS,nocc,nocc]
    call mem_alloc(tvirt,dimvirt(1),dimvirt(2),dimvirt(3),dimvirt(4) )  ! virt EOS amplitudes
    ! Special case: Last occupied index is both core+valence!
    call mem_alloc(gvirt,dimvirt(1),dimvirt(2),dimvirt(3),nocctot )  ! virt EOS integrals
    gvirt=0E0_realk
    tvirt=0E0_realk

    ! Arrays used for updating integrals used for first-order MP2 properties
    if(first_order_integrals) then
       call mem_alloc(VVVO,nvirt,nvirtEOS,nvirtEOS,nocc)
       call mem_alloc(OOOV,nocctot,noccEOS,noccEOS,nvirt)
       VVVO=0E0_realk
       OOOV=0E0_realk
    end if


    ! *************************************
    ! Get arrays for transforming integrals
    ! *************************************
    ! CDIAGocc, CDIAGvirt:  MO coefficients for basis where Fock matrix is diagonal
    ! Uocc, Uvirt: Transform from diagonal basis to local basis (and vice versa)
    ! Note: Uocc and Uvirt have indices (local,diagonal)
    call mem_alloc(EVocc,nocc)
    call mem_alloc(EVvirt,nvirt)
    if(fc) then
       ! Special routine for MP2 density/gradient using frozen core approx.
       ! Note: The occupied coefficient matrices CDIAGocc and Uocc and the eigenvalues EVocc
       !       refer only to the valence space, while
       !       LoccTALL, CDIAGoccTALL, and UoccALL refer to both core+valence!
       ! Special note: CDIAGoccTALL have valence orbitals BEFORE core orbitals, see
       !               get_MP2_integral_transformation_matrices_fc!
       call get_MP2_integral_transformation_matrices_fc(MyFragment,CDIAGocc, CDIAGvirt, Uocc, Uvirt, &
            & LoccTALL, CDIAGoccTALL, UoccALL, EVocc, EVvirt)
    else  ! not frozen core or simple energy calculation
       call get_MP2_integral_transformation_matrices(MyFragment,CDIAGocc, CDIAGvirt, Uocc, Uvirt, &
            & EVocc, EVvirt)
       LoccTALL = array2_init([nocc,nbasis])
       call mat_transpose(MyFragment%ypo,nbasis,nocc,LoccTALL%val)
    end if


    ! Extract occupied and virtual EOS indices from rows of Uocc and Uvirt
    call array2_extract_EOS(Uocc,MyFragment,'O','R',tmparray2)
    call mem_alloc(UoccEOS,tmparray2%dims(1),tmparray2%dims(2) )
    UoccEOS=tmparray2%val
    call array2_free(tmparray2)

    call array2_extract_EOS(Uvirt,MyFragment,'V','R',tmparray2)
    call mem_alloc(UvirtEOS,tmparray2%dims(1),tmparray2%dims(2) )
    UvirtEOS=tmparray2%val
    call array2_free(tmparray2)


    ! Extract occupied and virtual EOS indices from columns of MyFragment%ypo
    ! and MyFragment%ypv, i.e. the local EOS molecular orbital coefficients
    LoccEOS = array2_init_plain([nbasis,noccEOS])
    call extract_occupied_EOS_MO_indices(LoccEOS,MyFragment)
    LvirtEOS = array2_init_plain([nbasis,nvirtEOS])
    call extract_virtual_EOS_MO_indices(LvirtEOS,MyFragment)

    ! For efficiency when calling dgemm, save transposed matrices
    ! (Transposition itself is done below)
    call mem_alloc(CoccT,nocc,nbasis)
    call mem_alloc(CvirtT,nvirt,nbasis)
    call mem_alloc(UoccEOST,nocc,noccEOS)
    call mem_alloc(UvirtEOST,nvirt,nvirtEOS)
    call mem_alloc(LoccEOST,noccEOS,nbasis)
    call mem_alloc(LvirtEOST,nvirtEOS,nbasis)
    call mem_alloc(LvirtT,nvirt,nbasis)


    ! Determine optimal batchsizes with available memory
    ! **************************************************
    if(master) then
       call get_optimal_batch_sizes_for_mp2_integrals(MyFragment,first_order_integrals,bat,.true.)
    end if


    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************

    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchGamma,nbasis)
    call build_batchesofAOS(DECinfo%output,MyFragment%mylsitem%setting,bat%MaxAllowedDimGamma,&
         & nbasis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,orb2BatchGamma)

    if(master) write(DECinfo%output,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma

    ! Translate batchindex to orbital index
    ! -------------------------------------
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
    call build_batchesofAOS(DECinfo%output,MyFragment%mylsitem%setting,bat%MaxAllowedDimAlpha,&
         & nbasis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha)
    if(master) write(DECinfo%output,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha

    ! Translate batchindex to orbital index
    ! -------------------------------------
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


    ! **************************************************
    ! * Determine batch information for virtual batch  *
    ! **************************************************

    ! Vector for keeping track of virtual batches:
    nvbatches = ceiling(real(nvirt)/real(bat%virtbatch)) ! Number of virtual batches

    ! Structure of virtual batch vector for batch "i":
    ! V(1,i)=(i-1)*virtbatch + 1      ! start index of batch "i"
    ! V(2,i)=i*virtbatch              ! end index of batch "i"
    !
    ! Special case:  V(2,nvbatches) = nvirt
    !
    ! Thus, all batches have size virtbatch, except the last batch,
    ! which may be smaller (if nvirt/virtbatch is not an integer).

    call mem_alloc(V,2,nvbatches)
    V(1,1)=1
    do i=1,nvbatches
       V(1,i) = bat%virtbatch*(i-1) + 1
       V(2,i) = bat%virtbatch*i
    end do
    V(2,nvbatches)=nvirt
    if(master) write(DECinfo%output,*) 'BATCH: Number of virtual batches =', nvbatches


    ! *************************************************************
    ! *                    Start up MPI slaves                    *
    ! *************************************************************

#ifdef VAR_LSMPI

    ! Only use slave helper if there is at least two jobs AND
    ! there is at least one local slave available.
    if(nbatchesAlpha*nbatchesGamma >1 .and. infpar%lg_nodtot>1) then
       wakeslave=.true.
    else
       wakeslave=.false.
    end if

    ! Master starts up slave
    StartUpSlaves: if(wakeslave .and. master) then


       ! Wake up slaves to do the job: MP2 - integrals and amplitudes  (MP2INAMP)
       call ls_mpibcast(MP2INAMP,infpar%master,infpar%lg_comm)

       ! Sanity check
       if(.not. MyFragment%BasisInfoIsSet) then
          call lsquit('MP2_integrals_and_amplitudes_workhorse: &
               & Basis info for master is not set!',-1)
       end if

       ! Communicate fragment information to slaves
       call mpi_communicate_mp2_int_and_amp(MyFragment,bat,first_order_integrals,.true.)

    end if StartUpSlaves
    HSTATUS = 80
    CALL MPI_GET_PROCESSOR_NAME(HNAME,HSTATUS,IERR)

#endif



    ! Transpose matrices
    ! ******************
    call mat_transpose(CDIAGocc%val,nbasis,nocc,CoccT)
    call mat_transpose(CDIAGvirt%val,nbasis,nvirt,CvirtT)
    call mat_transpose(UoccEOS,noccEOS,nocc,UoccEOST)
    call mat_transpose(UvirtEOS,nvirtEOS,nvirt,UvirtEOST)
    call mat_transpose(LoccEOS%val,nbasis,noccEOS,LoccEOST)
    call mat_transpose(LvirtEOS%val,nbasis,nvirtEOS,LvirtEOST)
    call mat_transpose(MyFragment%ypv,nbasis,nvirt,LvirtT)


    ! ***************************************************************
    ! Make special CvirtT array to avoid passing elements which are
    ! non-consecutive in memory inside integral loop
    ! ***************************************************************
    call mem_alloc(CvirtTspecial,nvbatches,nbatchesAlpha)
    do alphaB = 1,nbatchesAlpha  ! AO batches
       dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
       AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
       AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch

       do Abat=1,nvbatches
          Astart = V(1,Abat)    ! first index in virtual batch
          Aend = V(2,Abat)      ! last index in virtual batch
          dimA = Aend-Astart+1  ! dimension of virtual batch

          call mem_alloc(CvirtTspecial(Abat,alphaB)%p,dimAlpha*dimA)
          counter=0
          do alpha=AlphaStart,AlphaEnd
             do A=Astart,Aend
                counter=counter+1
                CvirtTspecial(Abat,alphaB)%p(counter) = CvirtT(A,alpha)
             end do
          end do

       end do
    end do



    ! Set integral info
    ! *****************
    INTSPEC(1)='R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)='R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)='R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)='R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)='C' !C = Coulomb operator
    call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,MyFragment%mylsitem%setting,&
     &                           nbatches,nbatchesAlpha,nbatchesGamma,INTSPEC)
    IF(doscreen)then
       call II_getBatchOrbitalScreen(DecScreen,MyFragment%mylsitem%setting,&
            & nbasis,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,DECinfo%output,DECinfo%output)
    endif
    !setup LHS screening - the full AO basis is used so we can use the
    !                      full matrices:        FilenameCS and FilenamePS
    !Note that it is faster to calculate the integrals in the form
    !(dimAlpha,dimGamma,nbasis,nbasis) so the full AO basis is used on the RHS
    !but the integrals is stored and returned in (nbasis,nbasis,dimAlpha,dimGamma)
#ifdef VAR_OMP
nthreads=OMP_GET_MAX_THREADS()
if(DECinfo%PL>0) write(DECinfo%output,*) 'Starting DEC-MP2 integral/amplitudes - OMP. Number of threads: ', OMP_GET_MAX_THREADS()
#else
nthreads=1
num=1
if(DECinfo%PL>0) write(DECinfo%output,*) 'Starting DEC-MP2 integral/amplitudes - NO OMP!'
#endif


! For MPI: Get array defining which jobs are done by which ranks
#ifdef VAR_LSMPI
      call mem_alloc(decmpitasks,nbatchesAlpha*nbatchesGamma)
      if(wakeslave) then  ! share workload with slave(s)
         call distribute_mpi_jobs(decmpitasks,nbatchesAlpha,nbatchesGamma,&
              & batchdimAlpha,batchdimGamma,myload)
         if(DECinfo%PL>0) write(DECinfo%output,'(a,i6,a,i15)') 'Rank ', infpar%mynum, ' has load ', myload
      else ! master do all jobs
         decmpitasks=infpar%lg_mynum
      end if
#endif



      ! POINTER INITIALIZATION STUFF
      ! ****************************

      call mem_alloc(b1,nthreads)
      call mem_alloc(b2,nthreads)
      call mem_alloc(b3,nthreads)
      nullify(tmp1%p,tmp2%p,tmp3%p,tmp4%p)
      do i=1,nthreads
         nullify(b1(i)%p,b2(i)%p,b3(i)%p)
      end do

      ! Memory requirement for big array
      max1 = sum(bat%size1(1:3))  ! step 1 in integral/amplitude scheme
      max2 = bat%size2(4) + nthreads*sum(bat%size2(1:3)) ! step 2 in integral/amplitude scheme
      max3 = sum(bat%size3(1:2)) ! step 3 in integral/amplitude scheme

      maxdim=max(max1,max2,max3)
      ! Make maxdim extra large to ensure that all pointers start at 512+integer
      narrays = 1 + nthreads*3   ! number of arrays in step 2
      extra = narrays * 512    ! Extra size of maxdim to ensure this
      maxdim = maxdim + extra

      ! Print for statistics
      if(DECinfo%PL>0) write(DECinfo%output,'(a,4i14)') 'size1 ', bat%size1
      if(DECinfo%PL>0) write(DECinfo%output,'(a,4i14)') 'size2 ', bat%size2
      if(DECinfo%PL>0) write(DECinfo%output,'(a,4i14)') 'size3 ', bat%size3
      if(DECinfo%PL>0) write(DECinfo%output,'(a,3i14)') 'Tot sizes ', max1,max2,max3
      if(DECinfo%PL>0) write(DECinfo%output,*) 'Static array: elms/GB = ', maxdim, real(maxdim)*8.0e-9
      allocate(arr(maxdim),stat=ierr)
      arr=0.0E0_realk
      if(ierr == 0) then
#ifdef VAR_LSMPI
      if(DECinfo%PL>0) write(DECinfo%output,'(a,i7,i15)') 'MP2: Allocation OK for node/dim ', infpar%mynum,maxdim
#else
      if(DECinfo%PL>0) write(DECinfo%output,'(a,i15)') 'MP2: Allocation OK for dim ', maxdim
#endif
   else
#ifdef VAR_LSMPI
      write(DECinfo%output,'(a,i7,i15)') 'MP2: Error in allocation for node/dimm ', infpar%mynum,maxdim
#else
      write(DECinfo%output,'(a,i15)') 'MP2: Error in allocation for dim ', maxdim
#endif
      call lsquit('MP2: Something wrong for big array allocation!',-1)
   end if




      ! Pointers for step 1
      ! -------------------

      ! tmp1 starts pointing to element 1 in arr and has size bat%size1(1)
      start=1
      call mypointer_init(maxdim,arr,start,bat%size1(1),tmp1)
      CALL c_f_pointer(c_loc(arr(tmp1%start)),tmp1%p,[tmp1%N])
      ! tmp2 starts pointing to element tmp1%end+1 in arr and has size bat%size1(2)
      start=tmp1%end+1
      call mypointer_init(maxdim,arr,start,bat%size1(2),tmp2)
      CALL c_f_pointer(c_loc(arr(tmp2%start)),tmp2%p,[tmp2%N])
      ! tmp2 starts pointing to element tmp2%end+1 in arr and has size bat%size1(3)
      start=tmp2%end+1
      call mypointer_init(maxdim,arr,start,bat%size1(3),tmp3)
      CALL c_f_pointer(c_loc(arr(tmp3%start)),tmp3%p,[tmp3%N])


      ! Pointers for step 2
      ! -------------------
      ! Sanity check - size of tmp4 in step 2 cannot exceed size of tmp1+tmp2 in step 1
      if(bat%size2(4) > bat%size1(1) + bat%size2(2)) then
         call lsquit('MP2 integral/amplitudes: tmp4 is larger than tmp1+tmp2',-1)
      end if

      ! tmp4 starts pointing to element 1 in arr and has size bat%size2(4)
      start=1
      call mypointer_init(maxdim,arr,start,bat%size2(4),tmp4)
      CALL c_f_pointer(c_loc(arr(tmp4%start)),tmp4%p,[tmp4%N])
      start = tmp4%end +1

      do j=1,nthreads
         ! tmp array b1 inside OMP loop
         call mypointer_init(maxdim,arr,start,bat%size2(1),b1(j))
         CALL c_f_pointer(c_loc(arr(b1(j)%start)),b1(j)%p,[b1(j)%N])
         start = b1(j)%end + 1

         ! tmp array b2 inside OMP loop
         call mypointer_init(maxdim,arr,start,bat%size2(2),b2(j))
         CALL c_f_pointer(c_loc(arr(b2(j)%start)),b2(j)%p,[b2(j)%N])
         start = b2(j)%end + 1

         ! tmp array b3 inside OMP loop
         call mypointer_init(maxdim,arr,start,bat%size2(3),b3(j))
         CALL c_f_pointer(c_loc(arr(b3(j)%start)),b3(j)%p,[b3(j)%N])
         start = b3(j)%end + 1
      end do


      if(master) call LSTIMER('INIT MP2-INT',tcpu,twall,DECinfo%output)
      call LSTIMER('START',tcmpi1,twmpi1,DECinfo%output)
      if(.not. master) then  ! flop counting for slaves
         call start_flop_counter()
      end if

    ! Start looping over gamma and alpha batches and calculate integrals
    ! ******************************************************************
#ifdef VAR_LSMPI
      if(DECinfo%PL>0) write(DECinfo%output,'(a,g14.4,i7)') 'Memory (GB) available before loop/node ', &
           & (DECinfo%memory - 1.0E-9_realk*mem_allocated_global), infpar%mynum
#endif
      
    FullRHS = (nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)
    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
       dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
       GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
       GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch


    BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
       dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
       AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
       AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch


#ifdef VAR_LSMPI
       ! MPI: Only do (alpha,gamma) contribution if this is a task for this particular rank
       if(decmpitasks((alphaB-1)*nbatchesGamma+gammaB) /= infpar%lg_mynum) cycle
#endif


       ! *********************************************************************
       ! *                      STEP 1 IN INTEGRAL LOOP                      *
       ! *********************************************************************
       ! Step 1 is the calculation of AO integrals and transformation of three
       ! AO indices to MO indices. For the first-order property integrals
       ! all four indices are transformed in step 1.


       ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
       ! ************************************************************************************
       dim1 = i8*nbasis*nbasis*dimAlpha*dimGamma   ! dimension for integral array
       ! Store integral in tmp1(1:dim1) array in (beta,delta,alphaB,gammaB) order
       IF(doscreen) MyFragment%mylsitem%setting%LST_GAB_RHS => DECSCREEN%masterGabRHS
       IF(doscreen) MyFragment%mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p

       call LSTIMER('START',tcpu1,twall1,DECinfo%output)
       call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
            & MyFragment%mylsitem%setting, tmp1%p(1:dim1),batchindexAlpha(alphaB),batchindexGamma(gammaB),&
            & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),nbasis,nbasis,dimAlpha,dimGamma,FullRHS,&
            & nbatches,INTSPEC)

       call LSTIMER('START',tcpu2,twall2,DECinfo%output)

       ! Loop over each (alpha,gamma) within the (alphaB,gammaB) batch
       dim3=i8*nbasis*nocc*dimAlpha*dimGamma   ! new dimension of tmp3
       dim2=i8*nvirt*nocc*dimAlpha*dimGamma   ! new dimension of tmp2
       do i=1,dimAlpha*dimGamma

          ! Transform index delta to diagonal occupied index for each (alpha,gamma):
          ! tmp3(beta,J,alpha,gamma) = sum_{delta} tmp1(beta,delta,alpha,gamma) C_{delta J}

          ! NOTE!!! Due to lousy handling of pointers in Fortran it is better to make a
          ! small pointer (mini) which points to a specific part of a larger pointer (tmp)
          ! than to pass elements tmp(idx:idx2) to dgemm.
          ! Therefore, the code gets slightly uglier/more complicated.
          idx=i8*(i-1)*nbasis*nbasis+tmp1%start  ! start index for tmp1
          siz=nbasis*nbasis  ! size of (alpha,gamma) chunk of tmp1
          CALL c_f_pointer(c_loc(arr(idx)),mini1,[siz])  ! make mini1 point to this chunk of tmp1
          idx=i8*(i-1)*nbasis*nocc+tmp3%start ! start index for tmp3
          siz=nbasis*nocc  ! size of (alpha,gamma) chunk of tmp3
          CALL c_f_pointer(c_loc(arr(idx)),mini3,[siz]) ! make mini3 point to this chunk of tmp3
          call dec_simple_dgemm(nbasis,nbasis, nocc, mini1, &
               & CDIAGocc%val, mini3, 'n', 'n')

          ! tmp2(B,J,alpha,gamma) = sum_{beta} C^T_{B beta} tmp3(beta,J,alpha,gamma)
          idx=i8*(i-1)*nvirt*nocc + tmp2%start
          siz=nvirt*nocc
          CALL c_f_pointer(c_loc(arr(idx)),mini2,[siz])
          call dec_simple_dgemm(nvirt,nbasis, nocc, CvirtT,mini3,mini2, 'n', 'n')

       end do


       ! Integrals used for first-order MP2 properties
       ! *********************************************

       FirstOrder1: if(first_order_integrals) then


          ! (d a | b J) integrals stored as (d,a,b,J)
          ! =========================================


          ! Transform diagonal AOS index (B) to local EOS index (b):
          ! tmp1(b,J,alphaB,gammaB) = sum_{B} U_{bB} tmp2(B,J,alphaB,gammaB)
          n = nocc*dimAlpha*dimGamma
          dim1=i8*nvirtEOS*nocc*dimAlpha*dimGamma
          call dec_simple_dgemm(nvirtEOS,nvirt, n, UvirtEOS, tmp2%p(1:dim2), tmp1%p(1:dim1), 'n', 'n')


          ! tmp3(b,J,alphaB,a) = sum_{gamma in gammaB} tmp1(b,J,alphaB,gamma) L_{gamma a}
          ! (*) NOTE:  Even though tmp3 still has only four indices (b,J,alphaB,a),
          !            then behind the curtain it is also linked to a specific gammaB batch due
          !            due to the restricted summation "gamma in gammaB".
          !            Therefore, the "a" index in tmp1(B,J,alphaB,a) is not fully transformed,
          !            and only when the contributions are summed up at the end does it make
          !            since to talk about an virtual EOS index a.
          ! Also note that L_{gamma a} = (L^T_{a gamma})^T    [double transposition]
          ! By using this we can pass elements to dgemm which are stored consecutively in memory.
          m = nvirtEOS*nocc*dimAlpha
          dim3=i8*nvirtEOS*nocc*dimAlpha*nvirtEOS
          call dec_simple_dgemm(m,dimGamma, nvirtEOS, tmp1%p(1:dim1), &
               & LvirtEOST(1:nvirtEOS,GammaStart:GammaEnd), tmp3%p(1:dim3), 'n', 't')

          ! Reorder: tmp3(b,J,alphaB,a) --> tmp1(alphaB,a,b,J)
          dim1=dim3
          call mat_transpose(tmp3%p(1:dim3),nvirtEOS*nocc,dimAlpha*nvirtEOS,tmp1%p(1:dim1))

          ! Update: VVVO(d,a,b,J) += sum_{alpha in alphaB} L^T_{d alpha} tmp1(alpha,a,b,J)
          ! (Similarly to the comment above, the d index is only partly transformed by this,
          !  and we only have a true "d" index at the end when all contributions have been added).
          n = nvirtEOS*nvirtEOS*nocc
          call dec_simple_dgemm_update(nvirt,dimAlpha, n, LvirtT(1:nvirt,AlphaStart:AlphaEnd), &
               & tmp1%p(1:dim1), VVVO(:,:,:,:), 'n', 'n')



          ! (k i | j B) integrals stored as (k,i,j,B)
          ! =========================================


          ! Reorder: tmp2(B,J,alphaB,gammaB) --> tmp1(J,B,alphaB,gammaB)
          dim1=dim2
          do counter=1,dimGamma*dimAlpha
             idx=i8*(counter-1)*nvirt*nocc + tmp2%start
             siz = nvirt*nocc
             CALL c_f_pointer(c_loc(arr(idx)),mini2,[siz])
             idx=i8*(counter-1)*nvirt*nocc + tmp1%start
             CALL c_f_pointer(c_loc(arr(idx)),mini1,[siz])
             call mat_transpose(mini2,nvirt,nocc,mini1)
          end do

          ! tmp3(j,B,alphaB,gammaB) = sum_{J} U_{jJ} tmp1(J,B,alphaB,gammaB)
          n = nvirt*dimAlpha*dimGamma
          dim3=i8*noccEOS*nvirt*dimAlpha*dimGamma
          call dec_simple_dgemm(noccEOS,nocc, n, UoccEOS, tmp1%p(1:dim1), tmp3%p(1:dim3), 'n', 'n')

          ! tmp1(j,B,alphaB,i) = sum_{gamma in gammaB} tmp3(j,B,alphaB,gamma) L_{gamma i}
          m = noccEOS*nvirt*dimAlpha
          dim1 = i8*noccEOS*nvirt*dimAlpha*noccEOS
          ! Note: Use that L_{gamma i} = (L^T_{i gamma})^T to pass only elements stored
          !       consecutively in memory to dgemm.
          call dec_simple_dgemm(m,dimGamma,noccEOS, tmp3%p(1:dim3), &
               & LoccEOST(1:noccEOS,GammaStart:GammaEnd), tmp1%p(1:dim1), 'n', 't')

          ! Reorder: tmp1(j,B,alphaB,i) --> tmp3(alphaB,i,j,B)
          dim3=dim1
          call mat_transpose(tmp1%p(1:dim1),noccEOS*nvirt,dimAlpha*noccEOS,tmp3%p(1:dim3))

          ! Update: OOOV(k,i,j,B) += sum_{alpha in alphaB} L^T_{k alpha} tmp3(alpha,i,j,B)
          ! NOTE! "k" refers to BOTH core and valence, also for frozen core approximation
          n = noccEOS*noccEOS*nvirt
          call dec_simple_dgemm_update(nocctot,dimAlpha,n, LoccTALL%val(1:nocctot,AlphaStart:AlphaEnd), &
               & tmp3%p(1:dim3), OOOV(:,:,:,:), 'n', 'n')


       end if FirstOrder1



       ! Integrals used for MP2 energy
       ! *****************************

       ! tmp3(B,J,alphaB,I) = sum_{gamma in gammaB} tmp2(B,J,alphaB,gamma) C_{gamma I}
       ! (Same comment as (*) above)
       ! Note: C_{gamma I} = (C^T_{I gamma})^T  (double transposition)
       !       It is better to used the elements stored in the transposed matrix CoccT, since
       !       then we only pass elements which are stored consecutively in memory to dgemm.
       !       Note: For frozen core "J" is only valence, while "I" is core+valence
       !
       m = nvirt*nocc*dimAlpha
       dim3 = i8*nvirt*nocc*dimAlpha*nocctot  ! New dimension for tmp3
       if(fc) then
          call dec_simple_dgemm(m,dimGamma, nocctot, tmp2%p(1:dim2), &
               & CDIAGoccTALL%val(1:nocctot,GammaStart:GammaEnd), tmp3%p(1:dim3), 'n', 't')
       else
          call dec_simple_dgemm(m,dimGamma, nocctot, tmp2%p(1:dim2), &
               & CoccT(1:nocctot,GammaStart:GammaEnd), tmp3%p(1:dim3), 'n', 't')
       end if

       ! Transition from step 1 to step 2 in integral loop
       ! =================================================

       ! Reorder: tmp3(B,J,alphaB,I) --> tmp4(alphaB,B,J,I)
       dim4=i8*dimAlpha*nvirt*nocc*nocctot
       do counter=1,nocctot
             idx=i8*(counter-1)*nvirt*nocc*dimAlpha + tmp3%start
             siz = nvirt*nocc*dimAlpha
             CALL c_f_pointer(c_loc(arr(idx)),mini3,[siz])
             idx=i8*(counter-1)*nvirt*nocc*dimAlpha + tmp4%start
             CALL c_f_pointer(c_loc(arr(idx)),mini4,[siz])
             call mat_transpose(mini3,nvirt*nocc,dimAlpha,mini4)
       end do

       ! tmp4 will now be used in each step for each thread in step 2



       ! *********************************************************************
       ! *                      STEP 2 IN INTEGRAL LOOP                      *
       ! *********************************************************************
       ! Step 2 is the virtual batching where the final AO-->MO transformations
       ! are carried out and the MP2 amplitudes are determined.


call mem_TurnONThread_Memory()
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(Abat,Astart,Aend,dimA,m,n,siz,ts,&
!$OMP dim1,dim2,dim3,counter,alpha,A,B,i,j,idx,idx2,deltaeps,num,mini1,mini2,mini3,mini4) &
!$OMP SHARED(nvbatches,V,nvirt,nocc,nocctot,dimAlpha,alphaB,noccEOS,dim4,gocc,tocc, &
!$OMP gvirt,tvirt,nvirtEOS,b1,b2,b3,tmp4,CvirtTspecial,UoccEOS,UoccEOST,UvirtEOS, &
!$OMP EVocc, EVvirt,arr,fc)


call init_threadmemvar()

#ifdef VAR_OMP
num = OMP_GET_THREAD_NUM() +1 ! Start counting from 1
#else
num=1
#endif
ts=.true.

!$OMP DO SCHEDULE(dynamic,1)


       BatchA: do Abat=1,nvbatches
          Astart = V(1,Abat)    ! first index in virtual batch
          Aend = V(2,Abat)      ! last index in virtual batch
          dimA = Aend-Astart+1  ! dimension of virtual batch

          ! b1(Abat,B,J,I) = sum_{alpha in alphaB} C^T_{A alpha} tmp4(alpha,B,J,I)
          ! Note: Similarly to the comment above, even though b1 contains only four indices,
          ! behind the curtain it belongs to specific alphaB and gammaB batches, and therefore
          ! it formally has six indices.
          ! Note: For frozen core/gradient "I" is both core+valence here (but changes below).
          n = nvirt*nocc*nocctot
          dim1=i8*dimA*n
          ! Avoid passing elements which are non-consecutive in memory
          call dec_simple_dgemm(dimA,dimAlpha,n,CvirtTspecial(Abat,alphaB)%p, &
               & tmp4%p(1:dim4),b1(num)%p(1:dim1), 'n', 'n',use_thread_safe=ts)


          ! Transform from diagonal to local basis: Two-electron integrals, OCC partitioning
          ! ********************************************************************************
          ! Note: Here "I" is ONLY valence! (But that changes for VIRT partitioning below).
          ! This is the reason for the special ordering of occupied orbitals in CDIAGoccTALL!
          ! Now all core orbitals are listed LAST for the I-index in b1(Abat,B,J,I), and therefore
          ! we may simply access b(1:dimA*nvirt*nocc*nocc) to consider only valence orbitals.
          ! Set b1 dimension such that core orbitals are not considered:
          if(fc) then
             ! This only applies for first order integrals AND using frozen core approx
             ! (otherwise dim1 does not change)
             dim1=i8*dimA*nvirt*nocc*nocc
          end if

          ! Transform diagonal AOS index I to local EOS index i:
          ! b3(A,B,J,i) = sum_{I} b1(Abat,B,J,I) U^T_{Ii}
          m=dimA*nvirt*nocc
          dim3=i8*dimA*nvirt*nocc*noccEOS    ! dimension of b3
          call dec_simple_dgemm(m,nocc, noccEOS, b1(num)%p(1:dim1), UoccEOST, &
               & b3(num)%p(1:dim3), 'n', 'n',use_thread_safe=ts)

          ! Reorder: b3(Abat,B,J,i) --> b2(J,Abat,B,i)
          dim2=dim3
          do counter=1,noccEOS
             idx=i8*(counter-1)*dimA*nvirt*nocc + b3(num)%start
             siz = dimA*nvirt*nocc
             CALL c_f_pointer(c_loc(arr(idx)),mini3,[siz])
             idx=i8*(counter-1)*dimA*nvirt*nocc + b2(num)%start
             CALL c_f_pointer(c_loc(arr(idx)),mini2,[siz])
             call mat_transpose(mini3,dimA*nvirt,nocc,mini2)
          end do

          ! Transform diagonal AOS index J to local EOS index j:
          ! b3(j,Abat,B,i) = sum_{J} U_{jJ} b2(J,Abat,B,i)
          n=dimA*nvirt*noccEOS
          dim3=i8*noccEOS*dimA*nvirt*noccEOS    ! dimension of b3
          call dec_simple_dgemm(noccEOS,nocc, n, UoccEOS, b2(num)%p(1:dim2), &
               &  b3(num)%p(1:dim3), 'n', 'n',use_thread_safe=ts)


          ! Update gocc(B,i,j,Abat) += b3(j,Abat,B,i)
          ! -----------------------------------------
          idx=0

          do i=1,noccEOS
             do B=1,nvirt
                do A=Astart,Aend ! only over A batch
                   do j=1,noccEOS
                      idx=idx+1
                      gocc(B,i,j,A) = gocc(B,i,j,A) + b3(num)%p(idx)
                   end do
                end do
             end do
          end do



          ! Transform from diagonal to local basis: Two-electron integrals, VIRT partitioning
          ! *********************************************************************************
          ! Now again index "I" is both core+valence for first order integrals/frozen core:
          if(fc) then
             dim1=i8*dimA*nvirt*nocc*nocctot
          end if

          ! Reorder: b1(Abat,B,J,I) --> b3(B,Abat,J,I)
          dim3=dim1
          do counter=1,nocc*nocctot
             idx=i8*(counter-1)*dimA*nvirt + b1(num)%start
             siz = dimA*nvirt
             CALL c_f_pointer(c_loc(arr(idx)),mini1,[siz])
             idx=i8*(counter-1)*dimA*nvirt + b3(num)%start
             CALL c_f_pointer(c_loc(arr(idx)),mini3,[siz])
             call mat_transpose(mini1,dimA,nvirt,mini3)
          end do

          ! Transform diagonal AOS index B to local EOS index b:
          ! b2(b,Abat,J,I) = sum_{B} U_{bB} b3(B,Abat,J,I)
          n=dimA*nocc*nocctot
          dim2=i8*nvirtEOS*dimA*nocctot*nocc    ! dimension of b2
          call dec_simple_dgemm(nvirtEOS,nvirt, n, UvirtEOS, b3(num)%p(1:dim3), &
               & b2(num)%p(1:dim2), 'n', 'n',use_thread_safe=ts)

          ! Reorder: b2(b,Abat,J,I) --> b3(Abat,b,J,I)
          dim3=dim2
          do counter=1,nocc*nocctot
             idx=i8*(counter-1)*nvirtEOS*dimA + b2(num)%start
             siz=nvirtEOS*dimA
             CALL c_f_pointer(c_loc(arr(idx)),mini2,[siz])
             idx=i8*(counter-1)*nvirtEOS*dimA + b3(num)%start
             CALL c_f_pointer(c_loc(arr(idx)),mini3,[siz])
             call mat_transpose(mini2,nvirtEOS,dimA,mini3)
          end do

          ! Transform diagonal AOS index Abat to local EOS index a (but only inside A batch):
          ! (The A-->a transformation is only complete after all A batches are done).
          ! Update: b2(a,b,J,I) += sum_{A in Abat} U_{aA} b3(A,b,J,I)
          n=nvirtEOS*nocc*nocctot
          dim2=i8*nvirtEOS*nvirtEOS*nocc*nocctot
          call dec_simple_dgemm(nvirtEOS,dimA, n, UvirtEOS(1:nvirtEOS,Astart:Aend), &
               & b3(num)%p(1:dim3), b2(num)%p(1:dim2), 'n', 'n',use_thread_safe=ts)

          ! Update gvirt(a,b,J,I) += b2(a,b,J,I)
!$OMP CRITICAL (gvirtupdate)
          counter=0
          do I=1,nocctot
             do J=1,nocc
                do b=1,nvirtEOS
                   do a=1,nvirtEOS
                      counter=counter+1
                      gvirt(a,b,J,I) = gvirt(a,b,J,I) + b2(num)%p(counter)
                   end do
                end do
             end do
          end do
!$OMP END CRITICAL (gvirtupdate)


          ! Solve amplitude equation and transform amplitudes to EOS
          ! ********************************************************

          ! At this point b1 contains the two-electron integrals in the diagonal basis
          ! where the solution to the amplitude equation is trivial -
          ! The amplitudes in the diagonal basis are determined by dividing the
          ! two-electron integral with the corresponding elements of the diagonal Fock matrix:
          !
          ! b1(Abat,B,J,I) = b1(Abat,B,J,I) / (eI + eJ - eA - eB)
          ! where
          !       new b1: amplitudes in the diagonal basis
          !       old b1: two-electron integrals in the diagonal basis
          !       eI,eJ are occupied diagonal Fock matrix elements (orbital energies for full molecule)
          !       eA,eB are virtual diagonal Fock matrix elements (orbital energies for full molecule)
          !
          ! Note: For frozen core first order calculation, we skip core orbitals here!
          !       Since valence orbitals are ordered BEFORE core orbitals for the "I" index
          !       (and J is already only valence orbitals), we simply need to loop from
          !       1 to the number of valence orbitals nocc.
          idx=0
          do I=1,nocc ! only run over valence for frozen core 
             do J=1,nocc
                do B=1,nvirt
                   do A=Astart,Aend
                      idx=idx+1
                      deltaEPS = EVocc(I)+EVocc(J)-EVvirt(A)-EVvirt(B)
                      b1(num)%p(idx)=b1(num)%p(idx)/deltaEPS
                   end do
                end do
             end do
          end do

          ! The amplitudes may now be transformed to local EOS indices for
          ! both the occupied and virtual partitioning schemes - using exactly the
          ! same transformations as were used for the integrals above.
          dim1=i8*dimA*nvirt*nocc*nocc



          ! Transform from diagonal to local basis: Two-electron amplitudes, OCC partitioning
          ! *********************************************************************************

          ! Transform diagonal AOS index I to local EOS index i:
          ! b3(A,B,J,i) = sum_{I} b1(Abat,B,J,I) U^T_{Ii}
          m=dimA*nvirt*nocc
          dim3=i8*dimA*nvirt*nocc*noccEOS    ! dimension of b3
          call dec_simple_dgemm(m,nocc, noccEOS, b1(num)%p(1:dim1), UoccEOST, b3(num)%p(1:dim3), 'n', 'n',use_thread_safe=ts)

          ! Reorder: b3(Abat,B,J,i) --> b2(J,Abat,B,i)
          dim2=dim3
          do counter=1,noccEOS
             idx=i8*(counter-1)*dimA*nvirt*nocc + b3(num)%start
             siz = dimA*nvirt*nocc
             CALL c_f_pointer(c_loc(arr(idx)),mini3,[siz])
             idx=i8*(counter-1)*dimA*nvirt*nocc + b2(num)%start
             CALL c_f_pointer(c_loc(arr(idx)),mini2,[siz])
             call mat_transpose(mini3,dimA*nvirt,nocc,mini2)
          end do

          ! Transform diagonal AOS index J to local EOS index j:
          ! b3(j,Abat,B,i) = sum_{J} U_{jJ} b2(J,Abat,B,i)
          n=dimA*nvirt*noccEOS
          dim3=i8*noccEOS*dimA*nvirt*noccEOS    ! dimension of b3
          call dec_simple_dgemm(noccEOS,nocc, n, UoccEOS, b2(num)%p(1:dim2), b3(num)%p(1:dim3), 'n', 'n',use_thread_safe=ts)



          ! Update tocc(B,i,j,Abat) += b3(j,Abat,B,i)
          ! -----------------------------------------
          idx=0

          do i=1,noccEOS
             do B=1,nvirt
                do A=Astart,Aend ! only over A batch
                   do j=1,noccEOS
                      idx=idx+1
                      tocc(B,i,j,A) = tocc(B,i,j,A) + b3(num)%p(idx)
                   end do
                end do
             end do
          end do


          ! Transform from diagonal to local basis: Two-electron amplitudes, VIRT partitioning
          ! **********************************************************************************

          ! Reorder: b1(Abat,B,J,I) --> b3(B,Abat,J,I)
          dim3=dim1
          do counter=1,nocc*nocc
             idx=i8*(counter-1)*dimA*nvirt + b1(num)%start
             siz = dimA*nvirt
             CALL c_f_pointer(c_loc(arr(idx)),mini1,[siz])
             idx=i8*(counter-1)*dimA*nvirt + b3(num)%start
             CALL c_f_pointer(c_loc(arr(idx)),mini3,[siz])
             call mat_transpose(mini1,dimA,nvirt,mini3)
          end do

          ! Transform diagonal AOS index B to local EOS index b:
          ! b2(b,Abat,J,I) = sum_{B} U_{bB} b3(B,Abat,J,I)
          n=dimA*nocc*nocc
          dim2=i8*nvirtEOS*dimA*nocc*nocc    ! dimension of b2
          call dec_simple_dgemm(nvirtEOS,nvirt, n, UvirtEOS, b3(num)%p(1:dim3), b2(num)%p(1:dim2), 'n', 'n',use_thread_safe=ts)

          ! Reorder: b2(b,Abat,J,I) --> b3(Abat,b,J,I)
          dim3=dim2
          do counter=1,nocc*nocc
             idx=i8*(counter-1)*nvirtEOS*dimA + b2(num)%start
             siz = nvirtEOS*dimA
             CALL c_f_pointer(c_loc(arr(idx)),mini2,[siz])
             idx=i8*(counter-1)*nvirtEOS*dimA + b3(num)%start
             CALL c_f_pointer(c_loc(arr(idx)),mini3,[siz])
             call mat_transpose(mini2,nvirtEOS,dimA,mini3)
          end do

          ! Transform diagonal AOS index Abat to local EOS index a (but only inside A batch):
          ! (The A-->a transformation is only complete after all A batches are done).
          ! Update: b2(a,b,J,I) += sum_{A in Abat} U_{aA} b3(A,b,J,I)
          n=nvirtEOS*nocc*nocc
          dim2=i8*nvirtEOS*nvirtEOS*nocc*nocc
          call dec_simple_dgemm(nvirtEOS,dimA, n, UvirtEOS(1:nvirtEOS,Astart:Aend), &
               & b3(num)%p(1:dim3), b2(num)%p(1:dim2), 'n', 'n',use_thread_safe=ts)

          ! Update tvirt(a,b,J,I) += b2(a,b,J,I)
!$OMP CRITICAL (tvirtupdate)
          counter=0
          do I=1,nocc
             do J=1,nocc
                do b=1,nvirtEOS
                   do a=1,nvirtEOS
                      counter=counter+1
                      tvirt(a,b,J,I) = tvirt(a,b,J,I) + b2(num)%p(counter)
                   end do
                end do
             end do
          end do
!$OMP END CRITICAL (tvirtupdate)

       end do BatchA


!$OMP END DO NOWAIT


call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()

    end do BatchAlpha
 end do BatchGamma


 if(master) call LSTIMER('MP2-INT LOOP',tcpu,twall,DECinfo%output)
 call LSTIMER('START',tcmpi2,twmpi2,DECinfo%output)
 tmpidiff = twmpi2-twmpi1

#ifdef VAR_LSMPI
 if(DECinfo%PL>0) write(DECinfo%output,'(a,i6,i12,g18.8)') 'RANK, LOAD, TIME(s) ', infpar%mynum,myload, tmpidiff
 if(master) write(DECinfo%output,'(1X,a,g18.8)') 'TIME INTEGRALLOOP(s) = ', tmpidiff
#endif

if(.not. master) then
   ! effective time for slaves
   MyFragment%slavetime = tmpidiff
   ! FLOP count for integral loop for slaves
   call end_flop_counter(flops)
end if



#ifdef VAR_LSMPI
call mem_dealloc(decmpitasks)
#endif

 nullify(MyFragment%mylsitem%setting%LST_GAB_LHS)
 nullify(MyFragment%mylsitem%setting%LST_GAB_RHS)
 call free_decscreen(DECSCREEN)
 do alphaB = 1,nbatchesAlpha
    do Abat=1,nvbatches
       call mem_dealloc(CvirtTspecial(Abat,alphaB)%p)
       nullify(CvirtTspecial(Abat,alphaB)%p)
    end do
 end do
 call mem_dealloc(CvirtTspecial)

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


 call mem_dealloc(EVocc)
 call mem_dealloc(EVvirt)
 call mem_dealloc(V)




 ! **********************************************************************************
 ! *                         STEP 3 IN MP2-INTEGRAL SCHEME                          *
 ! **********************************************************************************
 ! Status: Currently gocc and tocc (gvirt and tvirt) contain the occupied (virtual)
 !         EOS integrals but the virtual (occupied) indices are still expressed
 !         in the diagonal basis. All that remains is to transform the virtual (occupied)
 !         diagonal indices for gocc and tocc (gvirt and tvirt) into the local basis
 !         using the unitary transformation matrix Uvirt (Uocc).


 ! Assign temporary arrays with dimensions for step 3
 ! --------------------------------------------------
 start = 1
 call mypointer_init(maxdim,arr,start,bat%size3(1),tmp1)
 CALL c_f_pointer(c_loc(arr(tmp1%start)),tmp1%p,[tmp1%N])
 start = tmp1%end+1
 call mypointer_init(maxdim,arr,start,bat%size3(2),tmp2)
 CALL c_f_pointer(c_loc(arr(tmp2%start)),tmp2%p,[tmp2%N])



 ! OCCUPIED PARTITIONING: Transform virtual diagonal indices to local basis
 ! ************************************************************************


 ! Integrals
 ! =========

 ! Transform: tmp2(d,i,j,A) = sum_{B} U_{dB} gocc(B,i,j,A)
 n=noccEOS*noccEOS*nvirt
 dim2=i8*nvirt*n    ! dimension of tmp2
 call dec_simple_dgemm(nvirt,nvirt,n, Uvirt%val, gocc(:,:,:,:), tmp2%p(1:dim2), 'n', 'n')
 call mem_dealloc(gocc)

 ! Transform: tmp1(d,i,j,c) = sum_{A} tmp2(d,i,j,A) U^T_{Ac}
 dim1=dim2
 m=nvirt*noccEOS*noccEOS
 call dec_simple_dgemm(m,nvirt,nvirt, tmp2%p(1:dim2), Uvirt%val, tmp1%p(1:dim1), 'n', 't')

 ! Put integrals into output array in the correct order
 dimocc = [nvirt,noccEOS,nvirt,noccEOS]   ! Output order
 goccEOS=array4_init(dimocc)
 idx=0
 do c=1,nvirt
    do j=1,noccEOS
       do i=1,noccEOS
          do d=1,nvirt
             idx=idx+1
             goccEOS%val(d,j,c,i) = tmp1%p(idx)
          end do
       end do
    end do
 end do



 ! Amplitudes
 ! ==========

 ! Transform: tmp2(d,i,j,A) = sum_{B} U_{dB} tocc(B,i,j,A)
 n=noccEOS*noccEOS*nvirt
 call dec_simple_dgemm(nvirt,nvirt,n, Uvirt%val, tocc(:,:,:,:), tmp2%p(1:dim2), 'n', 'n')
 call mem_dealloc(tocc)

 ! Transform: tmp1(d,i,j,c) = sum_{A} tmp2(d,i,j,A) U^T_{Ac}
 m=nvirt*noccEOS*noccEOS
 call dec_simple_dgemm(m,nvirt,nvirt, tmp2%p(1:dim2), Uvirt%val, tmp1%p(1:dim1), 'n', 't')


 ! Put amplitudes into output array in the correct order
 toccEOS=array4_init(dimocc)
 idx=0
 do c=1,nvirt
    do j=1,noccEOS
       do i=1,noccEOS
          do d=1,nvirt
             idx=idx+1
             toccEOS%val(d,j,c,i) = tmp1%p(idx)
          end do
       end do
    end do
 end do


 ! VIRTUAL PARTITIONING: Transform occupied diagonal indices to local basis
 ! ************************************************************************


 ! Integrals
 ! =========

 ! Transform: tmp1(a,b,J,k) = sum_{I} gvirt(a,b,J,I) U^T_{Ik}
 ! Note: For frozen core first order calc, "I" and "k" are both core+valence
 m = nvirtEOS*nvirtEOS*nocc
 dim1= i8*m*nocctot
 if(fc) then
    ! Recall that the orbitals for the I index in gvirt(a,b,J,I) are ordered with
    ! valence before core for convenience above.
    ! We therefore first need to put these orbitals back into the (core,valence) order...
    ! This is of course very ugly but it led to several simplifications above so it is
    ! worth it to do one ugly reordering loop here...
    call mem_alloc(gvirt2,nvirtEOS,nvirtEOS,nocc,nocctot)
    do I=1,ncore ! put core orbitals into right position
       gvirt2(:,:,:,I) = gvirt(:,:,:,I+nocc)
    end do
    do I=1,nocc  ! put valence orbitals into right position
       gvirt2(:,:,:,I+ncore) = gvirt(:,:,:,I)
    end do
    call mem_dealloc(gvirt)
    ! Do transformation tmp1(a,b,J,k) = sum_{I} gvirt2(a,b,J,I) U^T_{Ik}
    ! with I and k being core+valence
    call dec_simple_dgemm(m,nocctot,nocctot, gvirt2(:,:,:,:), UoccALL%val, tmp1%p(1:dim1), 'n', 't')
    call mem_dealloc(gvirt2)
 else
    call dec_simple_dgemm(m,nocc,nocc, gvirt(:,:,:,:), Uocc%val, tmp1%p(1:dim1), 'n', 't')
    call mem_dealloc(gvirt)
 end if


 ! Reorder: tmp1(a,b,J,k) --> tmp2(J,a,b,k)
 dim2=dim1
 do counter=1,nocctot
    idx=i8*(counter-1)*nvirtEOS*nvirtEOS*nocc + tmp1%start
    siz = nvirtEOS*nvirtEOS*nocc
    CALL c_f_pointer(c_loc(arr(idx)),mini1,[siz])
    idx=i8*(counter-1)*nvirtEOS*nvirtEOS*nocc + tmp2%start
    CALL c_f_pointer(c_loc(arr(idx)),mini2,[siz])
    call mat_transpose(mini1,nvirtEOS*nvirtEOS, nocc,mini2)
 end do

 ! Transform: tmp1(l,a,b,k) = sum_{J} U_{lJ} tmp2(J,a,b,k)
 n=nvirtEOS*nvirtEOS*nocctot
 call dec_simple_dgemm(nocc,nocc,n, Uocc%val, tmp2%p(1:dim2), tmp1%p(1:dim1), 'n', 'n')

 ! Put integrals into output array in the correct order
 dimvirt = [nvirtEOS,nocc,nvirtEOS,nocctot]   ! Output order
 gvirtEOS=array4_init(dimvirt)
 idx=0
 do k=1,nocctot
    do b=1,nvirtEOS
       do a=1,nvirtEOS
          do l=1,nocc
             idx=idx+1
             gvirtEOS%val(b,l,a,k) = tmp1%p(idx)
          end do
       end do
    end do
 end do


 ! Amplitudes
 ! ==========

  ! Transform: tmp1(a,b,J,k) = sum_{I} tvirt(a,b,J,I) U^T_{Ik}
 m = nvirtEOS*nvirtEOS*nocc
 call dec_simple_dgemm(m,nocc,nocc, tvirt(:,:,:,:), Uocc%val, tmp1%p(1:dim1), 'n', 't')
 call mem_dealloc(tvirt)

 ! Reorder: tmp1(a,b,J,k) --> tmp2(J,a,b,k)
 do counter=1,nocc
    idx=i8*(counter-1)*nvirtEOS*nvirtEOS*nocc + tmp1%start
    siz = nvirtEOS*nvirtEOS*nocc
    CALL c_f_pointer(c_loc(arr(idx)),mini1,[siz])
    idx=i8*(counter-1)*nvirtEOS*nvirtEOS*nocc + tmp2%start
    CALL c_f_pointer(c_loc(arr(idx)),mini2,[siz])
    call mat_transpose(mini1,nvirtEOS*nvirtEOS,nocc,mini2)
 end do

 ! Transform: tmp1(l,a,b,k) = sum_{J} U_{lJ} tmp2(J,a,b,k)
 n=nvirtEOS*nvirtEOS*nocc
 call dec_simple_dgemm(nocc,nocc,n, Uocc%val, tmp2%p(1:dim2), tmp1%p(1:dim1), 'n', 'n')

 ! Put amplitudes into output array in the correct order
 dimvirt = [nvirtEOS,nocc,nvirtEOS,nocc]   ! Output order
 tvirtEOS=array4_init(dimvirt)
 idx=0
 do k=1,nocc
    do b=1,nvirtEOS
       do a=1,nvirtEOS
          do l=1,nocc
             idx=idx+1
             tvirtEOS%val(b,l,a,k) = tmp1%p(idx)
          end do
       end do
    end do
 end do



 ! Finalize integrals used for first order MP2 integrals
 ! *****************************************************

 FirstOrder2: if(first_order_integrals) then


    ! (b l | a d) integrals in the order (b,l,a,d)
    ! ============================================

    ! Transform: tmp1(d,a,b,l) = sum_{J} VVVO(d,a,b,J) U^T_{Jl}
    m=nvirt*nvirtEOS*nvirtEOS
    dim1=i8*m*nocc
    call dec_simple_dgemm(m,nocc,nocc, VVVO(:,:,:,:),Uocc%val,tmp1%p(1:dim1), 'n', 't')
    call mem_dealloc(VVVO)

    ! Put amplitudes into output array in the correct order:
    ! (b l | a d) integrals in the order (b,l,a,d)
    dims=[nvirtEOS,nocc,nvirtEOS,nvirt]
    blad = array4_init(dims)
    idx=0
    do l=1,nocc
       do b=1,nvirtEOS
          do a=1,nvirtEOS
             do d=1,nvirt
                idx=idx+1
                blad%val(b,l,a,d) = tmp1%p(idx)
             end do
          end do
       end do
    end do


    ! (d j | i k) integrals in the order (d,j,i,k)
    ! ============================================

    ! Transform: tmp1(k,i,j,d) = sum_{J} OOOV(k,i,j,B) U^T_{Bd}
    m=nocctot*noccEOS*noccEOS
    dim1=i8*m*nvirt
    call dec_simple_dgemm(m,nvirt,nvirt, OOOV(:,:,:,:),Uvirt%val,tmp1%p(1:dim1), 'n', 't')
    call mem_dealloc(OOOV)

    ! Put amplitudes into output array in the correct order:
    ! (d j | i k) integrals in the order (d,j,i,k)
    ! "k" is both core and valence, also for frozen core approx!
    dims=[nvirt,noccEOS,noccEOS,nocctot]
    djik = array4_init(dims)
    idx=0
    do d=1,nvirt
       do j=1,noccEOS
          do i=1,noccEOS
             do k=1,nocctot
                idx=idx+1
                djik%val(d,j,i,k) = tmp1%p(idx)
             end do
          end do
       end do
    end do


 end if FirstOrder2


 ! Free stuff
 nullify(mini1,mini2,mini3,mini4)
 do i=1,nthreads
    nullify(b1(i)%p)
    nullify(b2(i)%p)
    nullify(b3(i)%p)
 end do
 nullify(tmp1%P)
 nullify(tmp2%P)
 nullify(tmp3%P)
 nullify(tmp4%P)
 deallocate(arr)
 call mem_dealloc(b1)
 call mem_dealloc(b2)
 call mem_dealloc(b3)
 call mem_dealloc(CoccT)
 call mem_dealloc(CvirtT)
 call mem_dealloc(UoccEOST)
 call mem_dealloc(UvirtEOST)
 call mem_dealloc(UoccEOS)
 call mem_dealloc(UvirtEOS)
 call mem_dealloc(LoccEOST)
 call mem_dealloc(LvirtEOST)
 call mem_dealloc(LvirtT)
 call array2_free(CDIAGocc)
 call array2_free(CDIAGvirt)
 call array2_free(Uocc)
 call array2_free(Uvirt)
 call array2_free(LoccEOS)
 call array2_free(LvirtEOS)
 call array2_free(LoccTALL)
if(fc) then
 call array2_free(CDIAGoccTALL)
 call array2_free(UoccALL)
end if



! MPI: Add arrays from master and all slaves to get final output arrays on master
! *******************************************************************************
#ifdef VAR_LSMPI

! If slaves were not invoked
! then we of course skip the addition of different components of the array.
 MPIcollect: if(wakeslave) then

    !PE: CHANGE IN ARGUMENT PASSING IN ACCORDANCE WITH REST OF MPI ROUTINES
    ! Add up contibutions to output arrays using MPI reduce
    !arrsize = goccEOS%dims(1)*goccEOS%dims(2)*goccEOS%dims(3)*goccEOS%dims(4)
    !call lsmpi_local_reduction(goccEOS%val(:,:,:,:),arrsize,0)
    call lsmpi_local_reduction(goccEOS%val(:,:,:,:),goccEOS%dims(1),&
          &goccEOS%dims(2),goccEOS%dims(3),goccEOS%dims(4),masterrank)

    !arrsize = toccEOS%dims(1)*toccEOS%dims(2)*toccEOS%dims(3)*toccEOS%dims(4)
    !call lsmpi_local_reduction(toccEOS%val(:,:,:,:),arrsize,0)
    call lsmpi_local_reduction(toccEOS%val(:,:,:,:),toccEOS%dims(1),&
          &toccEOS%dims(2),toccEOS%dims(3),toccEOS%dims(4),masterrank)

    !arrsize = gvirtEOS%dims(1)*gvirtEOS%dims(2)*gvirtEOS%dims(3)*gvirtEOS%dims(4)
    !call lsmpi_local_reduction(gvirtEOS%val(:,:,:,:),arrsize,0)
    call lsmpi_local_reduction(gvirtEOS%val(:,:,:,:),gvirtEOS%dims(1),&
          &gvirtEOS%dims(2),gvirtEOS%dims(3),gvirtEOS%dims(4),masterrank)

    !arrsize = tvirtEOS%dims(1)*tvirtEOS%dims(2)*tvirtEOS%dims(3)*tvirtEOS%dims(4)
    !call lsmpi_local_reduction(tvirtEOS%val(:,:,:,:),arrsize,0)
    call lsmpi_local_reduction(tvirtEOS%val(:,:,:,:),tvirtEOS%dims(1),&
          &tvirtEOS%dims(2),tvirtEOS%dims(3),tvirtEOS%dims(4),masterrank)

    if(first_order_integrals) then
       !arrsize = djik%dims(1)*djik%dims(2)*djik%dims(3)*djik%dims(4)
       !call lsmpi_local_reduction(djik%val(:,:,:,:),arrsize,0)
       call lsmpi_local_reduction(djik%val(:,:,:,:),djik%dims(1),&
             &djik%dims(2),djik%dims(3),djik%dims(4),masterrank)

       !arrsize = blad%dims(1)*blad%dims(2)*blad%dims(3)*blad%dims(4)
       !call lsmpi_local_reduction(blad%val(:,:,:,:),arrsize,0)
       call lsmpi_local_reduction(blad%val(:,:,:,:),blad%dims(1),&
             &blad%dims(2),blad%dims(3),blad%dims(4),masterrank)
   end if

   if(.not. master) then  ! SLAVE: Done with arrays and fragment
      call array4_free(goccEOS)
      call array4_free(toccEOS)
      call array4_free(gvirtEOS)
      call array4_free(tvirtEOS)
      if(first_order_integrals) call array4_free(djik)
      if(first_order_integrals) call array4_free(blad)
      call atomic_fragment_free(MyFragment)
   end if

   ! FLOP counting
   if(master) then
      flops=0.0E0_realk  ! we want to count only flops from slaves (these were set above)
   end if
   call lsmpi_reduction(flops,infpar%master,infpar%lg_comm)
   if(master) MyFragment%flops_slaves = flops ! save flops for local slaves (not local master)

   ! Total time for all slaves (not local master itself)
   if(master) MyFragment%slavetime=0.0E0_realk
   call lsmpi_reduction(MyFragment%slavetime,infpar%master,infpar%lg_comm)

end if MPIcollect

! Number of MPI tasks (=nalpha*ngamma)
MyFragment%ntasks= nbatchesAlpha*nbatchesGamma

#endif


if(master) then
 call LSTIMER('MP2-INT FIN',tcpu,twall,DECinfo%output)
 call LSTIMER('MP2-INT TOTAL',tcpuTOT,twallTOT,DECinfo%output)
 call LSTIMER('START',tcpu_end,twall_end,DECinfo%output)
end if

end subroutine MP2_integrals_and_amplitudes_workhorse




  !> \brief Get (a i | b j) integrals stored in the order (a,i,b,j).
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine get_VOVO_integrals(mylsitem,nbasis,nocc,nvirt,Cvirt,Cocc,VOVO)

    implicit none


    !> LS item
    type(lsitem), intent(inout) :: mylsitem
    !> Number of basis functions
    integer, intent(in) :: nbasis
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of virtual orbitals
    integer, intent(in) :: nvirt
    !> Occupied orbital coefficients
    type(array2), intent(in) :: Cocc
    !> Virtual orbital coefficients
    type(array2), intent(in) :: Cvirt
    !> (a i | b j) integrals stored in the order (a,i,b,j)
    type(array4),intent(inout) :: VOVO

    ! Get integrals (a i | b j) stored as (i,j,b,a)
    VOVO = array4_init([nocc,nocc,nvirt,nvirt])
    call get_ijba_integrals(mylsitem%setting,nbasis,nocc,nvirt,Cocc%val,Cvirt%val,VOVO%val)
    
    ! Reorder: (i,j,b,a) --> (a,i,b,j)
    call array4_reorder(VOVO,[4,1,3,2])

  end subroutine get_VOVO_integrals



  !> \brief Get (a i | b j) integrals stored in the order (i,j,b,a).
  !> No MPI here, intended to be a simple routine to be used for 
  !> (i) calculations not requiring MPI, (ii) debugging,
  !> (iii) starting point for more advanced routine giving other integrals.
  !> So: PLEASE DO NOT POLLUTE THIS SUBROUTINE, IT SHOULD BE KEPT AS AN EASILY ACCESIBLE
  !>     STARTING POINT FOR MORE ADVANCED ROUTINES!
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine get_ijba_integrals(MySetting,nbasis,nocc,nunocc,Cocc,Cunocc,ijba)

    implicit none

    !> Integrals settings
    type(lssetting), intent(inout) :: mysetting
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer,intent(in) :: nunocc
    !> Occupied MO coefficients
    real(realk),intent(in),dimension(nbasis,nocc) :: Cocc
    !> Unoccupied MO coefficients
    real(realk),intent(in),dimension(nbasis,nunocc) :: Cunocc
    !>  (a i | b j) integrals stored in the order (i,j,b,a)
    real(realk),intent(inout) :: ijba(nocc,nocc,nunocc,nunocc)
    integer :: alphaB,gammaB,dimAlpha,dimGamma,GammaStart, GammaEnd, AlphaStart, AlphaEnd
    real(realk),pointer :: tmp1(:),tmp2(:),CoccT(:,:), CunoccT(:,:)
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
    call mem_alloc(CoccT,nocc,nbasis)
    call mem_alloc(CunoccT,nunocc,nbasis)
    call mat_transpose(Cocc,nbasis,nocc,CoccT)
    call mat_transpose(Cunocc,nbasis,nunocc,CunoccT)



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
    INTSPEC(5)='C' !C = Coulomb operator

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
    ijba = 0.0_realk

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
       n = nocc                       ! # elements in second dimension of Cocc
       dim2 = i8*nocc*nbasis*dimAlpha*dimGamma  ! dimension of tmp2 array
       call mem_alloc(tmp2,dim2)
       call dec_simple_dgemm(m,k,n,tmp1,Cocc,tmp2, 't', 'n')
       call mem_dealloc(tmp1)


       ! Transform beta to unoccupied index "b".
       ! ***************************************
       ! tmp1(b,alphaB,gammaB,j) = sum_{delta} CunoccT(b,delta) tmp2(delta,alphaB,gammaB,j)
       ! Note: We have stored the transposed Cunocc matrix, so no need to transpose in
       ! the call to dgemm.
       m = nunocc
       k = nbasis
       n = dimAlpha*dimGamma*nocc
       dim1 = i8*nocc*nunocc*dimAlpha*dimGamma  ! dimension of tmp2 array
       call mem_alloc(tmp1,dim1)
       call dec_simple_dgemm(m,k,n,CunoccT,tmp2,tmp1, 'n', 'n')
       call mem_dealloc(tmp2)


       ! Transpose to make alphaB and gammaB indices available
       ! *****************************************************
       dim2=dim1
       call mem_alloc(tmp2,dim2)
       ! tmp2(gammaB,j,b,alphaB) = tmp1^T(b,alphaB;gammaB,j)
       m = nunocc*dimAlpha    ! dimension of "row" in tmp1 array (to be "column" in tmp2)
       n = nocc*dimGamma      ! dimension of "column" in tmp1 array (to be "row" in tmp2)
       call mat_transpose(tmp1,m,n,tmp2)
       call mem_dealloc(tmp1)


       ! Transform gamma batch index to occupied index
       ! *********************************************
       ! tmp1(i,j,b,alphaB) = sum_{gamma in gammaBatch} CoccT(i,gamma) tmp2(gamma,j,b,alphaB)
       m = nocc
       k = dimGamma
       n = nocc*nunocc*dimAlpha
       dim1 = i8*nocc*nocc*nunocc*dimAlpha
       call mem_alloc(tmp1,dim1)
       call dec_simple_dgemm(m,k,n,CoccT(:,GammaStart:GammaEnd),tmp2,tmp1, 'n', 'n')
       call mem_dealloc(tmp2)


       ! Transform alpha batch index to unoccupied index and update output integral
       ! **************************************************************************
       ! ijba(i,j,b,a) =+ sum_{alpha in alphaBatch} tmp1(i,j,b,alpha)  Cunocc(alpha,a)
       m = nocc*nocc*nunocc
       k = dimAlpha
       n = nunocc
       call dec_simple_dgemm_update(m,k,n,tmp1,CunoccT(:,AlphaStart:AlphaEnd),ijba, 'n', 't')
       call mem_dealloc(tmp1)
       ! Note: To have things consecutive in memory it is better to pass CunoccT to the dgemm
       ! routine and then transpose (insted of passing Cunocc and not transpose).

    end do BatchAlpha
 end do BatchGamma


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


 call mem_dealloc(CoccT)
 call mem_dealloc(CunoccT)

end subroutine Get_ijba_integrals





  !> \brief Calculate EOS integrals and EOS amplitudes for MP2 calculation -
  !> both for occupied and virtual partitioning schemes.
  !> See index convention in MP2_integrals_and_amplitudes_workhorse.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine MP2_integrals_and_amplitudes_energy(MyFragment,goccEOS, toccEOS,gvirtEOS, tvirtEOS)

    implicit none

    !> Atomic fragment (or pair fragment)
    type(ccatom), intent(inout) :: MyFragment
    !> Integrals for occ EOS: (d j|c i) in the order (d,j,c,i) [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: goccEOS
    !> Amplitudes for occ EOS in the order (d,j,c,i) [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: toccEOS
    !> Integrals for virt EOS: (b l|a k) in the order (b,l,a,k) [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: gvirtEOS
    !> Amplitudes for virt EOS in the order (b,l,a,k) [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: tvirtEOS
    type(array4) :: dummy1,dummy2
    type(mp2_batch_construction) :: bat
    logical :: first_order_integrals

    first_order_integrals=.false. ! just energy

    ! Calculate integrals and amplitudes
    call MP2_integrals_and_amplitudes_workhorse(MyFragment,goccEOS, toccEOS, &
         & gvirtEOS, tvirtEOS, dummy1, dummy2,bat,first_order_integrals)

  end subroutine MP2_integrals_and_amplitudes_energy


  !> \brief Workhorse for calculating EOS integrals and EOS amplitudes for MP2 calculation -
  !> both for occupied and virtual partitioning schemes.
  !> Furthermore, the EOS integrals for first order MP2 properties are calculated.
  !> See index convention in MP2_integrals_and_amplitudes_workhorse.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine MP2_integrals_and_amplitudes_energy_and_first_order(MyFragment,goccEOS, toccEOS, &
       & gvirtEOS, tvirtEOS, djik,blad)

    implicit none

    !> Atomic fragment (or pair fragment)
    type(ccatom), intent(inout) :: MyFragment
    !> Integrals for occ EOS: (d j|c i) in the order (d,j,c,i) [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: goccEOS
    !> Amplitudes for occ EOS in the order (d,j,c,i) [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: toccEOS
    !> Integrals for virt EOS: (b l|a k) in the order (b,l,a,k) [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: gvirtEOS
    !> Amplitudes for virt EOS in the order (b,l,a,k) [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: tvirtEOS
    !> Occ EOS integrals (d j | i k) in the order (d,j,i,k)  [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: djik
    !> Virt EOS integrals (b l | a d) in the order (b,l,a,d)  [see MP2_integrals_and_amplitudes_workhorse]
    type(array4),intent(inout) :: blad
    type(mp2_batch_construction) :: bat
    logical :: first_order_integrals

    first_order_integrals=.true. ! first order properties requested

    ! Calculate integrals and amplitudes
    call MP2_integrals_and_amplitudes_workhorse(MyFragment,goccEOS, toccEOS, &
         & gvirtEOS, tvirtEOS, djik,blad,bat,first_order_integrals)

  end subroutine MP2_integrals_and_amplitudes_energy_and_first_order





  !> \brief Get optimal batch sizes to be used in MP2_integrals_and_amplitudes
  !> for the given available memory.
  !> Note: We multiply the estimated available memory by 85% to be on the safe side.
  !> We calculate separate array sizes for three different steps in MP2_integrals_and_amplitudes.
  !> STEP 1: AO integral part and transformation of three indices
  !> STEP 2: Virtual batching part
  !> STEP 3: Final transformations (diagonal->local basis) after integral loop
  !> \author Kasper Kristensen
  !> \date December 2011
subroutine get_optimal_batch_sizes_for_mp2_integrals(MyFragment,first_order_integrals,bat,printstuff)

  implicit none

  !> Fragment info
  type(ccatom),intent(inout) :: MyFragment
  !> Are integrals needed for first-order properties also requested
  logical,intent(in) :: first_order_integrals
  !> Batch sizes used for MP2 integral/amplitude calculation (see mp2_batch_construction type)
  type(mp2_batch_construction),intent(inout) :: bat
  !> Print memory summary for local master?
  !> (If this subroutine is called by local slave we never print, 
  !> regardless of value of printstuff)
  logical,intent(in) :: printstuff
  real(realk) :: MemoryAvailable, GB, MemoryNeeded
  integer :: noccEOS,nocc,nvirtEOS,nvirt,nbasis,GammaOpt,AlphaOpt,step,nvbatches
  integer :: MaxAObatch, MinAOBatch, MaxVirtBatch, MinVirtBatch,gamma,alpha,A, nthreads
  logical :: doprint
#ifdef VAR_OMP
    integer, external :: OMP_GET_MAX_THREADS
#endif

doprint = printstuff
#ifdef VAR_LSMPI
! Only print for local master
if(infpar%lg_mynum/=0) doprint=.false.
#endif

#ifdef VAR_OMP
nthreads=OMP_GET_MAX_THREADS()
#else
nthreads=1
#endif
  if(DECinfo%PL>0) write(DECinfo%output,*) 'Estimating batch sizes for MP2 integrals/amplitudes.'


  ! Init stuff
  GB = 1.000E9_realk ! 1 GB

  ! For fragment with local orbitals where we really want to use the fragment-adapted orbitals
  ! we need to set nocc and nvirt equal to the fragment-adapted dimensions
  if(DECinfo%fragadapt .and. (.not. MyFragment%fragmentadapted) ) then
     nocc=MyFragment%noccFA
     nvirt=MyFragment%nunoccFA
  else
     nocc=MyFragment%noccAOS
     nvirt=MyFragment%nunoccAOS
  end if
  noccEOS=MyFragment%noccEOS
  nvirtEOS=MyFragment%nunoccEOS
  nbasis = MyFragment%number_basis



  ! Memory currently available
  ! **************************
  call get_currently_available_memory(MemoryAvailable)
  ! Note: We multiply by 85 % to be on the safe side!
  MemoryAvailable = 0.85*MemoryAvailable



  ! Maximum and minimum possible batch sizes
  ! ****************************************

  ! The largest possible AO batch is the number of basis functions
  MaxAObatch = nbasis

  ! The smallest possible AO batch depends on the basis set
  ! (More precisely, if all batches are made as small as possible, then the
  !  call below determines the largest of these small batches).
  call determine_maxBatchOrbitalsize(DECinfo%output,MyFragment%mylsitem%setting,MinAObatch)

  ! The smallest/largest possible virtual batch is simply 1/number of virtual orbitals.
  MinVirtBatch = 1
  MaxVirtBatch = nvirt



  ! Find optimal batch sizes for the available memory
  ! *************************************************

  ! Initialize batch sizes to be the minimum possible and then start increasing sizes below
  bat%MaxAllowedDimAlpha = MinAObatch
  bat%MaxAllowedDimGamma = MinAObatch
  bat%virtbatch = MinVirtBatch
  AlphaOpt=MinAObatch
  GammaOpt=MinAObatch


  ! *********************************************************************
  ! *                      STEP 1 IN INTEGRAL LOOP                      *
  ! *********************************************************************
  step=1

  ! Largest possible gamma batch size
  ! =================================
  GammaLoop: do gamma = MaxAObatch,MinAOBatch,-1

     call max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
          & bat%MaxAllowedDimAlpha, gamma, bat%virtbatch, step, nthreads, bat%size1, MemoryNeeded)

     if(MemoryNeeded < MemoryAvailable) then
        GammaOpt = gamma
        exit
     end if

  end do GammaLoop

  ! If gamma batch size was set manually we use that value instead
  if(DECinfo%ccsdGbatch/=0) then
     if(DECinfo%PL>0) write(DECinfo%output,*) 'Gamma batch size was set manually, use that value instead!'
     GammaOpt=DECinfo%ccsdGbatch
  end if

  ! The optimal gamma batch size is GammaOpt.
  ! We now find the maximum possible gamma batch size smaller than or equal to GammaOpt
  ! and store this number in bat%MaxAllowedDimGamma.
  call determine_MaxOrbitals(DECinfo%output,MyFragment%mylsitem%setting,GammaOpt,bat%MaxAllowedDimGamma)

  ! Max size with actual batchsizes
  call max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
       & bat%MaxAllowedDimAlpha,bat%MaxAllowedDimGamma,bat%virtbatch,step,nthreads,bat%size1,MemoryNeeded)
  if(DECinfo%PL>0) write(DECinfo%output,'(1X,a,2i8,g10.3)') 'Optimal/actual gamma size, memory (GB) =', &
       & GammaOpt,bat%MaxAllowedDimGamma,MemoryNeeded



  ! Largest possible alpha batch size
  ! =================================
  AlphaLoop: do alpha = MaxAObatch,MinAOBatch,-1

     call max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
          & alpha, bat%MaxAllowedDimGamma, bat%virtbatch,step,nthreads, bat%size1,MemoryNeeded)

     ! Sanity check: We must ensure that the sum of the sizes of tmp1 and tmp2
     ! in the first step is larger than the size of tmp4 in the second step.
     ! Other we get into memory problem when changing arrays sizes from
     ! size1 to size2 inside the integral loop.
     if(i8*alpha*nocc*nocc*nvirt > bat%size1(1)+bat%size1(2)) then
        if(DECinfo%PL>0) write(DECinfo%output,*) 'WARNING - reducing size1 because tmp4 is too large!'
        cycle
     end if

     if(MemoryNeeded < MemoryAvailable) then
        AlphaOpt = alpha
        exit
     end if

  end do AlphaLoop

  ! If alpha batch size was set manually we use that value instead
  if(DECinfo%ccsdAbatch/=0) then
     if(DECinfo%PL>0) write(DECinfo%output,*) 'Alpha batch size was set manually, use that value instead!'
     AlphaOpt=DECinfo%ccsdAbatch
  end if

  ! Find possible alpha batch size smaller than or equal to AlphaOpt
  call determine_MaxOrbitals(DECinfo%output,MyFragment%mylsitem%setting,AlphaOpt,bat%MaxAllowedDimAlpha)
  call max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
       & bat%MaxAllowedDimAlpha,bat%MaxAllowedDimGamma,bat%virtbatch, step,nthreads,bat%size1,MemoryNeeded)
  if(DECinfo%PL>0) write(DECinfo%output,'(1X,a,2i8,g10.3)') 'Optimal/actual alpha size, memory (GB) =', &
       & AlphaOpt,bat%MaxAllowedDimAlpha,MemoryNeeded



  ! *********************************************************************
  ! *                      STEP 2 IN INTEGRAL LOOP                      *
  ! *********************************************************************
  step=2

  ! Largest possible virtual batch size
  ! ===================================
  ALoop: do A = MaxVirtBatch,MinVirtBatch,-1

     ! Number of virtual batches must be as large as number of OMP threads
     ! to use OMP more efficiently.
     nvbatches = ceiling(real(nvirt)/real(A))
     if(nvbatches < nthreads ) cycle


     call max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
          & bat%MaxAllowedDimAlpha, bat%MaxAllowedDimGamma, A, step,nthreads,bat%size2,MemoryNeeded)

     if(MemoryNeeded < MemoryAvailable) then
        bat%virtbatch = A
        if(DECinfo%PL>0) write(DECinfo%output,'(1X,a,i8,g10.3)') 'Virtual batch size,  memory (GB) =', &
             & bat%virtbatch,MemoryNeeded
        if(DECinfo%PL>0) write(DECinfo%output,'(1X,a,i8)') 'Number of virtual batches =', nvbatches
        exit
     end if

  end do ALoop


  ! *********************************************************************
  ! *                      STEP 3 (AFTER INTEGRAL LOOP)                 *
  ! *********************************************************************

  step=3
     call max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
          & bat%MaxAllowedDimAlpha, bat%MaxAllowedDimGamma, bat%virtbatch, step,nthreads,bat%size3,MemoryNeeded)


  ! Print out and sanity check
  ! ==========================
if(doprint) then
  write(DECinfo%output,*)
  write(DECinfo%output,*)
  write(DECinfo%output,*) '======================================================================='
  write(DECinfo%output,*) '                  MP2 INTEGRALS/AMPLITUDES: MEMORY SUMMARY             '
  write(DECinfo%output,*) '======================================================================='
  write(DECinfo%output,*)
  write(DECinfo%output,'(1X,a,g10.3)') '85% of available memory (GB)            =', MemoryAvailable
  write(DECinfo%output,*)
  write(DECinfo%output,'(1X,a,i8)')    'Number of atomic basis functions        =', nbasis
  write(DECinfo%output,'(1X,a,2i8)')   'Number of occupied orbitals AOS/EOS     =', nocc, noccEOS
  write(DECinfo%output,'(1X,a,2i8)')   'Number of virtual  orbitals AOS/EOS     =', nvirt, nvirtEOS
  write(DECinfo%output,'(1X,a,i8)')    'Maximum alpha batch dimension           =', bat%MaxAllowedDimAlpha
  write(DECinfo%output,'(1X,a,i8)')    'Maximum gamma batch dimension           =', bat%MaxAllowedDimGamma
  write(DECinfo%output,'(1X,a,i8)')    'Maximum virtual batch dimension         =', bat%virtbatch
  write(DECinfo%output,'(1X,a,i8)')    'Number of OMP threads                   =', nthreads
  write(DECinfo%output,*)
end if

  step=1
  call max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
       & bat%MaxAllowedDimAlpha, bat%MaxAllowedDimGamma, bat%virtbatch, step,nthreads,bat%size1,MemoryNeeded)
if(MemoryNeeded > MemoryAvailable) then
   write(DECinfo%output,'(1X,a)') 'STEP 1 in integral loop'
   write(DECinfo%output,'(1X,a)') '-----------------------'
   write(DECinfo%output,'(1X,a,g10.3)') 'Tot memory required for tmp arrays (GB) =', MemoryNeeded
   write(DECinfo%output,'(1X,a,g10.3)') 'Memory required for tmp array 1 (GB)    =', realk*bat%size1(1)/GB
   write(DECinfo%output,'(1X,a,g10.3)') 'Memory required for tmp array 2 (GB)    =', realk*bat%size1(2)/GB
   write(DECinfo%output,'(1X,a,g10.3)') 'Memory required for tmp array 3 (GB)    =', realk*bat%size1(3)/GB
   write(DECinfo%output,*)
   call stats_mem(DECinfo%output)
   call lsquit('get_optimal_batch_sizes_for_mp2_integrals: Estimated array size is &
        & larger than the available memory!',DECinfo%output)
end if

  step=2
  call max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
       & bat%MaxAllowedDimAlpha, bat%MaxAllowedDimGamma, bat%virtbatch, step,nthreads,bat%size2,MemoryNeeded)
  if(MemoryNeeded > MemoryAvailable) then
     write(DECinfo%output,'(1X,a)') 'STEP 2 in integral loop'
     write(DECinfo%output,'(1X,a)') '-----------------------'
     write(DECinfo%output,'(1X,a,g10.3)') 'Tot memory required for tmp arrays (GB) =', MemoryNeeded
     write(DECinfo%output,'(1X,a,g10.3)') 'Memory required for tmp array 1 (GB)    =', realk*bat%size2(1)/GB
     write(DECinfo%output,'(1X,a,g10.3)') 'Memory required for tmp array 2 (GB)    =', realk*bat%size2(2)/GB
     write(DECinfo%output,'(1X,a,g10.3)') 'Memory required for tmp array 3 (GB)    =', realk*bat%size2(3)/GB
     write(DECinfo%output,'(1X,a,g10.3)') 'Memory required for tmp array 4 (GB)    =', realk*bat%size2(4)/GB
     write(DECinfo%output,*)
     call stats_mem(DECinfo%output)
     call lsquit('get_optimal_batch_sizes_for_mp2_integrals: Estimated array size is &
          & larger than the available memory!',DECinfo%output)
  end if

  step=3
  call max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
       & bat%MaxAllowedDimAlpha, bat%MaxAllowedDimGamma, bat%virtbatch, step,nthreads,bat%size3,MemoryNeeded)

  if(MemoryNeeded > MemoryAvailable) then
     write(DECinfo%output,'(1X,a)') 'STEP 3 in integral loop'
     write(DECinfo%output,'(1X,a)') '-----------------------'
     write(DECinfo%output,'(1X,a,g10.3)') 'Tot memory required for tmp arrays (GB) =', MemoryNeeded
     write(DECinfo%output,'(1X,a,g10.3)') 'Memory required for tmp array 1 (GB)    =', realk*bat%size3(1)/GB
     write(DECinfo%output,'(1X,a,g10.3)') 'Memory required for tmp array 2 (GB)    =', realk*bat%size3(2)/GB
     write(DECinfo%output,*)
     call stats_mem(DECinfo%output)
     call lsquit('get_optimal_batch_sizes_for_mp2_integrals: Estimated array size is &
          & larger than the available memory!',DECinfo%output)
  end if


end subroutine get_optimal_batch_sizes_for_mp2_integrals



  !> \brief Get maximum size of each of the four-dimensional arrays used in
  !> MP2_integrals_and_amplitudes for given values of virtual batch size,
  !> alpha batch size, and gamma batch size.
  !> This is done for a given step in the MP2 integral loop.
  !> (See MP2_integrals_and_amplitudes_workhorse for details)
  !> \author Kasper Kristensen
  !> \date February 2011
subroutine max_arraysize_for_mp2_integrals(MyFragment,first_order_integrals,&
     & BatchDimAlpha, BatchdimGamma, BatchdimA, step, nthreads, maxsize, MemoryNeeded)

  implicit none

  !> Fragment info
  type(ccatom),intent(inout) :: MyFragment
  !> Are integrals needed for first-order properties also requested
  logical,intent(in) :: first_order_integrals
  !> Maximum size of AO batch Alpha
  integer, intent(in) :: BatchDimAlpha
  !> Maximum size of AO batch Gamma
  integer, intent(in) :: BatchDimGamma
  !> Maximum size of virtual batch A
  integer, intent(in) :: BatchDimA
  !> Which step in integral loop (see MP2_integrals_and_amplitudes_workhorse for details)
  integer,intent(in) :: step
  !> Number of OMP threads (1 if OMP is not used)
  integer,intent(in) :: nthreads
  !> Maximum array size for each of the four arrays tmp1,tmp2,tmp3, and tmp4
  integer(kind=long),intent(inout) :: maxsize(4)
  !> Memory needed for each of the arrays tmp1,tmp2,tmp3, and tmp4 (measured in GB)
  real(realk),intent(inout) :: MemoryNeeded
  integer :: nocc,nvirt,nbasis,noccEOS,nvirtEOS,i,nocctot
  real(realk) :: GB
  integer(kind=long) :: dim1,dim2,dim3,dim4


  GB = 1.000E9_realk ! 1 GB


  ! Set the relevant dimensions for easy reference
  ! **********************************************
  nocc=MyFragment%noccAOS  ! occupied AOS (only valence for frozen core)
  nvirt = MyFragment%nunoccAOS   ! virtual AOS
  nbasis = MyFragment%number_basis      ! number of basis functions in atomic extent
  noccEOS = MyFragment%noccEOS  ! occupied EOS
  nvirtEOS = MyFragment%nunoccEOS  ! virtual EOS
  if(DECinfo%frozencore .and. first_order_integrals) then
     nocctot = MyFragment%nocctot     ! core+valence
  else
     nocctot = nocc   ! core+valence without frozen core, only valence for frozen core energy calc.
  end if
  ! We know go through the different steps in the integral/amplitude scheme in
  ! MP2_integrals_and_amplitudes to find the largest array size.
  ! Currently this is hardcoded and there are 9 possible candidates for the largest array


  ! Set hardcoded array sizes
  ! *************************

  ! This is currently hardcoded according to the precise structure of MP2_integrals_and_amplitudes_workhorse
  ! --- better solution is under investigation....

  dim1=0
  dim2=0
  dim3=0
  dim4=0


  WhichStep: select case(step)

     ! STEP 1
  case(1)

     dim1=max(dim1,i8*nbasis*nbasis*BatchDimAlpha*BatchDimGamma)
     dim3=max(dim3,i8*nbasis*nocc*BatchDimAlpha*BatchDimGamma)
     dim2=max(dim2,i8*nvirt*nocc*BatchDimAlpha*BatchDimGamma)

     FirstOrder1: if(first_order_integrals) then
        dim3=max(dim3,i8*nvirtEOS*nocc*BatchDimAlpha*nvirtEOS)
        dim3=max(dim3,i8*noccEOS*nvirt*BatchDimAlpha*noccEOS)
     end if FirstOrder1
     dim3=max(dim3,i8*nvirt*nocc*BatchDimAlpha*nocctot)

     ! Memory needed is the sum of the array sizes (no tmp4 here)
     MemoryNeeded=real(dim1) + real(dim2) + real(dim3)

     ! STEP 2
  case(2)

     dim4=i8*BatchDimAlpha*nvirt*nocc*nocctot
     dim1=max(dim1,i8*BatchDimA*nvirt*nocc*nocctot)
     dim3=max(dim3,i8*BatchDimA*nvirt*nocc*nocctot)
     dim2=max(dim2,i8*BatchDimA*nvirt*nocc*noccEOS)
     dim2=max(dim2,i8*nvirtEOS*BatchDimA*nocc*nocctot)
     dim2=max(dim2,i8*nvirtEOS*nvirtEOS*nocc*nocctot)

     ! Memory needed is the sum of the array sizes
     ! However, tmp1,tmp2, and tmp3 are used for EACH thread.
     ! Thus, we multiply these by the number of threads
     MemoryNeeded=real(nthreads)*(real(dim1) + real(dim2) + real(dim3)) +real(dim4)

     ! After loop
  case(3)

     dim1=max(dim1,i8*nvirt*noccEOS*noccEOS*nvirt)
     dim1=max(dim1,i8*nvirtEOS*nvirtEOS*nocc*nocctot)
     dim2=dim1

     if(first_order_integrals) then
        dim1=max(dim1,i8*nvirt*nvirtEOS*nvirtEOS*nocc)
        dim1=max(dim1,i8*nocctot*noccEOS*noccEOS*nvirt)
     end if

     ! Memory needed is the sum of the array sizes
     MemoryNeeded=real(dim1) + real(dim2)

     case default
        call lsquit('max_arraysize_for_mp2_integrals: &
             & step must be 1,2, or 3',DECinfo%output)

  end select WhichStep


  ! Convert total memory requirement to GB
  MemoryNeeded = realk*MemoryNeeded/GB

  ! Set output vector
  maxsize(1) = dim1
  maxsize(2) = dim2
  maxsize(3) = dim3
  maxsize(4) = dim4


end subroutine max_arraysize_for_mp2_integrals



!> \brief Get memory used for updating arrays in MP2_integrals_and_amplitudes.
!> \author Kasper Kristensen
!> \date January 2012
subroutine mp2_integrals_memory_for_updating_arrays(noccEOS,nvirtEOS,noccAOS,nvirtAOS,&
     & MemUpdate,first_order_integrals)
  implicit none
  !> Number of occupied EOS orbitals
  integer,intent(in) :: noccEOS
  !> Number of occupied AOS orbitals
  integer,intent(in) :: noccAOS
  !> Number of virtual EOS orbitals
  integer,intent(in) :: nvirtEOS
  !> Number of virtual AOS orbitals
  integer,intent(in) :: nvirtAOS
  !> Memory used for updating arrays (measured in GB)
  real(realk),intent(inout) :: MemUpdate
  !> Are first order integrals (MP2 density) required?
  logical,intent(in) :: first_order_integrals
  real(realk) :: GB

  GB = 1.0E9_realk

  ! Occupied partitioning (amplitudes + integrals)
  MemUpdate = 2.0E0_realk*real(noccEOS)*real(noccEOS)*real(nvirtAOS)*real(nvirtAOS)
  ! Virtual partitioning (amplitudes + integrals)
  MemUpdate = MemUpdate + 2.0E0_realk*real(noccAOS)*real(noccAOS)*real(nvirtEOS)*real(nvirtEOS)
  ! Two integrals used for first order properties
  if(first_order_integrals) then
     MemUpdate = MemUpdate + real(noccEOS)*real(noccEOS)*real(noccAOS)*real(nvirtAOS)
     MemUpdate = MemUpdate + real(noccAOS)*real(nvirtAOS)*real(nvirtEOS)*real(nvirtEOS)
  end if

  ! Memory used during update loop in GB
  MemUpdate = realk*MemUpdate/GB


end subroutine mp2_integrals_memory_for_updating_arrays

!> \brief Get MO integrals (A I | B J), stored as (A,I,B,J), from full AO integrals
!> Only to be used for testing purposes for full molecule.
!> \author Kasper Kristensen
!> \date May 2012
subroutine get_VOVO_from_full_AO(nbasis,nocc,nvirt,Cocc,Cvirt,gao,gmo)
  implicit none
  !> Number of basis functions
  integer, intent(in) :: nbasis
  !> Number of occupied orbitals
  integer, intent(in) :: nocc
  !> Number of virtual orbitals
  integer, intent(in) :: nvirt
  !> Occupied MO coefficients
  real(realk),dimension(nbasis,nocc),intent(in) :: Cocc
  !> Virtual MO coefficients
  real(realk),dimension(nbasis,nvirt),intent(in) :: Cvirt
  !> AO integrals
  real(realk),intent(in) :: gao(nbasis,nbasis,nbasis,nbasis)
  !> MO integrals (A I | B J), stored as (A,I,B,J)
  real(realk),intent(inout) :: gmo(nvirt,nocc,nvirt,nocc)
  real(realk),pointer :: tmp1(:), tmp2(:)
  integer :: m,n,dim1,dim2


  ! Transform: tmp1(mu,nu,rho,I) = sum_{sigma} gao(mu,nu,rho,sigma) C(sigma,I)
  dim1=nbasis*nbasis*nbasis*nocc
  call mem_alloc(tmp1,dim1)
  m = nbasis**3
  n = nocc
  call dec_simple_dgemm(m,nbasis,n, gao,Cocc, tmp1, 'n', 'n')

  ! Transform: tmp2(B,nu,rho,I) = sum_{mu} C^T_{B mu} tmp1(mu,nu,rho,I)
  dim2=nbasis*nbasis*nvirt*nocc
  call mem_alloc(tmp2,dim2)
  m = nvirt
  n = nbasis*nbasis*nocc
  call dec_simple_dgemm(m,nbasis,n, Cvirt, tmp1, tmp2, 't', 'n')

  ! Reorder: tmp2(B,nu,rho,I) --> tmp1(rho,I,B,nu)
  dim1=dim2
  call mat_transpose(tmp2(1:dim2),nvirt*nbasis,nbasis*nocc,tmp1(1:dim1))


  ! Transform: tmp2(rho,I,B,J) = sum_{nu} tmp1(rho,I,B,nu) C(nu,J)
  dim2=nbasis*nocc*nvirt*nocc
  m = nbasis*nocc*nvirt
  n = nocc
  call dec_simple_dgemm(m,nbasis,n, tmp1(1:dim1),Cocc, tmp2(1:dim2), 'n', 'n')
  call mem_dealloc(tmp1)

  ! Transform: gmo(A,I,B,J) = sum_{mu} C^T_{A rho} tmp2(rho,I,B,J)
  m = nvirt
  n = nocc*nvirt*nocc
  call dec_simple_dgemm(m,nbasis,n, Cvirt, tmp2(1:dim2), gmo(:,:,:,:), 't', 'n')
  call mem_dealloc(tmp2)




end subroutine get_VOVO_from_full_AO


#ifdef VAR_LSMPI

!> \brief Get array defining - for each (alpha,gamma) batch
!> in MP2_integrals_and_amplitudes_workhorse - which rank is supposed to do
!> (alpha,gamma) part of the calculation..
!> \author Kasper Kristensen
!> \March 2012
subroutine get_mpi_tasks_for_MP2_int_and_amp(nalpha,ngamma,decmpitasks)

  implicit none

  !> Number of alpha batches
  integer,intent(in) :: nalpha
  !> Number of gamma batches
  integer,intent(in) :: ngamma
  integer,dimension(nalpha,ngamma),intent(inout) :: decmpitasks
  integer :: alpha,gamma,rank

  ! Example - if nalpha=5 and ngamma=2, and there is one local master and two local slave nodes:
  !
  !                | 0   2 |
  !                | 1   0 |
  ! decmpitasks =  | 2   1 |
  !                | 0   2 |
  !                | 1   0 |
  !
  ! Thus the (alpha,gamma)=(1,1) calculations is done by rank 0 (master)
  !      the (alpha,gamma)=(2,1) calculations is done by rank 1 (first slave)
  !      the (alpha,gamma)=(3,1) calculations is done by rank 2 (second slave)
  ! etc.

  rank=0
  do gamma=1,ngamma
     do alpha=1,nalpha
        decmpitasks(alpha,gamma) = rank
        rank=rank+1

        ! Start counting from local master rank again if number of nodes is exceeded
        if(rank > infpar%lg_nodtot -1) rank=0
     end do
  end do



end subroutine get_mpi_tasks_for_MP2_int_and_amp

#endif


end module mp2_module


#ifdef VAR_LSMPI
!> \brief MPI Slave routine for MP2_integrals_and_amplitudes_workhorse.
!> The slave gets fragment information and other information from master rank,
!> then calls MP2_integrals_and_amplitudes_workhorse to do its specific components
!> of the alpha and gamma loops.
!> \author Kasper Kristensen
!> \March 2012
subroutine MP2_integrals_and_amplitudes_workhorse_slave()

  use precision
  use dec_typedef_module

  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
  use decmpi_module, only: mpi_communicate_mp2_int_and_amp
  use mp2_module,only: MP2_integrals_and_amplitudes_workhorse

  implicit none

  !> Fragment information
  type(ccatom) :: MyFragment
  !> Batch sizes
  type(mp2_batch_construction) :: bat
  !> Calculate intgrals for first order MP2 properties?
  logical :: first_order_integrals
  !> Integrals for occ EOS: (d j|c i) in the order (d,j,c,i) [see notation inside]
  type(array4) :: goccEOS
  !> Amplitudes for occ EOS in the order (d,j,c,i) [see notation inside]
  type(array4) :: toccEOS
  !> Integrals for virt EOS: (b l|a k) in the order (b,l,a,k) [see notation inside]
  type(array4) :: gvirtEOS
  !> Amplitudes for virt EOS in the order (b,l,a,k) [see notation inside]
  type(array4) :: tvirtEOS
  !> Occ EOS integrals (d j | i k) in the order (d,j,i,k)  [only if first_order_integrals is true]
  type(array4) :: djik
  !> Virt EOS integrals (b l | a d) in the order (b,l,a,d)  [only if first_order_integrals is true]
  type(array4) :: blad


  ! Receive fragment structure and other information from master rank
  ! *****************************************************************
  call mpi_communicate_mp2_int_and_amp(MyFragment,bat,first_order_integrals,.true.)

  ! Calculate contribution to integrals/amplitudes for slave
  ! ********************************************************
  call MP2_integrals_and_amplitudes_workhorse(MyFragment,goccEOS, toccEOS, &
       & gvirtEOS, tvirtEOS, djik,blad,bat,first_order_integrals)

end subroutine MP2_integrals_and_amplitudes_workhorse_slave

#endif