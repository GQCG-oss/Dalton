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

module rif12_integrals_module 

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
  use molecule_module
  use IntegralInterfaceMod!, only: II_getBatchOrbitalInfo
  use IntegralInterfaceDEC!, only: II_precalc_DECScreenMat,&

  !  use orbital_operations

  use ccintegrals!, only: get_full_AO_integrals,get_AO_hJ,get_AO_K,get_AO_Fock

  ! Patricks mat_transpose routine 
  use reorder_frontend_module!, only: mat_transpose(rows,column,pref1,A,pref2,AT)

  ! Thomas free_cabs() for aa free MO_CABS_save_created, CMO_RI_save_created
  use CABS_operations

  ! *********************************************
  !   DEC DEPENDENCIES (within deccc directory) 
  ! *********************************************
#ifdef VAR_MPI
  use decmpi_module !, only: mpi_communicate_mp2_int_and_amp
#endif

  ! Yangs F12 routines
  use decf12_routines_module
  use f12_routines_module
  use f12ri_util_module,only: GeneralTwo4CenterDECF12RICoef1112,&
       & GeneralTwo4CenterDECF12RICoef1223,&
       & DECF12RIB4,DECF12RIB4MPI
  
  use ri_util_module
  !#endif 


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

  public :: get_rif12_fragment_energy

  private

contains

  !> Brief: Gives the single and pair fragment energy for RI-MP2-F12
  !> Author: Yang M. Wang
  !> Date: April 2015
  subroutine get_rif12_fragment_energy(MyFragment, Taibj, Tai, fragcase, Fragment1, Fragment2)
    implicit none

    !> Atomic fragment to be determined (Single or Pair fragment)
    type(decfrag), intent(inout) :: MyFragment
    !> t2EOS amplitudes stored in the order T(a,i,b,j)
    real(realk), intent(in), pointer :: Taibj(:,:,:,:) 
    !> t1EOS amplitudes stored in the order T(a,i)
    real(realk), intent(in), pointer :: Tai(:,:) 
    !> Case MODEL
    integer :: fragcase
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(in), optional :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(in), optional :: Fragment2
    !> Logical variable to check if this is a pair fragment
    logical :: dopair
    !Pairfragment
    logical,pointer :: dopair_occ(:,:)

    !========================================================
    !   Allocating integer space sizes
    !========================================================
    !> number of AO orbitals
    integer :: nbasis
    !> number of occupied MO orbitals in EOS 
    integer :: noccEOS, nvirtAOS
    !> number of occupied MO orbitals in AOS 
    integer :: noccAOS
    !> number of occupied MO orbitals in AOS tot 
    integer :: noccAOStot
    !> number of occupied MO + virtual MO orbitals in AOS 
    integer :: nocvAOS
    !> number of CABS AO orbitals
    integer :: ncabsAO
    !> number of CABS MO orbitals
    integer :: ncabsMO

    !> MO coefficient matrix for the CABS MOs
    real(realk), pointer :: CMO_CABS(:,:)
    !> MO coefficient matrix for the RI MOs
    real(realk), pointer :: CMO_RI(:,:)

    integer :: ix, iy, i, j, m, n, k, l, p, q, c, r, s, t, a, b
    real(realk) :: E_21, E_21C, E_22, E_23, E_F12
    real(realk) :: tmp, energy, tmp2
    real(realk) :: temp
    real(realk) :: mp2_energy
    real(realk) :: mp2f12_energy
    !========================================================
    !  Additional variables
    !========================================================
    real(realk),pointer :: gmo(:,:,:,:)
    real(realk),pointer :: gao(:,:,:,:)
    real(realk),pointer :: Fac(:,:)
    real(realk),pointer :: Fpp(:,:)
    real(realk) :: eps

    !========================================================
    ! Additional variables
    !========================================================
    real(realk),pointer :: CoEOS(:,:)
    real(realk),pointer :: CoAOS(:,:)
    real(realk),pointer :: CoAOStot(:,:)
    real(realk),pointer :: CvAOS(:,:)
    !> MO coefficient matrix for the OCC + VIRT
    real(realk),pointer :: Cfull(:,:)
    real(realk),pointer :: Fkj(:,:)
    real(realk),pointer :: Fab(:,:)
    real(realk),pointer :: Fnm(:,:)

    !Handling of doapir 
    integer :: noccpair
    integer,pointer :: Kval(:,:)
    !========================================================
    !  RI variables
    !========================================================
    real(realk),pointer :: CalphaR(:),CalphaG(:),CalphaF(:),CalphaD(:),CalphaCvirt(:),CalphaCocc(:)
    real(realk),pointer :: CalphaP(:), CalphaCoccT(:), CalphaX2(:)
    real(realk),pointer :: CalphaRcabsMO(:),CalphaGcabsAO(:),CalphaX(:),CalphaCcabs(:), CalphaCcabsT(:), CalphaT(:)
    real(realk),pointer :: CalphaGcabsMO(:),CalphaXcabsAO(:)
    real(realk),pointer :: ABdecompR(:,:),ABdecompG(:,:),ABdecompC(:,:)
    real(realk),pointer :: ABdecompF(:,:),Umat(:,:),Rtilde(:,:),ABdecompX(:,:),ABdecompX2(:,:)
    logical :: ABdecompCreateR,ABdecompCreateG,ABdecompCreateF,ABdecompCreateC,ABdecompCreateX2
    logical :: FORCEPRINT,use_bg_buf,LS,ABdecompCreateX

    character :: intspec(5)
    real(realk) :: TS,TE,TS2,TE2
    real(realk) :: ExchangeF12V1,CoulombF12V1
    real(realk) :: CoulombF12X1,ExchangeF12X1 
    real(realk) :: EV1,EV2,EV3,EV4,EV5,EX1,EX2,EX3,EX4
    real(realk) :: EB1,EB2,EB3,EB4,EB5,EB6,EB7,EB8,EB9
    integer :: lupri

    !> Lsitem structure
    !type(lsitem) :: mylsitem

    !> Timings
    real(realk) :: tcpu,twall
    logical :: Collaborate,DoBasis
    integer :: n1,n2,n3,n4,Tain1,Tain2

    !> bat
    type(mp2_batch_construction) :: bat

    !> Number of Auxilliary Basis functions
    integer :: NBA, NBA2
    !> MPI related variables
    integer(kind=ls_mpik) :: node
    integer :: mynum,numnodes
    logical :: wakeslaves, master
    !> Indexing
    integer :: nAtoms,naux,natomsAux,MinAuxBatch 
    !> Frozen core
    integer :: offset
    !> Special treatment of Aux
    integer(kind=long) :: nSize,nsize1,nsize2,nsize3
    logical :: ChangedDefault
    integer :: oldAOdfAux, oldAORegular, nbasis2

    !========================================================
    !  MPI related stuff
    !========================================================
    integer :: nbuf1, inode
    integer,pointer :: nAuxMPI(:)
    integer :: AuxMPIstartMy,iAuxMPIextraMy,AuxMPIstartMPI,iAuxMPIextraMPI
    integer :: numnodesstd
    real(realk),pointer :: UmatTmp(:,:)
    real(realk) :: factor, EV2tmp, EX2tmp, EV3tmp, EV4tmp
    real(realk),pointer :: CalphaMPI(:),CalphaMPI2(:)

#ifdef VAR_MPI
    real(realk) :: lsmpibufferRIMP2(20)
    lsmpibufferRIMP2=0.0E0_realk
#endif

    !mylsitem = MyFragment%MyLsitem

    lupri = DECinfo%output
#ifdef VAR_TIME
    FORCEPRINT = .TRUE.
#else
    FORCEPRINT = .FALSE.
#endif    
    call LSTIMER('START ',TS,TE,DECinfo%output,ForcePrint)

#ifdef VAR_MPI
    master= (infpar%lg_mynum == infpar%master)
    mynum = infpar%lg_mynum
    numnodes = infpar%lg_nodtot
    wakeslaves = infpar%lg_nodtot.GT.1
    if(.NOT.master)lupri = 6
#else
    ! If MPI is not used, consider the single node to be "master"
    master=.true.
    mynum = 0
    numnodes = 1
    wakeslaves = .false.
#endif
    numnodesstd = numnodes
    call LSTIMER('START ',TS2,TE2,DECinfo%output,ForcePrint)

    !========================================================
    !  Init stuff
    !========================================================
    natoms    = MyFragment%natoms
    nbasis    = MyFragment%nbasis

    nvirtAOS  = MyFragment%nvirtAOS
    noccEOS   = MyFragment%noccEOS
    noccAOS   = MyFragment%noccAOS
    
    !noccfull  = noccEOS
    ! For frozen core: noccAOS in only valence. But
    ! noccAOStot is core+valence, and nocvAOStot is core+valence+virtual,
    ! both with and without frozen core.
    noccAOStot = MyFragment%nocctot  !core + valens
    nocvAOS = noccAOStot + nvirtAOS

    ! Offset: Used for frozen core
    if(DECinfo%frozencore) then
       offset = MyFragment%ncore
    else
       offset = 0
    end if
    
    ncabsAO = size(MyFragment%Ccabs,1)    
    ncabsMO = size(MyFragment%Ccabs,2)

    nbasis2 = 0
    natomsaux = 0
    naux = 0
    E_21C = 0.0E0_realk

    !IF(DECinfo%frozencore)call lsquit('DEC RI-/MP2-F12 frozen core not implemented',-1)

    !========================================================
    !  Special treatment of Aux
    !========================================================
    call determine_maxBatchOrbitalsize(DECinfo%output,MyFragment%mylsitem%SETTING,MinAuxBatch,'D')

    IF(DECinfo%AuxAtomicExtent)THEN
       call getMolecularDimensions(MyFragment%mylsitem%INPUT%AUXMOLECULE,nAtomsAux,nBasis2,nAux)
    ELSE
       call getMolecularDimensions(MyFragment%mylsitem%SETTING%MOLECULE(1)%p,nAtomsAux,nBasis2,nAux)
       if(natoms.NE.natomsAux)call lsquit('Error in DEC RIMP2F12 natoms dim mismatch',-1)
    ENDIF

    IF(nAux.EQ.0)THEN
       WRITE(DECinfo%output,'(1X,A)')'DEC RIMP2F12 MEM: Warning no Aux basis have been chosen for RIMP2, Using Regular'
       ChangedDefault = .TRUE.
       call get_default_AOs(oldAORegular,oldAOdfAux) !the current values for Regular and Aux Basis 
       call set_default_AOs(oldAORegular,oldAORegular) !change to use Regular for Aux 
       call getMolecularDimensions(MyFragment%mylsitem%SETTING%MOLECULE(1)%p,nAtoms,nBasis2,nAux)
    ENDIF

    IF(nAux.EQ.0)THEN
       WRITE(DECinfo%output,'(1X,A)')'DEC RIMP2F12 MEM: Warning no Aux basis have been chosen for RIMP2, Using Regular'
       ChangedDefault = .TRUE.
       call get_default_AOs(oldAORegular,oldAOdfAux) !the current values for Regular and Aux Basis 
       call set_default_AOs(oldAORegular,oldAORegular) !change to use Regular for Aux 
       call getMolecularDimensions(MyFragment%mylsitem%SETTING%MOLECULE(1)%p,nAtoms,nBasis2,nAux)
    ENDIF

    !========================================================
    !  Sanity Check if we do the pair calculation
    !========================================================
    if((present(Fragment1) .AND. (.NOT. present(Fragment2))) .OR. &
         & (present(Fragment2) .AND. (.NOT. present(Fragment1)))) then
       call lsquit("get_f12_fragment_energy: Missing optional arguments Fragment1 and Fragment2")
    endif

    dopair = .FALSE.
    if(present(Fragment1) .AND. present(Fragment2)) then
       dopair = .TRUE.
    endif

    ! ***********************************************************
    !   Constructing Coefficient matrices 
    ! ***********************************************************   
    if(DECinfo%F12debug .and. master) then
       print *, "-------------------------------------------------"
       print *, "     F12-integrals.F90                           "
       print *, "-------------------------------------------------"
       print *, "nbasis:       ", nbasis
       print *, "noccEOS:      ", noccEOS
       print *, "nvirtAOS:     ", nvirtAOS
       print *, "-------------------------------------------------"
       print *, "noccAOS       ", noccAOS 
       print *, "noccAOStot    ", noccAOStot
       print *, "ncabsAO       ", ncabsAO
       print *, "ncabsMO       ", ncabsMO
       print *, "-------------------------------------------------"
    end if

     print*,'F12debug Mynum', mynum

     ! Creating a CoccEOS matrix 
    call mem_alloc(CoEOS, MyFragment%nbasis, noccEOS)
    do i=1, MyFragment%noccEOS
       ix = MyFragment%idxo(i)
       CoEOS(:,i) = MyFragment%Co(:,ix)
    end do

     print*,'CoEOS Mynum', mynum

    ! Creating a CoccAOS matrix (always core+valence, also for frozen core)
    call mem_alloc(CoAOS, MyFragment%nbasis, noccAOS)
    ! Only for frozen core: Include core MOs explicitly
    ! (without frozen core: core MOs are included in MyFragment%Co already)
    do i=1, MyFragment%noccAOS
       CoAOS(:,i) = MyFragment%Co(:,i)
    end do

    !call sleep(mynum*2)  
     print*,'CoAOS Mynum', mynum
   ! call ls_output(CoAOS,1,nbasis,1,noccAOS,nbasis,noccAOS,1,6)


    call mem_alloc(CoAOStot, MyFragment%nbasis, noccAOStot)
    do i=1,offset
       CoAOStot(:,i) = MyFragment%CoreMO(:,i)
    end do
    do i=1, MyFragment%noccAOS
       CoAOStot(:,i+offset) = MyFragment%Co(:,i)
    end do

     print*,'CvAOS Mynum', mynum

    ! Creating a CvirtAOS matrix 
    call mem_alloc(CvAOS, MyFragment%nbasis, nvirtAOS)
    do i=1, nvirtAOS
       CvAOS(:,i) = MyFragment%Cv(:,i)
    end do

     print*,'CvAOS Mynum', mynum

    ! Creating a CovAOS matrix 
    !call mem_alloc(CovAOS, MyFragment%nbasis, nocvAOS)
    !do i=1,noccAOS
    !   CovAOS(:,i) = CoAOS(:,i)
    !end do
    !do i=1,nvirtAOS
    !   CovAOS(:,i+noccAOS) = MyFragment%Cv(:,i)
    !end do

    ! ***********************************************************
    ! Creating the F matrix 
    ! ***********************************************************
    ! Creating a Fkj MO matrix occ EOS

    print *, "Creating the F matrix"
    call mem_alloc(Fkj, noccAOS, noccEOS)
    Fkj = 0E0_realk
    do j=1, noccEOS
       iy = MyFragment%idxo(j)
       Fkj(:,j) = MyFragment%ppfock(:,iy)
    end do

     print*,'Fkj Mynum', mynum

    ! Note that Fkj contains only valence orbitals since F(core,valence)=0.

    !========================================================
    !  Creating Coeff-matrices 
    !========================================================
    print *,"Creating Coeff-matrices"
    call mem_alloc(CMO_Cabs, ncabsAO, ncabsMO)
    do i=1, ncabsMO
       CMO_Cabs(:,i) = MyFragment%Ccabs(:,i)
    end do

    !call sleep(mynum*2)  
     print*,'CMO_Cabs Mynum', mynum
    !call ls_output(CoAOS,1,nbasis,1,noccAOS,nbasis,noccAOS,1,6)
    

    call mem_alloc(CMO_RI, ncabsAO, ncabsAO)
    do i=1, ncabsAO
       CMO_RI(:,i) = MyFragment%Cri(:,i)
    end do


    !call sleep(mynum*2)  
     print*,'CMO_RI Mynum', mynum
    !call ls_output(CMO_RI,1,ncabsAO,1,ncabsAO,ncabsAO,ncabsAO,1,6)


    ! Creating a CocvAOStot matrix 
    call mem_alloc(Cfull, MyFragment%nbasis, nocvAOS)
    do i=1,noccAOStot                                                                                                        
       Cfull(:,i) = CoAOStot(:,i)                                                                                     
    end do                                                                                                                   
    do i=1,nvirtAOS
       Cfull(:,i+noccAOStot) = MyFragment%Cv(:,i)                                                                       
    end do   
   
     print*,'Cfull Mynum', mynum
 
    call mem_alloc(Fac, nvirtAOS, ncabsMO)
    do a=1, nvirtAOS
       do c=1, ncabsMO
          Fac(a,c) = Myfragment%Fcp(c,a+noccAOStot)
       enddo
    enddo

     print*,'Fac Mynum', mynum

    !Slow implementation can be done with 2dgemms in B6
    call mem_alloc(Fpp,nvirtAOS+noccAOStot,nvirtAOS+noccAOStot)
    Fpp = 0.0E0_realk
    do p=1,offset
       do q=1,offset
          Fpp(p,q) = MyFragment%ccfock(p,q)
       enddo
    enddo

    do p=1+offset, noccAOStot
       do q=1+offset, noccAOStot
          Fpp(p,q) = MyFragment%ppfock(p-offset,q-offset)
       enddo
    enddo

    do p=noccAOStot+1, nvirtAOS+noccAOStot
       do q=noccAOStot+1, nvirtAOS+noccAOStot
          Fpp(p,q) = MyFragment%qqfock(p-noccAOStot,q-noccAOStot)
       enddo
     enddo

     print*,'Fpp Mynum', mynum

    call mem_alloc(Fnm,noccAOStot,noccAOStot)
    Fnm = 0.0E0_realk
    do n=1,offset
       do m=1,offset
          Fnm(n,m) = MyFragment%ccfock(n,m)
       enddo
    enddo
    do n=1+offset, noccAOStot
       do m=1+offset, noccAOStot
          Fnm(n,m) = MyFragment%ppfock(n-offset,m-offset)
       enddo
    enddo

    !call sleep(mynum*2)  
    print*,'Fnm Mynum', mynum
    !call ls_output(Fnm,1,noccAOStot,1,noccAOStot,noccAOStot,noccAOStot,1,6)


  !   Fpp = 0.0E0_realk

  !    do p=1, noccAOStot
  !      do q=1, noccAOStot
  !         Fpp(p,q) = MyFragment%ppfock(p,q)
  !     enddo
  !  enddo
  !  do p=noccAOStot+1, nvirtAOS+noccAOStot
  !     do q=noccAOStot+1, nvirtAOS+noccAOStot
  !        Fpp(p,q) = MyFragment%qqfock(p-noccAOStot,q-noccAOStot)
  !     enddo
  !  enddo

  !=================================================================
  != Step 0: Creating of dopair_occ                                =
  !=================================================================
#ifdef VAR_MPI

   print *, "Before waking up the slaves: ", mynum

  !=================================================================
  != Step 1: Wake up the slaves                                    =
  !=================================================================
  if(wakeslaves .and. master) then
     print *, "We wake up the slaves"
     call ls_mpibcast(DECRIMP2F12,infpar%master,infpar%lg_comm)

     call mpi_communicate_MyFragment_f12(MyFragment,Taibj,fragcase,dopair)

     call mem_alloc(dopair_occ,noccEOS,noccEOS)
     dopair_occ = .FALSE.

     print *,"master dopair", dopair
     if(dopair) then 
        call which_pairs_occ(Fragment1,Fragment2,MyFragment,dopair_occ)
        call ls_mpibcast(dopair_occ,noccEOS,noccEOS,infpar%master,infpar%lg_comm)
     else
        dopair_occ = .TRUE.
     endif

  else
     call mem_alloc(dopair_occ,noccEOS,noccEOS)
     dopair_occ = .FALSE.
      print *,"slave dopair", dopair
     if(dopair) then 
        call ls_mpibcast(dopair_occ,noccEOS,noccEOS,infpar%master,infpar%lg_comm)
        !call which_pairs_occ(Fragment1,Fragment2,MyFragment,dopair_occ)
     else
        dopair_occ = .TRUE.
     endif 
  endif

#else

  call mem_alloc(dopair_occ,noccEOS,noccEOS)
  dopair_occ = .FALSE.
  if(dopair) then 
     call which_pairs_occ(Fragment1,Fragment2,MyFragment,dopair_occ)
  else
     dopair_occ = .TRUE.
  endif

#endif

    !=================================================================
    != Step 2:  Fijkl,Xijkl,Dijkl                                    =
    !=          corresponding to V1,X1,B1                            =
    !=================================================================
    !normally I do not like to allocate things at the beginning but 
    !due to an analysis of memory heap performance and the 
    !background buffer features of ordered allocations this is beneficial
    call mem_alloc(ABdecompR,nAux,nAux)
    call mem_alloc(ABdecompF,nAux,nAux)
    call mem_alloc(ABdecompX,nAux,nAux)
    ABdecompCreateF = .TRUE.
    ABdecompCreateX = .TRUE.
    ABdecompCreateR = .TRUE.
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'R' !Regular AO basis function on center 4

    ! Calculate the Fitting Coefficients (alpha|F|ij)
    use_bg_buf = .FALSE.
    mp2f12_energy = 0.0E0_realk 
    intspec(4) = 'F' !The Gaussian geminal divided by the Coulomb operator g/r12 (GGemCouOperator)
    intspec(5) = 'F' !The Gaussian geminal divided by the Coulomb operator g/r12 (GGemCouOperator)
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CoEOS,noccEOS,&
         & mynum,numnodes,CalphaF,NBA,ABdecompF,ABdecompCreateF,intspec,use_bg_buf)
    ABdecompCreateF = .FALSE.
    !perform this suborutine on the GPU (async)  - you do not need to wait for the results
    call ContractOne4CenterF12IntegralsRIV1_dec(NBA,noccEOS,CalphaF,CoulombF12V1,ExchangeF12V1,dopair_occ)
    EV1 = -1.0E0_realk*((5.0E0_realk*0.25E0_realk)*CoulombF12V1-ExchangeF12V1*0.25E0_realk)

#ifdef VAR_MPI 
    print *, "EV1, mynum", EV1, mynum
    lsmpibufferRIMP2(2)=EV1
#else
    mp2f12_energy = mp2f12_energy  + EV1
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V1,RI) = ', EV1
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V1,RI) = ', EV1
#endif         

    !Calculate the Fitting Coefficients (alpha|g^2|ij) 
    intspec(4) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
    intspec(5) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CoAOS,noccAOS,&
         & mynum,numnodes,CalphaX,NBA,ABdecompX,ABdecompCreateX,intspec,use_bg_buf)
    ABdecompCreateX = .FALSE.

    !perform this suborutine on the GPU (async)  - you do not need to wait for the results
    nsize = NBA*noccEOS*noccEOS
    call mem_alloc(CalphaT,nsize)     ! G_ij = C_ik F_kj 
    M = NBA*noccEOS         !rows of Output Matrix
    K = noccAOS             !summation dimension
    N = noccEOS             !columns of Output Matrix
    !call dgemm('N','N',M,N,K,1.0E0_realk,CalphaX,M,Fmm%elms,K,0.0E0_realk,CalphaT,M)
    call dgemm('N','N',M,N,K,1.0E0_realk,CalphaX,M,Fkj,K,0.0E0_realk,CalphaT,M)
    call ContractOne4CenterF12IntegralsRIX1_dec(NBA,noccEOS,CalphaX,CalphaT,EX1,dopair_occ)
    call mem_dealloc(CalphaT)

#ifdef VAR_MPI 
    print *, "EX1, mynum", EX1, mynum
    lsmpibufferRIMP2(3)=EX1
#else
    mp2f12_energy = mp2f12_energy  + EX1
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X1,RI) = ', EX1
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X1,RI) = ', EX1
#endif         

    !Calculate the Fitting Coefficients (alpha|[[T,g],g]|ij) 
    intspec(4) = 'D' !The double commutator [[T,g],g] 
    intspec(5) = 'C' !The metric operator = Regular Coulomb operator 1/r12
    !Build the R coefficient of Eq. 90 of J Comput Chem 32: 2492–2513, 2011
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CoEOS,noccEOS,&
         & mynum,numnodes,CalphaD,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
    ABdecompCreateR = .FALSE.
    !We need CalphaR(NBA,noccEOS,noccEOS) but this is a subset of the CalphaR(NBA,noccEOS,nbasis)
    !so we calculate the full CalphaR(NBA,noccEOS,nbasis)
    intspec(4) = 'C' !Regular Coulomb operator 1/r12
    intspec(5) = 'C' !The metric operator = Regular Coulomb operator 1/r12
    !Build the G coefficient of Eq. 90 of J Comput Chem 32: 2492–2513, 2011
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,Cfull,nocvAOS,&
         & mynum,numnodes,CalphaR,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
    !Build the U matrix in Eq. 88 of J Comput Chem 32: 2492–2513, 2011
    call mem_alloc(Umat,nAux,nAux)
    !perform this suborutine on the GPU (Async)
    call Build_RobustERImatU(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CoEOS,noccEOS,&
         & mynum,numnodes,ABdecompR,'D',Umat)
  
#ifdef VAR_MPI
      !Build the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 2011
      !We need to do:
      !CalphaD(NBA,nocc,nocc) = -0.5*CalphaD(NBA,nocc,nocc) + Loop_NBA2 Umat(NBA,NBA2)*CalphaR(NBA2,nocc,nocc)
      !Where NBA is the Auxiliary assigned to this node, while NBA2 can be assigned to another node
      IF(wakeslaves)THEN
         nbuf1=numnodes
         call mem_alloc(nAuxMPI,nbuf1)
         call BuildnAuxMPIUsedRI(nAux,numnodesstd,nAuxMPI) !3 50     
         call BuildnAuxMPIUsedRIinfo(nAux,numnodesstd,mynum,AuxMPIstartMy,iAuxMPIextraMy)
         DO inode = 1,numnodes
            nbuf1 = nAuxMPI(inode)
            NBA2 = nAuxMPI(inode)
            nsize = nbuf1*noccEOS*noccEOS
            !nbuf1 = NBA2
            IF(inode.EQ.1) THEN
               factor = -0.5E0_realk
            ELSE
               factor = 1.0E0_realk
            ENDIF
            call BuildnAuxMPIUsedRIinfo(nAux,numnodesstd,inode-1,AuxMPIstartMPI,iAuxMPIextraMPI)
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(UmatTmp,NBA*i8,NBA2*i8)
            ELSE
               call mem_alloc(UmatTmp,NBA,NBA2)
            ENDIF
            Call BuilDUmatTmpRIF12(Umat,nAux,UmatTmp,NBA,NBA2,AuxMPIstartMy,iAuxMPIextraMy,&
                 & AuxMPIstartMPI,iAuxMPIextraMPI)
            IF(mynum.EQ.inode-1)THEN
               !I Bcast My Own CalphaR (the first NBA*nocc*nocc part )
               node = mynum
               call ls_mpibcast(CalphaR,nsize,node,infpar%lg_comm)
               M = NBA                !rows of Output Matrix
               N = noccEOS*noccEOS    !columns of Output Matrix
               K = NBA                !summation dimension
               call dgemm('N','N',M,N,K,1.0E0_realk,UmatTmp,M,CalphaR,K,factor,CalphaD,M)
            ELSE
               node = inode-1
               !recieve
               IF(use_bg_buf)THEN
                  call mem_pseudo_alloc(CalphaMPI,nsize)
               ELSE
                  call mem_alloc(CalphaMPI,nsize)
               ENDIF
               call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
               M = NBA                !rows of Output Matrix
               N = noccEOS*noccEOS    !columns of Output Matrix
               K = NBA2               !summation dimension
               call dgemm('N','N',M,N,K,1.0E0_realk,UmatTmp,M,CalphaMPI,K,factor,CalphaD,M)
               IF(use_bg_buf)THEN
                  call mem_pseudo_dealloc(CalphaMPI)
               ELSE
                  call mem_dealloc(CalphaMPI)
               ENDIF
            ENDIF
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(UmatTmp)
            ELSE
               call mem_dealloc(UmatTmp)
            ENDIF
         ENDDO
      ELSE

         M = NBA          !rows of Output Matrix
         K = NBA          !summation dimension
         N = noccEOS*noccEOS    !columns of Output Matrix
         call dgemm('N','N',M,N,K,1.0E0_realk,Umat,M,CalphaR(1+offset*NBA*noccEOS),K,-0.5E0_realk,CalphaD,M)

      ENDIF
#else
      !Build the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 2011
      !perform this suborutine on the GPU (Async)
      !note CalphaR is actual of dimensions (NBA,nocc,nbasis) but here we only access
      !the first part (NBA,nocc,nocc) 
      M = NBA          !rows of Output Matrix
      K = NBA          !summation dimension
      N = noccEOS*noccEOS    !columns of Output Matrix
      call dgemm('N','N',M,N,K,1.0E0_realk,Umat,M,CalphaR(1+offset*NBA*noccEOS),K,-0.5E0_realk,CalphaD,M)
#endif 
    call mem_dealloc(Umat)
    !CalphaD is now the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 20
    !perform this suborutine on the GPU (Async)
    call ContractOne4CenterF12IntegralsRobustRIB1_dec(nBA,offset,noccEOS,nocvAOS,CalphaD,CalphaR,EB1,dopair_occ)

#ifdef VAR_MPI 
    print *, "EB1, mynum", EB1, mynum
    lsmpibufferRIMP2(1)=EB1
#else
    mp2f12_energy = mp2f12_energy  + EB1
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B1,RI) = ', EB1
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B1,RI) = ', EB1
#endif         

    !The minus is due to the Valeev factor
    !mp2f12_energy = mp2f12_energy + EV1
    !WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V1,RI) = ',EV1
    !WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V1,RI) = ',EV1
    !minus is due to the overall minus from equation (41) and (42) due to
    !contribution from the \bar{B}_{ij}^{ij}
    !EX1 = -1.0E0_realk*(0.21875E0_realk*CoulombF12X1 + 0.03125E0_realk*ExchangeF12X1)
    !mp2f12_energy = mp2f12_energy + EX1
    !WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X1,RI) = ',EX1
    !WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X1,RI) = ',EX1
    !mp2f12_energy = mp2f12_energy + EB1
    !WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B1,RI) = ', EB1
    !WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B1,RI) = ', EB1

    call mem_dealloc(CalphaD)
    call LSTIMER('FULLRIMP2:Step1',TS2,TE2,DECinfo%output,ForcePrint)

    !==============================================================
    !=  B2: sum_c' (ic'|f12^2|jj) hJ_ic' - (jc'|f12^2|ij) hJ_ic'  =
    !=  B3: sum_c' (ii|f12^2|jc') hJ_jc' - (ij|f12^2|ic') hJ_ic'  =
    !==============================================================
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'C' !CABS AO basis function on center 4
    intspec(4) = '2' !The f12 Operator
    intspec(5) = '2' !The f12 Operator 
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CMO_RI,ncabsAO,&
         & mynum,numnodes,CalphaXcabsAO,NBA,ABdecompX,ABdecompCreateX,intspec,use_bg_buf)

    call ContractOne4CenterF12IntegralsRIB23_dec(NBA,noccEOS,ncabsAO,CalphaXcabsAO,CalphaX,&
         & MyFragment%hJir,1.0E0_realk,EB2,EB3,dopair_occ)

#ifdef VAR_MPI 
print *, "EB2, mynum", EB2, mynum
print *, "EB3, mynum", EB3, mynum
   lsmpibufferRIMP2(4)=EB2       !we need to perform a MPI reduction at the end 
   lsmpibufferRIMP2(5)=EB3       !we need to perform a MPI reduction at the end 
#else
   !1.0E0_realk because that term has an overall pluss in Eqs. 25-26
   mp2f12_energy = mp2f12_energy  + EB2
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B2,RI) = ',EB2
   mp2f12_energy = mp2f12_energy  + EB3
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B3,RI) = ',EB3
   WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B2,RI) = ',EB2
   WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B3,RI) = ',EB3
#endif

    call mem_dealloc(CalphaXcabsAO)
    call mem_dealloc(CalphaX)
    call mem_dealloc(CalphaF)
    call mem_dealloc(ABdecompX)
    call mem_dealloc(ABdecompF)
    !ABdecompR still exists
    call LSTIMER('FULLRIMP2:B2B3',TS2,TE2,DECinfo%output,ForcePrint)

    !=================================================================
    != Step 2: Ripjq*Gipjq+Rimjc*Gimjc+Rjmic*Gjmic                   =
    !=        +Gipjq*Gipjq+Gimjc*Gimjc+Gjmic*Gjmic                   =
    !=         corresponding to V2,V3,V4,X2,X3,X4                    =
    != These are special since                                       =
    != 1. the have (more or less) the same structure                 =
    != 2. They can be built from same intermediates                  =
    !=    CalphaR(alpha,i,p),CalphaR(alpha,i,c)                      =
    !=    CalphaG(alpha,i,p),CalphaG(alpha,i,c)                      =
    != 3. The intermediates are only used once                       =
    != 4. Due to noccEOS,noccEOS,noccEOS,noccEOS very small mem req  =
    !=================================================================

    call mem_alloc(ABdecompG,nAux,nAux)
    ABdecompCreateG = .TRUE.

    !==========================================================
    !=                                                        =
    != V2: Ripjq*Gipjq                                        =
    != The Coulomb Operator Int multiplied with               =
    != The Gaussian geminal operator g                        =
    != Dim(noccEOS,nbasis,noccEOS,nbasis)                     =
    !=                                                        =
    !==========================================================
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'R' !Regular AO basis function on center 4
    intspec(4) = 'G' !The Gaussian geminal operator g
    intspec(5) = 'G' !The Gaussian geminal operator g
    !Unique: CalphaG(NBA,nocvAOS,noccEOS)
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI, &
       & FORCEPRINT,wakeslaves,Cfull,nocvAOS,CoEOS,noccEOS, &
       & mynum,numnodes,CalphaG,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
    ABdecompCreateG = .FALSE.

    !X2
    !==========================================================
    !=                                                        =
    != X2: Gipjq*Gipjq                                        =
    != The Gaussian geminal operator Int multiplied with      =
    != The Gaussian geminal operator g                        =
    != Dim(nocc,nbasis,nocc,nbasis)                           =
    !=                                                        =
    !==========================================================
    !NB No we need this term since its 
    call mem_alloc(ABdecompX2,nAux,nAux)
    ABdecompCreateX2 = .TRUE.
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI, &
       & FORCEPRINT,wakeslaves,Cfull,nocvAOS,CoAOS,noccAOS, &
       & mynum,numnodes,CalphaX2,NBA,ABdecompX2,ABdecompCreateX2,intspec,use_bg_buf)
    call mem_dealloc(ABdecompX2)

    !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
#ifdef VAR_MPI
    !X2 
    nsize = NBA*nocvAOS*noccEOS
    call mem_alloc(CalphaT,nsize)  !G_qj = C_qk F_kj 
    M = nocvAOS*NBA                !rows of Output Matrix
    K = noccAOS                    !summation dimension
    N = noccEOS                    !columns of Output Matrix
    call dgemm('N','N',M,N,K,1.0E0_realk,CalphaX2,M,Fkj,K,0.0E0_realk,CalphaT,M)

    IF(wakeslaves)THEN
       EV2 = 0.0E0_realk
       EX2 = 0.0E0_realk
       DO inode = 1,numnodes
          nbuf1 = nAuxMPI(inode)
          NBA2 =  nAuxMPI(inode)
          nsize = nbuf1*nocvAOS*noccEOS
          IF(mynum.EQ.inode-1)THEN
             !I Bcast My Own CalphaG
             node = mynum
             IF(size(CalphaG).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 A',-1)
             call ls_mpibcast(CalphaG,nsize,node,infpar%lg_comm)
             call ContractTwo4CenterF12IntegralsRIV2_dec(nBA,nBA2,noccEOS,nocvAOS,CalphaR,CalphaG,EV2tmp,dopair_occ)
             call ContractTwo4CenterF12IntegralsRIX2_dec(nBA,nBA2,noccEOS,nocvAOS,CalphaG,CalphaG,CalphaT,EX2tmp,dopair_occ)  
          ELSE
             node = inode-1
             !recieve
             IF(use_bg_buf)THEN
                call mem_pseudo_alloc(CalphaMPI,nsize)
             ELSE
                call mem_alloc(CalphaMPI,nsize)
             ENDIF
             call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
             call ContractTwo4CenterF12IntegralsRIV2_dec(nBA,nBA2,noccEOS,nocvAOS,CalphaR,CalphaMPI,EV2tmp,dopair_occ)
             call ContractTwo4CenterF12IntegralsRIX2_dec(nBA,nBA2,noccEOS,nocvAOS,CalphaG,CalphaMPI,CalphaT,EX2tmp,dopair_occ)
             IF(use_bg_buf)THEN
                call mem_pseudo_dealloc(CalphaMPI)
             ELSE
                call mem_dealloc(CalphaMPI)
             ENDIF
          ENDIF
          EV2 = EV2 + EV2tmp
          EX2 = EX2 + EX2tmp
       ENDDO
    ELSE
       call ContractTwo4CenterF12IntegralsRIV2_dec(nBA,nBA,noccEOS,nocvAOS,CalphaR,CalphaG,EV2,dopair_occ)
       call ContractTwo4CenterF12IntegralsRIX2_dec(nBA,nBA,noccEOS,nocvAOS,CalphaG,CalphaG,CalphaT,EX2,dopair_occ)
    ENDIF
    
    print *, "EV2, mynum", EV2, mynum
    print *, "EX2, mynum", EX2, mynum
    lsmpibufferRIMP2(6)=EV2      !we need to perform a MPI reduction at the end 
    lsmpibufferRIMP2(7)=EX2      !we need to perform a MPI reduction at the end 
    call mem_dealloc(CalphaT)

#else

    call ContractTwo4CenterF12IntegralsRIV2_dec(nBA,nBA,noccEOS,nocvAOS,CalphaR,CalphaG,EV2,dopair_occ)
    mp2f12_energy = mp2f12_energy + EV2
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V2,RI) = ',EV2
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V2,RI) = ',EV2


    !X2 
    nsize = NBA*nocvAOS*noccEOS
    call mem_alloc(CalphaT,nsize)  !G_qj = C_qk F_kj 
    M = nocvAOS*NBA                !rows of Output Matrix
    K = noccAOS                    !summation dimension
    N = noccEOS                    !columns of Output Matrix
    call dgemm('N','N',M,N,K,1.0E0_realk,CalphaX2,M,Fkj,K,0.0E0_realk,CalphaT,M)

    call mem_dealloc(CalphaX2)
    call ContractTwo4CenterF12IntegralsRIX2_dec(nBA,nBA,noccEOS,nocvAOS,CalphaG,CalphaG,CalphaT,EX2,dopair_occ)
    call mem_dealloc(CalphaT)

    mp2f12_energy = mp2f12_energy  + EX2
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X2,RI) = ',EX2
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X2,RI) = ',EX2

#endif

    !==========================================================
    !=                                                        =
    != V3: Rimjc*Gimjc                                        =
    != V4: Rjmic*Gjmic                                        =
    != The Coulomb Operator Int multiplied with               =
    != The Gaussian geminal operator g                        =
    != Dim: (noccEOS,noccAOS,noccEOS,ncabsMO)  need 4 Calphas =
    !=                                                        =
    !==========================================================

    !   We need CalphaRocc(NBA,noccEOS,noccEOS) but this is a subset of CalphaR(NBA,noccEOS,nbasis)
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'C' !CABS AO basis function on center 4
    intspec(4) = 'C' !The Coulomb Operator
    intspec(5) = 'C' !The Coulomb Operator
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CMO_CABS,ncabsMO,&
         & mynum,numnodes,CalphaRcabsMO,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)

    !   We need CalphaGocc(NBA,noccEOS,noccEOS) but this is a subset of CalphaG(NBA,noccEOS,nbasis)
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'C' !CABS AO basis function on center 4
    intspec(4) = 'G' !The Gaussian geminal operator g
    intspec(5) = 'G' !The Gaussian geminal operator g
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CMO_CABS,ncabsMO,&
         & mynum,numnodes,CalphaGcabsMO,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)

#ifdef VAR_MPI 
   IF(wakeslaves)THEN
      EV3 = 0.0E0_realk
      EV4 = 0.0E0_realk
      DO inode = 1,numnodes
         nbuf1 = nAuxMPI(inode)
         NBA2 = nAuxMPI(inode)
         nsize = nbuf1*noccEOS*ncabsMO  !CalphaRcabsMO(NBA,nocc,ncabsMO)
         nsize2 = nbuf1*noccEOS*nocvAOS    !CalphaR(NBA,nocc,nocv)
         IF(mynum.EQ.inode-1)THEN
            !I Bcast My Own CalphaG
            node = mynum
            IF(size(CalphaRcabsMO).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 C1',-1)
            call ls_mpibcast(CalphaRcabsMO,nsize,node,infpar%lg_comm)
            IF(size(CalphaR).NE.nsize2)call lsquit('MPI Bcast error in Full RIMP2F12 C2',-1)
            call ls_mpibcast(CalphaR,nsize2,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIV34_dec(NBA,NBA,noccEOS,noccAOS,ncabsMO,nocvAOS,&
                 & CalphaRcabsMO,CalphaGcabsMO,CalphaR,CalphaG,EV3tmp,EV4tmp)
         ELSE
            node = inode-1
            !recieve
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(CalphaMPI,nsize)
               call mem_pseudo_alloc(CalphaMPI2,nsize2)
            ELSE
               call mem_alloc(CalphaMPI,nsize)
               call mem_alloc(CalphaMPI2,nsize2)
            ENDIF
            call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
            call ls_mpibcast(CalphaMPI2,nsize2,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIV34_dec(NBA,NBA2,noccEOS,noccAOS,ncabsMO,nocvAOS,&
                 & CalphaMPI,CalphaGcabsMO,CalphaMPI2,CalphaG,EV3tmp,EV4tmp)
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(CalphaMPI2)
               call mem_pseudo_dealloc(CalphaMPI)
            ELSE
               call mem_dealloc(CalphaMPI)
               call mem_dealloc(CalphaMPI2)
            ENDIF
         ENDIF
         EV3 = EV3 + EV3tmp
         EV4 = EV4 + EV4tmp
      ENDDO
   ELSE
      call ContractTwo4CenterF12IntegralsRIV34_dec(NBA,NBA,noccEOS,noccAOS,ncabsMO,nocvAOS,&
           & CalphaRcabsMO,CalphaGcabsMO,CalphaR,CalphaG,EV3,EV4)
   ENDIF

    print *, "EV3, mynum", EV3, mynum
    print *, "EV4, mynum", EV4, mynum

   lsmpibufferRIMP2(8)=EV3      !we need to perform a MPI reduction at the end 
   lsmpibufferRIMP2(9)=EV4      !we need to perform a MPI reduction at the end 
#else
   !Do on GPU (Async)
   call ContractTwo4CenterF12IntegralsRIV34_dec(NBA,NBA,noccEOS,noccAOS,ncabsMO,nocvAOS,&
        & CalphaRcabsMO,CalphaGcabsMO,CalphaR,CalphaG,EV3,EV4)
   mp2f12_energy = mp2f12_energy  + EV3 + EV4
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V3,RI) = ',EV3
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V4,RI) = ',EV4
   WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V3,RI) = ',EV3
   WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V4,RI) = ',EV4
#endif
   call mem_dealloc(ABdecompR)
   call mem_dealloc(CalphaR)
   call mem_dealloc(CalphaRcabsMO)

    !==========================================================
    != V5: Caibj = (Gcibj*Fac + Gcjai*Fcb)*Taibj              =
    !========================================================== 
    call mem_alloc(ABdecompC,nAux,nAux)
    ABdecompCreateC = .TRUE.
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 4 
    intspec(3) = 'C' !Cabs AO basis function on center 3
    intspec(4) = 'G' !The Gaussian geminal operator g
    intspec(5) = 'G' !The Gaussian geminal operator g

    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CMO_CABS,ncabsMO,&
         & mynum,numnodes,CalphaCcabs,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)

    m = NBA*noccEOS       
    k = ncabsMO         ! C_mn = A_mk B_kn
    n = nvirtAOS

    !C(alpha*i,ncabsMO)*F(cabsMO,nvirt)
    nsize = nBA*noccEOS*ncabsMO
    call mem_alloc(CalphaD, nsize)
    call dgemm('N','T',m,n,k,1.0E0_realk,CalphaCcabs,m,Fac,n,0.0E0_realk,CalphaD,m)

    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'R' !Regular AO basis function on center 4
    intspec(4) = 'G' !The Gaussian geminal operator g
    intspec(5) = 'G' !The Gaussian geminal operator g
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CvAOS,nvirtAOS,&
         & mynum,numnodes,CalphaCvirt,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)

    call ContractTwo4CenterF12IntegralsRIV5_dec(nBA,noccEOS,nvirtAOS,CalphaCvirt,CalphaD,Taibj,EV5,dopair_occ)


#ifdef VAR_MPI 
    print *, "EV5, mynum", EV5, mynum
    lsmpibufferRIMP2(10)=EV5
#else
    mp2f12_energy = mp2f12_energy  + EV5
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V5,RI) = ', EV5
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V5,RI) = ', EV5
#endif  


    ABdecompCreateG = .FALSE.
    call mem_dealloc(CalphaD)
    call mem_dealloc(ABdecompC)
    call mem_dealloc(CalphaCcabs)
    call mem_dealloc(CalphaCvirt)

    !==========================================================
    !=                                                        =
    != X3:         Step 3  Gimjc*Gimjc                        =
    != X4:         Step 4  Gjmic*Gjmic                        =
    != The Coulomb Operator Int multiplied with               =
    != The Gaussian geminal operator g                        =
    != Dim: (noccEOS,noccAOS,noccEOS,ncabsMO)  need 4 Calphas =
    !=                                                        =
    !==========================================================
    !   We need CalphaGocc(NBA,noccEOS,noccEOS) but this is a subset of CalphaG(NBA,noccEOS,nbasis)
     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 3
     intspec(3) = 'R' !Regular AO basis function on center 4
     intspec(4) = 'G' !The Gaussian geminal operator g
     intspec(5) = 'G' !The Gaussian geminal operator g
     
     call mem_alloc(ABdecompC,nAux,nAux)
     call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,CoAOStot,noccAOStot,CoAOS,noccAOS,mynum,numnodes,CalphaCocc,&
        & NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)
     call mem_dealloc(ABdecompC)

     !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
     nsize = NBA*noccAOStot*noccEOS
     call mem_alloc(CalphaT,nsize)
     M = NBA*noccAOStot        
     K = noccAOS            !G_mj = C_mk F_kj  
     N = noccEOS             
     call dgemm('N','N',M,N,K,1.0E0_realk,CalphaCocc,M,Fkj,K,0.0E0_realk,CalphaT,M)

     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'C' !Regular AO basis function on center 4 
     intspec(3) = 'R' !Cabs AO basis function on center 3
     intspec(4) = 'G' !The Gaussian geminal operator g
     intspec(5) = 'G' !The Gaussian geminal operator g
     
     call mem_alloc(ABdecompC,nAux,nAux)
     call Build_CalphaMO2(MyFragment%MyLsitem,master,ncabsAO,nbasis,nAux,LUPRI, &
        & FORCEPRINT,wakeslaves,CMO_CABS,ncabsMO,CoAOS,noccAOS, &
        & mynum,numnodes,CalphaCcabsT,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)
     call mem_dealloc(ABdecompC)
     
     nsize = NBA*noccEOS*ncabsMO
     call mem_alloc(CalphaP,nsize)
     M = NBA*ncabsMO        
     K = noccAOS               !P_cj = C_ck F_kj   
     N = noccEOS          
     call dgemm('N','N',M,N,K,1.0E0_realk,CalphaCcabsT,M,Fkj,K,0.0E0_realk,CalphaP,M)
     !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
     call ContractTwo4CenterF12IntegralsRIX34_dec(NBA,noccEOS,noccAOStot,ncabsMO,&
        & CalphaGcabsMO,CalphaCocc,CalphaT,CalphaP,EX3,EX4,dopair_occ)

     call mem_dealloc(CalphaCcabsT)
     call mem_dealloc(CalphaCocc)
     call mem_dealloc(CalphaT)
     call mem_dealloc(CalphaP)

     mp2f12_energy = mp2f12_energy  + EX3 + EX4
     WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X3,RI) = ',EX3
     WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X4,RI) = ',EX4
     WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X3,RI) = ',EX3
     WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X4,RI) = ',EX4

   call LSTIMER('FULLRIMP2:Step2',TS2,TE2,DECinfo%output,ForcePrint)

    !=================================================================
    != Step 3: The remaining B terms                                 =
    !=================================================================

    !==============================================================
    !=  B4: (ir|f12|jt)Kst(ir|f12|js)      (r,s,t=CabsAO)         =
    !==============================================================
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'C' !Regular AO basis function on center 4
    intspec(4) = 'G' !The f12 Operator
    intspec(5) = 'G' !The f12 Operator
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CMO_RI,ncabsAO,&
         & mynum,numnodes,CalphaGcabsAO,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)

    nsize = nBA*noccEOS*ncabsAO
    call mem_alloc(CalphaD, nsize)
    m =  nBA*noccEOS                    ! D_js = C_jr K_rs
    k =  ncabsAO
    n =  ncabsAO

    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,MyFragment%Krs,k,0.0E0_realk,CalphaD,m)
    call BuildKval(dopair_occ,noccEOS,Kval,noccpair,dopair)
    call GeneralTwo4CenterDECF12RICoef1112(nBA,CalphaGcabsAO,noccEOS,ncabsAO,CalphaD,noccEOS,ncabsAO,EB4,&
         & noccEOS,wakeslaves,use_bg_buf,numnodes,nAuxMPI,mynum,Kval,noccpair,DECF12RIB4,DECF12RIB4MPI)
    call mem_dealloc(KVAL)

    mp2f12_energy = mp2f12_energy + EB4
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B4,RI) = ',EB4
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B4,RI) = ', EB4

    !==============================================================
    !=  B5: (ir|f12|jm)Frs(si|f12|mj)        (r,s=CabsAO)         =
    !==============================================================   
    !   We need CalphaGocc(NBA,noccEOS,noccEOS) but this is a subset of CalphaG(NBA,noccEOS,nbasis)
     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 3
     intspec(3) = 'R' !Regular AO basis function on center 4
     intspec(4) = 'G' !The Gaussian geminal operator g
     intspec(5) = 'G' !The Gaussian geminal operator g
    
     !NB! Can be optimized! 
     call mem_alloc(ABdecompC,nAux,nAux)
     call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,CoAOStot,noccAOStot,CoEOS,noccEOS,mynum,numnodes,CalphaCocc,&
        & NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)
     call mem_dealloc(ABdecompC)

   !We need CalphaG(NBA,nocc,noccAOS) but this is a subset of                                                               
   !CalphaG(NBA,nocc,nbasis) which we already have                                                                           
   !> Dgemm 
   nsize = nBA*noccEOS*ncabsAO                                                                                                  
   IF(size(CalphaD).NE.nsize)call lsquit('dim mismatch CalphaD',-1)                                                          
   m =  nBA*noccEOS                 
   k =  ncabsAO               ! D_js = C_jr F_rs                                                                       
   n =  ncabsAO
               
   !Do on GPU (Async)
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,MyFragment%Frs,k,0.0E0_realk,CalphaD,m)                                    
   !Do on GPU (Async)
   call ContractTwo4CenterF12IntegralsRIB5_dec(nBA,noccEOS,ncabsAO,noccAOStot,nocvAOS,CalphaGcabsAO,CalphaCocc, &
                                              & CalphaD,EB5,dopair_occ)         
               
   mp2f12_energy = mp2f12_energy  + EB5
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B5,RI) = ',EB5                                       
   WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B5,RI) = ', EB5                                                   
   call mem_dealloc(CalphaD)              

    !==============================================================
    !=  B6: (ip|f12|ja)Fqp(qi|f12|aj)                             =
    !==============================================================
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 1 
    intspec(3) = 'R' !Regular AO basis function on center 2 
    intspec(4) = 'G' !The Gaussian geminal operator g
    intspec(5) = 'G' !The Gaussian geminal operator g

    nsize = NBA*noccEOS*nocvAOS
    ! call mem_alloc(ABdecompG,nAux,nAux)
    call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
       & FORCEPRINT,wakeslaves,CoEOS,noccEOS,Cfull,nocvAOS,&
       & mynum,numnodes,CalphaP,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
    !call mem_dealloc(ABdecompG)

    nsize = nBA*noccEOS*nocvAOS
    call mem_alloc(CalphaD,nsize)

    m =  nBA*noccEOS               ! D_jq =P_jp F_pq
    k =  nocvAOS
    n =  nocvAOS
    !Do on GPU (Async)
    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaP,m,Fpp,k,0.0E0_realk,CalphaD,m)
    call mem_dealloc(CalphaP)
    !Do on GPU (Async)
    call ContractTwo4CenterF12IntegralsRIB6_dec(nBA,noccEOS,nvirtAOS,nocvAOS,noccAOStot,&
                                              & CalphaG,CalphaD,EB6,dopair_occ)
    call mem_dealloc(CalphaD)

    mp2f12_energy = mp2f12_energy  + EB6
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B6,RI) = ',EB6
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B6,RI) = ', EB6


    !==============================================================
    !=  B7: (ic|f12|jm)Fnm(ci|F12|nj)                             =
    !==============================================================
  
   !   We need CalphaGocc(NBA,nocc,nocc) but this is a subset of CalphaG(NBA,nocc,nbasis)
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'R' !Regular AO basis function on center 4
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g

   call mem_alloc(ABdecompC,nAux,nAux)
   ABdecompCreateC = .TRUE.
   call Build_CalphaMO2(MyFragment%MyLsitem,master,nbasis,nbasis,nAux,LUPRI,&
      & FORCEPRINT,wakeslaves,CoEOS,noccEOS,CoAOStot,noccAOStot,mynum,numnodes,CalphaCoccT,&
      & NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)
   call mem_dealloc(ABdecompC)

   !> Dgemm 
   nsize = nBA*noccEOS*noccAOStot
   call mem_alloc(CalphaD, nsize)
   m =  nBA*noccEOS                       ! D_jm = C_jn F_nm
   k =  noccAOStot
   n =  noccAOStot

   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaCoccT,m,Fnm,k,0.0E0_realk,CalphaD,m)
   !call dgemm('N','N',m,n,k,1.0E0_realk,CalphaCoccT,m,Fpp,k,0.0E0_realk,CalphaD,m)
   !May be optimized !!!
   call ContractTwo4CenterF12IntegralsRIB7_dec(nBA,noccEOS,noccAOStot,ncabsMO,CalphaGcabsMO, &
                                                & CalphaCoccT,CalphaD,EB7,dopair_occ)
   call mem_dealloc(CalphaCoccT)
   call mem_dealloc(CalphaD)
   mp2f12_energy = mp2f12_energy  + EB7
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B7,RI) = ',EB7
   WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B7,RI) = ', EB7

    !==============================================================
    !=  B8: (ic|f12|jm)Frm(ci|f12|rj)                             =
    !==============================================================

    !> Dgemm 
    nsize = nBA*noccEOS*noccAOStot
    call mem_alloc(CalphaD, nsize)
    m =  nBA*noccEOS                    ! D_jm = C_jr F_rm
    k =  ncabsAO
    n =  noccAOStot
    
    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,MyFragment%Frm,k,0.0E0_realk,CalphaD,m)
    call ContractTwo4CenterF12IntegralsRIB8_dec(nBA,noccEOS,ncabsMO,nocvAOS,noccAOStot,CalphaGcabsMO,CalphaG,CalphaD,EB8,dopair_occ)

    mp2f12_energy = mp2f12_energy  + EB8
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B8,RI) = ', EB8
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B8,RI) = ', EB8

    call mem_dealloc(CalphaGcabsAO)
    call mem_dealloc(CalphaD)

    !==============================================================
    !=  B9: (ip|f12|ja)Fcp(ci|f12|aj)                             =
    !==============================================================
    !> Dgemm 
    nsize = nBA*noccEOS*nocvAOS
    call mem_alloc(CalphaD, nsize)
    m =  nBA*noccEOS                    ! D_jp = C_jc F_cp
    k =  ncabsMO   
    n =  nocvAOS

    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsMO,m,MyFragment%Fcp,k,0.0E0_realk,CalphaD,m)
    
    call mem_dealloc(CalphaGcabsMO)
    
    call ContractTwo4CenterF12IntegralsRIB9_dec(nBA,noccEOS,nvirtAOS,nocvAOS,noccAOStot,CalphaG,&
         & CalphaD,EB9,dopair_occ)
!    call BuildKval(dopair_occ,noccEOS,Kval,noccpair,dopair)
!    call GeneralTwo4CenterDECF12RICoef1112(nBA,CalphaG,nocvAOS,noccEOS,CalphaD,noccEOS,nocvAOS,EB9,&
!         & noccfull,wakeslaves,use_bg_buf,numnodes,nAuxMPI,mynum,Kval,noccpair,DECF12RIB9,DECF12RIB9MPI)
!    call mem_dealloc(KVAL)

    mp2f12_energy = mp2f12_energy  + EB9
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B9,RI) = ', EB9
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B9,RI) = ', EB9

    ! ***********************************************************
    !    Free Memory
    ! ***********************************************************
    call mem_dealloc(dopair_occ)
    call mem_dealloc(CalphaCocc)

    call mem_dealloc(CalphaG)
    call mem_dealloc(ABdecompG)
    call mem_dealloc(CalphaD)
    
    call mem_dealloc(CoEOS)
    call mem_dealloc(CoAOS)
    call mem_dealloc(CoAOStot)
    call mem_dealloc(CvAOS)

    call mem_dealloc(Fac)
    call mem_dealloc(Fpp)
    call mem_dealloc(Fkj)
    call mem_dealloc(Fnm)

    call mem_dealloc(Cfull)
    call mem_dealloc(CMO_CABS)
    call mem_dealloc(CMO_RI)

    !> Need to be free to avoid memory leak for the type(matrix) CMO_RI in CABS.F90
    ! call free_cabs()

    call LSTIMER('DEC RIMP2:Step3',TS2,TE2,DECinfo%output,ForcePrint)
    call LSTIMER('DEC RIMP2F12',TS,TE,DECinfo%output,ForcePrint)

#ifdef VAR_MPI
   nbuf1 = 20
   CALL lsmpi_reduction(lsmpibufferRIMP2,nbuf1,infpar%master,infpar%lg_comm)
  
   if(master) then
      EV1=lsmpibufferRIMP2(2)
      EX1=lsmpibufferRIMP2(3)
      EB1=lsmpibufferRIMP2(1)
      EB2=lsmpibufferRIMP2(4)
      EB3=lsmpibufferRIMP2(5)
      print *,"After Reduction EV1:", EV1
      print *,"After Reduction EX1:", EX1
      print *,"After Reduction EB1:", EB1
      print *,"After Reduction EB2:", EB2
      print *,"After Reduction EB3:", EB3
   endif
#endif

    E_21 = 0.0E0_realk
    E_21 = EV1 + EV2 + EV3 + EV4 +EV5 

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, ' E21 V term                             '
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E21_CC_term:  ", E_21C
       write(*,'(1X,a,g25.16)') " E21_V_term1:  ", EV1
       write(*,'(1X,a,g25.16)') " E21_V_term2:  ", EV2
       write(*,'(1X,a,g25.16)') " E21_V_term3:  ", EV3
       write(*,'(1X,a,g25.16)') " E21_V_term4:  ", EV4
       write(*,'(1X,a,g25.16)') " E21_V_term5:  ", EV5
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E21_Vsum:     ", E_21

       write(DECinfo%output,'(1X,a,g25.16)') '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') ' E21 V term                             '
       write(DECinfo%output,'(1X,a,g25.16)') '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E21_CC_term:  ", E_21C
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term1:  ", EV1
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term2:  ", EV2
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term3:  ", EV3
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term4:  ", EV4
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term5:  ", EV5
       write(DECinfo%output,'(1X,a,g25.16)') '----------------------------------------'
       write(DECinfo%output,'(1X,a,f15.16)') " E21_Vsum:     ", E_21
    end if

    E_22 = 0.0E0_realk
    E_22 = EX1 + EX2 + EX3 + EX4 

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, ' E_22 X term                            '
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E22_X_term1: ", EX1
       write(*,'(1X,a,g25.16)') " E22_X_term2: ", EX2
       write(*,'(1X,a,g25.16)') " E22_X_term3: ", EX3
       write(*,'(1X,a,g25.16)') " E22_X_term4: ", EX4
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E22_Xsum:    ", E_22

       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E_22 X term                            '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E22_X_term1: ", EX1
       write(DECinfo%output,'(1X,a,g25.16)') " E22_X_term2: ", EX2
       write(DECinfo%output,'(1X,a,g25.16)') " E22_X_term3: ", EX3
       write(DECinfo%output,'(1X,a,g25.16)') " E22_X_term4: ", EX4
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E22_Xsum: ", E_22
    end if

    ! ***********************************************************
    !   Creating the B matrix 
    ! ***********************************************************

    E_23 = 0.0E0_realk
    E_23 = EB1 + EB2 + EB3 + EB4 + EB5 + EB6 + EB7 + EB8 + EB9

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, ' E_22 B term                            '
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E23_B_term1: ", EB1
       write(*,'(1X,a,g25.16)') " E23_B_term2: ", EB2 
       write(*,'(1X,a,g25.16)') " E23_B_term3: ", EB3   
       write(*,'(1X,a,g25.16)') " E23_B_term4: ", EB4   
       write(*,'(1X,a,g25.16)') " E23_B_term5: ", EB5   
       write(*,'(1X,a,g25.16)') " E23_B_term6: ", EB6  
       write(*,'(1X,a,g25.16)') " E23_B_term7: ", EB7   
       write(*,'(1X,a,g25.16)') " E23_B_term8: ", EB8   
       write(*,'(1X,a,g25.16)') " E23_B_term9: ", EB9  
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E23_B_sum:   ", E_23

       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E_22 B term                            '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term1: ", EB1
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term2: ", EB2   
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term3: ", EB3  
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term4: ", EB4   
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term5: ", EB5   
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term6: ", EB6   
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term7: ", EB7  
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term8: ", EB8  
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term9: ", EB9   
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_sum:   ", E_23
    end if

    E_F12 = 0.0E0_realk
    E_F12 = E_21 + E_22 + E_23

    !> MP2-energy from an MP2-calculation
    MP2_energy = Myfragment%energies(FRAGMODEL_OCCMP2)
    print *, "MP2_energy: ", MP2_energy

    if(DECinfo%F12debug) then
       print *,   '----------------------------------------------------------------'
       print *,   '                   DEC-MP2-F12 CALCULATION                      '
       print *,   '----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: MP2 CORRELATION ENERGY (For CC) =  ', MP2_energy
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E_21
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E_22
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E_23
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 CORRECTION TO ENERGY =         ', E_F12
       print *, '-------------------------------------------------------'
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: TOTAL CORRELATION ENERGY (For CC) =', MP2_energy+E_F12
    end if

    write(DECinfo%output,'(1X,a,f20.10)') '----------------------------------------------------------------'
    write(DECinfo%output,'(1X,a,f20.10)') '                  WANGY DEC-MP2-F12 CALCULATION                 '
    write(DECinfo%output,'(1X,a,f20.10)') '----------------------------------------------------------------'
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: MP2 CORRELATION ENERGY (For CC) =  ', MP2_energy
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E_21
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E_22
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E_23
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 CORRECTION TO ENERGY =         ', E_F12
    write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: MP2-F12 CORRELATION ENERGY (CC) =  ', MP2_energy+E_F12

    !> Setting the MP2-F12 correction
    Myfragment%energies(FRAGMODEL_RIMP2f12) = E_F12


  end subroutine get_rif12_fragment_energy

  subroutine BuildKval(dopair_occ,noccEOS,Kval,noccpair,dopair)
    implicit none
    integer,intent(in) :: noccEOS
    logical,intent(in) :: dopair_occ(noccEOS,noccEOS),dopair
    integer,pointer    :: KVAL(:,:)
    integer,intent(inout) :: noccpair
    !local variables
    integer :: nK,J,I
    IF(dopair)THEN
       noccpair=COUNT(dopair_occ)/2
       print*,'PAIR: noccpair=COUNT(dopair_occ)/2=',COUNT(dopair_occ)/2
    ELSE
       nK=NINT(SQRT(real(COUNT(dopair_occ))))
       noccpair=(nK*(nK+1))/2
       print*,'ATOMIC: noccpair=',noccpair
    ENDIF
    call mem_alloc(KVAL,3,noccpair)
    nK=0
    DO j=1,noccEOS 
       IF(dopair_occ(J,J)) THEN
          nK=nK+1
          KVAL(1,nK)=J
          KVAL(2,nK)=J
          KVAL(3,nK)=1.0E0_realk
       ENDIF
       DO i=j+1,noccEOS 
          IF(dopair_occ(I,J)) THEN
             nK=nK+1
             KVAL(1,nK)=I
             KVAL(2,nK)=J
             KVAL(3,nK)=2.0E0_realk
          ENDIF
       ENDDO
    ENDDO
    IF(nK.NE.noccpair)call lsquit('BuildKval: nK.NE.noccpair',-1)
  end subroutine BuildKval
  
end module rif12_integrals_module

subroutine get_rif12_fragment_energy_slave()
   use rif12_integrals_module
#ifdef VAR_MPI   
   use infpar_module
   use lsmpi_type
#endif
   use precision
   use dec_typedef_module
   !use memory_handling
   use lsparameters
   ! *********************************************
   !   DEC DEPENDENCIES (within deccc directory) 
   ! *********************************************
#ifdef VAR_MPI
   use decmpi_module !, only: mpi_communicate_mp2_int_and_amp
#endif
   implicit none

  !> Atomic fragment to be determined (Single or Pair fragment)
   type(decfrag) :: MyFragment
   !> t2EOS amplitudes stored in the order T(a,i,b,j)
   real(realk), pointer :: Taibj(:,:,:,:) 
   !> t1EOS amplitudes stored in the order T(a,i)
   real(realk), pointer :: Tai(:,:) 
   !> Case MODEL
   integer :: fragcase
   !> Fragment 1 in the pair fragment
   type(decfrag) :: Fragment1
   !> Fragment 2 in the pair fragment
   type(decfrag) :: Fragment2
   !> Logical variable to check if this is a pair fragment
   logical :: dopair

   !print *, "I am slave nr: ", infpar%lg_mynum 
#ifdef VAR_MPI
   call mpi_communicate_MyFragment_f12(MyFragment,Taibj,fragcase,dopair)
   if(dopair) then
      call get_rif12_fragment_energy(MyFragment, Taibj, Tai, fragcase, Fragment1, Fragment2) 
   else
      call get_rif12_fragment_energy(MyFragment, Taibj, Tai, fragcase)
   endif
#endif

end subroutine get_rif12_fragment_energy_slave
