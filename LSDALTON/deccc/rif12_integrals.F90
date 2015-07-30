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

  !  use orbital_operations

  use ccintegrals!, only: get_full_AO_integrals,get_AO_hJ,get_AO_K,get_AO_Fock

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

  !#ifdef MOD_UNRELEASED
  ! Yangs F12 routines
  use f12_routines_module!, only: MO_transform_AOMatrix, matrix_print, norm4D, norm2D, get_mp2f12_MO
  use rimp2_module
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
#endif

contains
#ifdef MOD_UNRELEASED 

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
    integer, intent(in) :: fragcase
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
    integer :: nocc, nvirt, noccfull
    !> number of occupied MO orbitals in AOS 
    integer :: noccAOS
    !> number of occupied MO orbitals in AOS tot 
    integer :: noccAOStot
    !> number of occupied MO + virtual MO orbitals in AOS 
    integer :: nocv  
    !> number of CABS AO orbitals
    integer :: ncabsAO
    !> number of CABS MO orbitals
    integer :: ncabsMO

    !> MO coefficient matrix for the CABS MOs
    real(realk), pointer :: CMO_CABS(:,:)
    !> MO coefficient matrix for the RI MOs
    real(realk), pointer :: CMO_RI(:,:)
    !> MO coefficient matrix for the OCC + VIRT
    real(realk), pointer :: Cfull(:,:)

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
    real(realk) :: eps
    real(realk),pointer :: Fac(:,:)
    real(realk),pointer :: Fpp(:,:)
    !========================================================
    !  RI variables
    !========================================================
    real(realk),pointer :: CalphaR(:),CalphaG(:),CalphaF(:),CalphaD(:),CalphaCvirt(:)
    real(realk),pointer :: CalphaRcabsMO(:),CalphaGcabsAO(:),CalphaX(:),CalphaCcabs(:)
    real(realk),pointer :: CalphaGcabsMO(:),CalphaXcabsAO(:)
    real(realk),pointer :: ABdecompR(:,:),ABdecompG(:,:),ABdecompC(:,:)
    real(realk),pointer :: ABdecompF(:,:),Umat(:,:),Rtilde(:,:),ABdecompX(:,:)
    logical :: ABdecompCreateR,ABdecompCreateG,ABdecompCreateF,ABdecompCreateC
    logical :: FORCEPRINT,use_bg_buf,LS,ABdecompCreateX

    character :: intspec(5)
    real(realk) :: TS,TE,TS2,TE2
    real(realk) :: ExchangeF12V1,CoulombF12V1
    real(realk) :: CoulombF12X1,ExchangeF12X1 
    real(realk) :: EV1,EV2,EV3,EV4,EV5,EX1,EX2,EX3,EX4
    real(realk) :: EB1,EB2,EB3,EB4,EB5,EB6,EB7,EB8,EB9
    integer :: lupri

    !> Lsitem structure
    type(lsitem) :: mylsitem

    !> Timings
    real(realk) :: tcpu,twall
    logical :: Collaborate,DoBasis
    integer :: n1,n2,n3,n4,Tain1,Tain2

    !> bat
    type(mp2_batch_construction) :: bat

    !> Number of Auxilliary Basis functions
    integer :: NBA
    !> MPI related variables
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

    mylsitem = MyFragment%MyLsitem

    lupri = DECinfo%output
#ifdef VAR_TIME
    FORCEPRINT = .TRUE.
#else
    FORCEPRINT = .FALSE.
#endif    
    call LSTIMER('START ',TS,TE,DECinfo%output,ForcePrint)

#ifdef VAR_MPI
    master= (infpar%mynum == infpar%master)
    mynum = infpar%mynum
    numnodes = infpar%nodtot
    wakeslaves = infpar%nodtot.GT.1
    if(.NOT.master)lupri = 6
#else
    ! If MPI is not used, consider the single node to be "master"
    master=.true.
    mynum = 0
    numnodes = 1
    wakeslaves = .false.
#endif
    call LSTIMER('START ',TS2,TE2,DECinfo%output,ForcePrint)

    !========================================================
    !  Init stuff
    !========================================================
    natoms = MyFragment%natoms
    nbasis = MyFragment%nbasis
    nvirt  = MyFragment%nvirtAOS
    nocc   = MyFragment%noccEOS
    noccAOStot = MyFragment%nocctot
    noccAOS = MyFragment%noccAOS
    ncabsAO = size(MyFragment%Ccabs,1)
    ncabsMO = size(MyFragment%Ccabs,2)

    nocv = noccAOS + nvirt

    nbasis2 = 0
    natomsaux = 0
    naux = 0
    E_21C = 0.0E0_realk

    ! Offset: Used for frozen core
    if(DECinfo%frozencore) then
       offset = MyFragment%ncore
    else
       offset=0
    end if
    noccfull = nocc

    IF(DECinfo%frozencore)call lsquit('DEC RI-/MP2-F12 frozen core not implemented',-1)

    !========================================================
    !  Special treatment of Aux
    !========================================================
    call determine_maxBatchOrbitalsize(DECinfo%output,MyFragment%mylsitem%SETTING,MinAuxBatch,'D')

    IF(DECinfo%AuxAtomicExtent)THEN
       call getMolecularDimensions(MyFragment%mylsitem%INPUT%AUXMOLECULE,nAtomsAux,nBasis2,nAux)
    ELSE
       print *, "natomsaux nbasis2 naux",nAtomsAux,nBasis2,nAux
       call getMolecularDimensions(MyFragment%mylsitem%SETTING%MOLECULE(1)%p,nAtomsAux,nBasis2,nAux)
       print *, "natomsaux nbasis2 naux",nAtomsAux,nBasis2,nAux
       if(natoms.NE.natomsAux)call lsquit('Error in DEC RIMP2F12 natoms dim mismatch',-1)
    ENDIF
    print *, "END Special treatment of Aux"

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

    !========================================================
    !  Printing Input variables 
    !========================================================
    if(DECinfo%F12debug) then
       print *, "-------------------------------------------------"
       print *, "     F12-integrals.F90                           "
       print *, "-------------------------------------------------"
       print *, "nbasis:       ", nbasis
       print *, "noccEOS:      ", nocc
       print *, "nvirtAOS:     ", nvirt
       print *, "-------------------------------------------------"
       print *, "noccAOS       ", noccAOS 
       print *, "noccAOStot    ", noccAOStot
       print *, "nocc+nvirt    ", nocc+nvirt   
       print *, "noccAOS+nvirt ", noccAOS+nvirt
       print *, "ncabsAO       ", ncabsAO
       print *, "ncabsMO       ", ncabsMO
    end if

    !========================================================
    !  Creating Coeff-matrices 
    !========================================================
    call mem_alloc(CMO_Cabs, ncabsAO, ncabsMO)
    do i=1, ncabsMO
       CMO_Cabs(:,i) = MyFragment%Ccabs(:,i)
    end do

    call mem_alloc(CMO_RI, ncabsAO, ncabsAO)
    do i=1, ncabsAO
       CMO_RI(:,i) = MyFragment%Cri(:,i)
    end do

    call mem_alloc(Cfull,nbasis,nocv)
    do J=1,noccAOS
       do I=1,nbasis
          Cfull(I,J) = MyFragment%Co(I,J)
       enddo
    enddo
    do P=1,nvirt
       do I=1,nbasis
          Cfull(I,noccAOS+P) = MyFragment%Cv(I,P)
       enddo
    enddo

    call mem_alloc(Fac, nvirt, ncabsMO)
    do a=1, nvirt
       do c=1, ncabsMO
          Fac(a,c) = Myfragment%Fcp(c,a+noccAOStot)
       enddo
    enddo

    !Slow implementation can be done with 2dgemms in B6
    call mem_alloc(Fpp, noccAOS+nvirt, noccAOS+nvirt)
    Fpp = 0.0E0_realk
    do p=1, noccAOS
       do q=1, noccAOS
          Fpp(p,q) = MyFragment%ppfock(p,q)
       enddo
    enddo
    do p=noccAOS+1, noccAOS+nvirt
       do q=noccAOS+1, noccAOS+nvirt
          Fpp(p,q) = MyFragment%qqfock(p-noccAOS,q-noccAOS)
       enddo
    enddo

    !=================================================================
    != Step 0: Creating of dopair_occ                                =
    !=================================================================
    print *, "Step 0: Creating of dopair_occ"
    print *, "dopair: ", dopair
    call mem_alloc(dopair_occ,nocc,nocc)
    dopair_occ = .FALSE.
    if(dopair) then 
       call which_pairs_occ(Fragment1,Fragment2,MyFragment,dopair_occ)
    else
       dopair_occ = .TRUE.
    endif
    !=================================================================
    != Step 1:  Fijkl,Xijkl,Dijkl                                    =
    !=          corresponding to V1,X1,B1                            =
    !=================================================================
    print *, "Step 1:  Fijkl,Xijkl,Dijkl "
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
    call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,MyFragment%Co,nocc,&
         & mynum,numnodes,CalphaF,NBA,ABdecompF,ABdecompCreateF,intspec,use_bg_buf)
    ABdecompCreateF = .FALSE.
    !perform this suborutine on the GPU (async)  - you do not need to wait for the results
    call ContractOne4CenterF12IntegralsRI(NBA,nocc,CalphaF,CoulombF12V1,ExchangeF12V1,dopair_occ)
    
    !Calculate the Fitting Coefficients (alpha|g^2|ij) 
    intspec(4) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
    intspec(5) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
    call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,MyFragment%Co,nocc,&
         & mynum,numnodes,CalphaX,NBA,ABdecompX,ABdecompCreateX,intspec,use_bg_buf)
    ABdecompCreateX = .FALSE.

    !perform this suborutine on the GPU (async)  - you do not need to wait for the results
    call ContractOne4CenterF12IntegralsRI2(NBA,nocc,CalphaX,MyFragment%ppfock,CoulombF12X1,ExchangeF12X1,dopair_occ)

    !Calculate the Fitting Coefficients (alpha|[[T,g],g]|ij) 
    intspec(4) = 'D' !The double commutator [[T,g],g] 
    intspec(5) = 'C' !The metric operator = Regular Coulomb operator 1/r12
    !Build the R coefficient of Eq. 90 of J Comput Chem 32: 2492–2513, 2011
    call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,MyFragment%Co,nocc,&
         & mynum,numnodes,CalphaD,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
    ABdecompCreateR = .FALSE.
    !We need CalphaR(NBA,nocc,nocc) but this is a subset of the CalphaR(NBA,nocc,nbasis)
    !so we calculate the full CalphaR(NBA,nocc,nbasis)
    intspec(4) = 'C' !Regular Coulomb operator 1/r12
    intspec(5) = 'C' !The metric operator = Regular Coulomb operator 1/r12
    !Build the G coefficient of Eq. 90 of J Comput Chem 32: 2492–2513, 2011
    call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,Cfull,nbasis,&
         & mynum,numnodes,CalphaR,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
    !Build the U matrix in Eq. 88 of J Comput Chem 32: 2492–2513, 2011
    call mem_alloc(Umat,nAux,nAux)
    !perform this suborutine on the GPU (Async)
    call Build_RobustERImatU(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,MyFragment%Co,nocc,&
         & mynum,numnodes,ABdecompR,'D',Umat)
    !Build the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 2011
    M = NBA          !rows of Output Matrix
    N = nocc*nocc    !columns of Output Matrix
    K = NBA          !summation dimension
    !perform this suborutine on the GPU (Async)
    !note CalphaR is actual of dimensions (NBA,nocc,nbasis) but here we only access
    !the first part (NBA,nocc,nocc) 
    call dgemm('N','N',M,N,K,1.0E0_realk,Umat,M,CalphaR,K,-0.5E0_realk,CalphaD,M)
    call mem_dealloc(Umat)
    !CalphaD is now the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 20
    !perform this suborutine on the GPU (Async)
    call ContractOne4CenterF12IntegralsRobustRI(nAux,nocc,nbasis,CalphaD,CalphaR,EB1,dopair_occ)

    !The minus is due to the Valeev factor
    EV1 = -1.0E0_realk*((5.0E0_realk*0.25E0_realk)*CoulombF12V1-ExchangeF12V1*0.25E0_realk)
    mp2f12_energy = mp2f12_energy + EV1
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V1,RI) = ',EV1
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V1,RI) = ',EV1
    !minus is due to the overall minus from equation (41) and (42) due to
    !contribution from the \bar{B}_{ij}^{ij}
    EX1 = -1.0E0_realk*(0.21875E0_realk*CoulombF12X1 + 0.03125E0_realk*ExchangeF12X1)
    mp2f12_energy = mp2f12_energy + EX1
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X1,RI) = ',EX1
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X1,RI) = ',EX1
    mp2f12_energy = mp2f12_energy + EB1
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B1,RI) = ', EB1
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B1,RI) = ', EB1

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
    call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,CMO_RI,ncabsAO,&
         & mynum,numnodes,CalphaXcabsAO,NBA,ABdecompX,ABdecompCreateX,intspec,use_bg_buf)

    call ContractOne4CenterF12IntegralsRIB23(nBA,nocc,ncabsAO,CalphaXcabsAO,CalphaX,&
         & MyFragment%hJir,1.0E0_realk,EB2,EB3,dopair_occ)
    !1.0E0_realk because that term has an overall pluss in Eqs. 25-26
    mp2f12_energy = mp2f12_energy + EB2
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B2,RI) = ',EB2
    mp2f12_energy = mp2f12_energy + EB3
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B3,RI) = ',EB3
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B2,RI) = ',EB2
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B3,RI) = ',EB3

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
    != Dim(nocc,nbasis,nocc,nbasis)                           =
    !=                                                        =
    !==========================================================
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'R' !Regular AO basis function on center 4
    intspec(4) = 'G' !The Gaussian geminal operator g
    intspec(5) = 'G' !The Gaussian geminal operator g

    call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,Cfull,nocv,&
         & mynum,numnodes,CalphaG,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
    ABdecompCreateG = .FALSE.

    !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.

    call ContractTwo4CenterF12IntegralsRI(nBA,nocc,nocv,CalphaR,CalphaG,EV2,dopair_occ)
    mp2f12_energy = mp2f12_energy + EV2
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V2,RI) = ',EV2       
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V2,RI) = ',EV2

    !==========================================================
    !=                                                        =
    != X2: Gipjq*Gipjq                                        =
    != The Gaussian geminal operator Int multiplied with      =
    != The Gaussian geminal operator g                        =
    != Dim(nocc,nbasis,nocc,nbasis)                           =
    !=                                                        =
    !==========================================================
    !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
    call ContractTwo4CenterF12IntegralsRIX(nBA,nocc,nocv,CalphaG,MyFragment%ppfock,EX2,dopair_occ)
    mp2f12_energy = mp2f12_energy + EX2
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X2,RI) = ',EX2       
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(X2,RI) = ',EX2

    !==========================================================
    !=                                                        =
    != V3: Rimjc*Gimjc                                        =
    != V4: Rjmic*Gjmic                                        =
    != The Coulomb Operator Int multiplied with               =
    != The Gaussian geminal operator g                        =
    != Dim: (nocc,noccfull,nocc,ncabsMO)  need 4 Calphas      =
    !=                                                        =
    !==========================================================

    !   We need CalphaRocc(NBA,nocc,nocc) but this is a subset of CalphaR(NBA,nocc,nbasis)
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'C' !CABS AO basis function on center 4
    intspec(4) = 'C' !The Coulomb Operator
    intspec(5) = 'C' !The Coulomb Operator
    call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,CMO_CABS,ncabsMO,&
         & mynum,numnodes,CalphaRcabsMO,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)

    !   We need CalphaGocc(NBA,nocc,nocc) but this is a subset of CalphaG(NBA,nocc,nbasis)
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'C' !CABS AO basis function on center 4
    intspec(4) = 'G' !The Gaussian geminal operator g
    intspec(5) = 'G' !The Gaussian geminal operator g
    call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,CMO_CABS,ncabsMO,&
         & mynum,numnodes,CalphaGcabsMO,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)

    !Do on GPU (Async)
    call ContractTwo4CenterF12IntegralsRI2V3V4(NBA,nocc,noccfull,ncabsMO,nbasis,&
         & CalphaRcabsMO,CalphaGcabsMO,CalphaR,CalphaG,EV3,EV4,dopair_occ)

    mp2f12_energy = mp2f12_energy + EV3 + EV4
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V3,RI) = ',EV3       
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V4,RI) = ',EV4       
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V3,RI) = ',EV3       
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V4,RI) = ',EV4       

    call mem_dealloc(ABdecompR)
    call mem_dealloc(CalphaR)
    call mem_dealloc(CalphaRcabsMO)

    !==========================================================
    != V5: Caibj = (Gcibj*Fac + Gcjai*Fcb)*Taibj              =
    !========================================================== 

    ! **********************
    !  Get all AO integrals 
    ! **********************
    call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')
    ! Transform AO integrals to MO integrals (A I | B J)
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyFragment%Co, MyFragment%Cv,'aiai',gAO,gMO)
    call mem_dealloc(gao)

    call mem_alloc(ABdecompC,nAux,nAux)
    ABdecompCreateC = .TRUE.
    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 4 
    intspec(3) = 'C' !Cabs AO basis function on center 3
    intspec(4) = 'G' !The Gaussian geminal operator g
    intspec(5) = 'G' !The Gaussian geminal operator g

    call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,CMO_CABS,ncabsMO,&
         & mynum,numnodes,CalphaCcabs,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)

    m = NBA*nocc       
    k = ncabsMO         ! C_mn = A_mk B_kn
    n = nvirt  

    !C(alpha*i,ncabsMO)*F(cabsMO,nvirt)
    nsize = nBA*nocc*ncabsMO
    call mem_alloc(CalphaD, nsize)
    call dgemm('N','T',m,n,k,1.0E0_realk,CalphaCcabs,m,Fac,n,0.0E0_realk,CalphaD,m)

    !m = nvirt      
    !k = ncabsMO         ! C_mn = A_mk B_kn
    !n = nocc*NBA  
    !call dgemm('N','T',m,n,k,1.0E0_realk,Fac%elms,m,CalphaCcabs,k,0.0E0_realk,CalphaD,n)

    intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
    intspec(2) = 'R' !Regular AO basis function on center 3
    intspec(3) = 'R' !Regular AO basis function on center 4
    intspec(4) = 'G' !The Gaussian geminal operator g
    intspec(5) = 'G' !The Gaussian geminal operator g
    call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,MyFragment%Cv,nvirt,&
         & mynum,numnodes,CalphaCvirt,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)

    call ContractTwo4CenterF12IntegralsRIC(nBA,nocc,nvirt,CalphaCvirt,CalphaD,Taibj,EV5,dopair_occ)
    mp2f12_energy = mp2f12_energy + EV5
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V5,RI) = ', EV5
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(V5,RI) = ', EV5

    ABdecompCreateG = .FALSE.
    call mem_dealloc(CalphaD)
    call mem_dealloc(ABdecompC)
    call mem_dealloc(CalphaCcabs)
    call mem_dealloc(CalphaCvirt)

    !Additional 
    call mem_dealloc(gMO)

    !==========================================================
    !=                                                        =
    != X3:         Step 3  Gimjc*Gimjc                        =
    != X4:         Step 4  Gjmic*Gjmic                        =
    != The Coulomb Operator Int multiplied with               =
    != The Gaussian geminal operator g                        =
    != Dim: (nocc,noccfull,nocc,ncabsMO)  need 4 Calphas      =
    !=                                                        =
    !==========================================================
    !   We need CalphaGocc(NBA,nocc,nocc) but this is a subset of CalphaG(NBA,nocc,nbasis)

    !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
    call ContractTwo4CenterF12IntegralsRI2X(NBA,nocc,noccAOS,ncabsMO, &
         & CalphaGcabsMO,CalphaG,MyFragment%ppfock,EX3,EX4,dopair_occ)

    mp2f12_energy = mp2f12_energy + EX3 + EX4
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
    call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
         & FORCEPRINT,wakeslaves,MyFragment%Co,nocc,CMO_RI,ncabsAO,&
         & mynum,numnodes,CalphaGcabsAO,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)

    nsize = nBA*nocc*ncabsAO
    call mem_alloc(CalphaD, nsize)
    m =  nBA*nocc                    ! C_mn = A_mk B_kn
    k =  ncabsAO
    n =  ncabsAO

    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,MyFragment%Krs,k,0.0E0_realk,CalphaD,m)
    call ContractTwo4CenterF12IntegralsRIB4(nBA,nocc,ncabsAO,CalphaGcabsAO,CalphaD,EB4,dopair_occ)

    mp2f12_energy = mp2f12_energy + EB4
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B4,RI) = ',EB4
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B4,RI) = ', EB4

    !==============================================================
    !=  B5: (ir|f12|jm)Fsr(si|f12|mj)        (r,s=CabsAO)         =
    !==============================================================   
    !We need CalphaG(NBA,nocc,noccfull) but this is a subset of 
    !CalphaG(NBA,nocc,nbasis) which we already have
    !> Dgemm 
    nsize = nBA*nocc*ncabsAO
    IF(size(CalphaD).NE.nsize)call lsquit('dim mismatch CalphaD',-1)
    m =  nBA*nocc                    ! D_jq = C_jp F_qp
    k =  ncabsAO
    n =  ncabsAO

    !Do on GPU (Async)
    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,Myfragment%Frs,k,0.0E0_realk,CalphaD,m)
    !Do on GPU (Async)
    call ContractTwo4CenterF12IntegralsRIB5(nBA,nocc,ncabsAO,noccAOS,CalphaGcabsAO,CalphaG,CalphaD,EB5,dopair_occ)

    mp2f12_energy = mp2f12_energy + EB5
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B5,RI) = ',EB5
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B5,RI) = ', EB5

    call mem_dealloc(CalphaD)

    !==============================================================
    !=  B6: (ip|f12|ja)Fqp(qi|f12|aj)                             =
    !==============================================================
    !> Dgemm 
    nsize = nBA*nocc*(nvirt+noccAOS)
    call mem_alloc(CalphaD, nsize)
    m =  nBA*nocc              ! D_jq = C_jp F_pq
    k =  nvirt+noccAOS
    n =  nvirt+noccAOS

    !Do on GPU (Async)
    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaG,m,Fpp,k,0.0E0_realk,CalphaD,m)   
    !Do on GPU (Async)
    call ContractTwo4CenterF12IntegralsRIB6(nBA,nocc,nvirt,nocv,CalphaG,CalphaD,EB6,dopair_occ)

    mp2f12_energy = mp2f12_energy  + EB6
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B6,RI) = ',EB6
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B6,RI) = ', EB6

    call mem_dealloc(CalphaD)

    !==============================================================
    !=  B7: (ic|f12|jm)Fnm(ci|F12|nj)                             =
    !==============================================================
    !We need CalphaG(NBA,nocc,noccfull) but this is a subset of CalphaG(NBA,nocc,nbasis) 
    !that we already have

    !> Dgemm 
    nsize = nBA*nocc*noccAOS
    call mem_alloc(CalphaD, nsize)
    m =  nBA*nocc               ! D_jn = C_jn F_nm
    k =  noccAOS
    n =  noccAOS

    !NB! Changed T to N, dont think it will matter but...
    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaG,m,MyFragment%ppfock,k,0.0E0_realk,CalphaD,m)   
    ! Do not compute T.K. needs to explain this...
    !call ContractOccCalpha(NBA,nocc,noccAOS,nbasis,CalphaG,MyFragment%ppfock,CalphaD,dopair_occ)
    call ContractTwo4CenterF12IntegralsRIB7(nBA,nocc,ncabsMO,noccAOS,CalphaGcabsMO,CalphaG,CalphaD,EB7,dopair_occ)

    mp2f12_energy = mp2f12_energy  + EB7
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B7,RI) = ',EB7
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B7,RI) = ', EB7

    call mem_dealloc(CalphaD)

    !==============================================================
    !=  B8: (ic|f12|jm)Frm(ci|f12|rj)                             =
    !==============================================================

    !> Dgemm 
    nsize = nBA*nocc*nocc
    call mem_alloc(CalphaD, nsize)
    m =  nBA*nocc                    ! D_jm = C_jp F_pm
    n =  nocc   
    k =  ncabsAO
    !NB! Changed T to N, dont think it will matter but...
    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,MyFragment%Frm,k,0.0E0_realk,CalphaD,m)

    !we need CalphaG(NBA,nocc,noccfull) but this is a subset of CalphaG(NBA,nocc,nbasis)
    call ContractTwo4CenterF12IntegralsRIB8(nBA,nocc,ncabsMO,nocv,CalphaGcabsMO,CalphaG,CalphaD,EB8,dopair_occ)

    mp2f12_energy = mp2f12_energy  + EB8
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B8,RI) = ', EB8
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B8,RI) = ', EB8

    call mem_dealloc(CalphaGcabsAO)
    call mem_dealloc(CalphaD)

    !==============================================================
    !=  B9: (ip|f12|ja)Fcp(ci|f12|aj)                             =
    !==============================================================
    !> Dgemm 
    nsize = nBA*nocc*nocv
    call mem_alloc(CalphaD, nsize)
    m =  nBA*nocc                    ! D_jp = C_jc F_cp
    n =  nocv
    k =  ncabsMO   
    call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsMO,m,MyFragment%Fcp,k,0.0E0_realk,CalphaD,m)
    call mem_dealloc(CalphaGcabsMO)

    call ContractTwo4CenterF12IntegralsRIB9(nBA,nocc,nvirt,nocv,CalphaG,CalphaD,EB9,dopair_occ)

    mp2f12_energy = mp2f12_energy  + EB9
    WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B9,RI) = ', EB9
    WRITE(*,'(A50,F20.13)')'DEC RIMP2F12 Energy contribution: E(B9,RI) = ', EB9

    ! ***********************************************************
    !    Free Memory
    ! ***********************************************************
    call mem_dealloc(dopair_occ)

    call mem_dealloc(Fac)
    call mem_dealloc(Fpp)

    call mem_dealloc(CalphaG)
    call mem_dealloc(ABdecompG)
    call mem_dealloc(CalphaD)

    call mem_dealloc(Cfull)
    call mem_dealloc(CMO_CABS)
    call mem_dealloc(CMO_RI)

    !> Need to be free to avoid memory leak for the type(matrix) CMO_RI in CABS.F90
    ! call free_cabs()

    call LSTIMER('DEC RIMP2:Step3',TS2,TE2,DECinfo%output,ForcePrint)
    call LSTIMER('DEC RIMP2F12',TS,TE,DECinfo%output,ForcePrint)

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

#else

  subroutine wangy_dummy_module()
    implicit none
  end subroutine wangy_dummy_module

#endif

end module rif12_integrals_module

