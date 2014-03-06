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

  ! Yangs F12 routines
  use f12_routines_module!, only: MO_transform_AOMatrix, matrix_print

  ! Patricks mat_transpose routine 
  use reorder_frontend_module!, only: mat_transpose(rows,column,pref1,A,pref2,AT)

  ! Thomas free_cabs() for aa free MO_CABS_save_created, CMO_RI_save_created
  use CABS_operations

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

  public :: get_f12_fragment_energy, get_f12_pair_fragment_energy, matrix_print_4d, matrix_print_2d, get_mp2f12_sf_E21

  private

  !> Coefficient Type
  TYPE ctype
     real(realk), pointer :: cmat(:,:)
     integer :: n1
     integer :: n2
  END TYPE ctype

contains
  !> Brief: Gives the single fragment energy for MP2F12
  !> Author: Yang M. Wang
  !> Date: April 2013
  subroutine get_f12_fragment_energy(MyFragment,Fragment1,Fragment2,natoms)
    implicit none

    !> Atomic fragment to be determined (Single or Pair fragment)
    type(decfrag), intent(inout) :: MyFragment
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(inout), optional :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(inout), optional :: Fragment2
    !> Number of atoms for full molecule
    integer, intent(in), optional :: natoms
    !> Logical variable to check if this is a pair fragment
    logical :: dopair

    ! ***********************************************************
    !   Allocating for Coefficient matrix
    ! ***********************************************************
    !> MO coefficient matrix for the occupied EOS
    real(realk), pointer :: CoccEOS(:,:)
    !> MO coefficient matrix for the occupied AOS
    real(realk), pointer :: CoccAOS(:,:)
    !> MO coefficient matrix for the virtual AOS
    real(realk), pointer :: CvirtAOS(:,:)
    !> MO coefficient matrix for the occupied + virtual AOS
    real(realk), pointer :: CocvAOS(:,:)
    !> MO coefficient matrix for the CABS MOs
    real(realk), pointer :: Ccabs(:,:)
    !> MO coefficient matrix for the RI MOs
    real(realk), pointer :: Cri(:,:)

    ! ***********************************************************
    !   Allocating for V matrix
    ! ***********************************************************
    !> F12 integrals for the V1_ijkl <ij|f12*r^-1|kl>
    real(realk), pointer :: V1ijkl(:,:,:,:) 
    !> F12 integrals for the V2_ijkl sum_pq <ij|r^-1|pq> * <pq|f12|kl>
    real(realk), pointer :: V2ijkl(:,:,:,:)   
    real(realk), pointer :: Gijpq(:,:,:,:) ! <ij|r^-1|pq>
    real(realk), pointer :: Rijpq(:,:,:,:) ! <ij|f12|pq>
    !> F12 integrals for the V3_term sum_mc <ij|r^-1|mc> * <mc|f12|kl>
    real(realk), pointer :: V3ijkl(:,:,:,:)   
    real(realk), pointer :: Gijmc(:,:,:,:) ! <ij|r^-1|ma'>
    real(realk), pointer :: Rijmc(:,:,:,:) ! <ij|f12|ma'>
    !> F12 integrals for the V4_term sum_mc <ji|r^-1|cm> * <cm|f12|lk>
    real(realk), pointer :: V4ijkl(:,:,:,:) ! Not necessary 

    ! ***********************************************************
    !   Allocating for C matrix
    ! ***********************************************************
    !> Fock Fij 
    real(realk), pointer :: Fij(:,:)
    !> Fock Fmn 
    real(realk), pointer :: Fmn(:,:)
    !> Fock Fab
    real(realk), pointer :: Fab(:,:)
    !> Fock Fpq
    real(realk), pointer :: Fpq(:,:)
    !> Cijab
    real(realk), pointer :: Cijab(:,:,:,:)

    ! ***********************************************************
    !   Allocating for X matrix
    ! ***********************************************************
    !> F12 integral X1ijkl
    real(realk), pointer :: X1ijkl(:,:,:,:)
    !> F12 integrals for the X2_term sum_pq <ij|g|pq> * <pq|g|kl>
    real(realk), pointer :: X2ijkl(:,:,:,:) 
    !> F12 integrals for the X3_term sum_pq <ij|g|mc> * <mc|g|kl>
    real(realk), pointer :: X3ijkl(:,:,:,:) 
    !> F12 integrals for the X4_term sum_pq <ij|g|cm> * <cm|g|kl>
    real(realk), pointer :: X4ijkl(:,:,:,:)    

    ! ***********************************************************
    !   Allocating for B matrix
    ! ***********************************************************
    !> F12 integrals for the B1_term <ij|[[T,f12],f12]|kl> 
    real(realk), pointer :: B1ijkl(:,:,:,:)   
    !> F12 integrals for the B2_term <ij|f12^2|rk>  r = RI MO   
    real(realk), pointer :: B2ijkl(:,:,:,:)     
    real(realk), pointer :: R2rlij(:,:,:,:)
    !> F12 integrals for the B3_term <ij|f12^2|kr>  r = RI MO         
    real(realk), pointer :: B3ijkl(:,:,:,:)
    real(realk), pointer :: R2ijkr(:,:,:,:)
    !> F12 integrals for the B4_term
    real(realk), pointer :: B4ijkl(:,:,:,:)
    real(realk), pointer :: Rijrs(:,:,:,:)
    !> F12 integrals for the B5_term
    real(realk), pointer :: B5ijkl(:,:,:,:)
    real(realk), pointer :: Rijrm(:,:,:,:)
    !> F12 integrals for the B6_term
    real(realk), pointer :: B6ijkl(:,:,:,:)
    real(realk), pointer :: Rijpa(:,:,:,:)
    !> F12 integrals for the B7_term
    real(realk), pointer :: B7ijkl(:,:,:,:)
    real(realk), pointer :: Rijcm(:,:,:,:)
    !> F12 integrals for the B8_term
    real(realk), pointer :: B8ijkl(:,:,:,:)
    !  real(realk), pointer :: Rijcm(:,:,:,:)
    real(realk), pointer :: Rijcr(:,:,:,:) 
    !> F12 integrals for the B9_term   
    real(realk), pointer :: B9ijkl(:,:,:,:)                                                                     
    !  real(realk), pointer :: Rijap(:,:,:,:)
    real(realk), pointer :: Rijca(:,:,:,:)

    ! ***********************************************************
    !   Allocating integer space sizes
    ! ***********************************************************
    !> number of AO orbitals
    integer :: nbasis
    !> number of occupied MO orbitals in EOS 
    integer :: noccEOS, nunoccEOS, noccfull
    !> number of occupied + virtual MO orbitals in EOS 
    integer :: nocvAOS  
    !> number of virtual MO orbitals in AOS 
    integer :: nvirtAOS
    !> number of occupied MO orbitals in AOS 
    integer :: noccAOS

    !> number of CABS AO orbitals
    integer :: ncabsAO
    !> number of CABS MO orbitals
    integer :: ncabsMO

    integer :: ix, iy, i, j, m, n, k, l, p, q, c, r, s, t, a, b

    real(realk) :: V1energy, V2energy, V3energy, V4energy
    real(realk) :: X1energy, X2energy, X3energy, X4energy
    real(realk) :: B1energy, B2energy, B3energy, B4energy
    real(realk) :: B5energy, B6energy, B7energy, B8energy, B9energy  
    real(realk) :: E_21, E_22, E_23, E_F12
    real(realk) :: tmp, energy, tmp2
    real(realk) :: temp
    real(realk) :: MP2energy

    ! ***********************************************************
    !   Check if we do the pair calculation
    ! ***********************************************************
    if((present(Fragment1) .AND. (.NOT. present(Fragment2))) .OR. &
         & (present(Fragment2) .AND. (.NOT. present(Fragment1)))) then
       call lsquit("get_f12_fragment_energy: Missing optional arguments Fragment1 and Fragment2")
    endif

    dopair = .FALSE.
    if(present(Fragment1) .AND. present(Fragment2)) then
       print *, "------------------------------"
       print *, " Do pair fragment calculation "
       print *, "------------------------------"
       dopair = .TRUE.
    endif

    if(.NOT. present(natoms) .AND. dopair) then
       call lsquit("get_f12_fragment_energy: Missing optional argument natoms")
    endif

    nbasis   = MyFragment%nbasis
    noccEOS  = MyFragment%noccEOS
    nunoccEOS = MyFragment%nunoccEOS
    noccfull = noccEOS

    nocvAOS = MyFragment%noccAOS + MyFragment%nunoccAOS
    noccAOS = MyFragment%noccAOS
    nvirtAOS = MyFragment%nunoccAOS

    ncabsAO = size(MyFragment%Ccabs,1)    
    ncabsMO = size(MyFragment%Ccabs,2)


    ! ***********************************************************
    !   Printing Input variables 
    ! ***********************************************************
    
    if(DECinfo%F12debug) then
       print *, "--------------------------"
       print *, "F12-integrals single fragment energy"
       print *, "--------------------------"
       print *, "nbasis: ", nbasis
       print *, "noccEOS: ", noccEOS
       print *, "nunoccEOS: ", nunoccEOS
       print *, "--------------------------"
       print *, "nocvAOS", nocvAOS
       print *, "noccAOS", noccAOS
       print *, "nvirtAOS", nvirtAOS
       print *, "ncabsAO", ncabsAO
       print *, "ncabsMO", ncabsMO
       print *, "--------------------------"
    end if

    ! ***********************************************************
    !   Allocating memory for V matrix
    ! ***********************************************************
    call mem_alloc(V1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)

    call mem_alloc(V2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)  
    call mem_alloc(Gijpq,  noccEOS, noccEOS, nocvAOS, nocvAOS)    
    call mem_alloc(Rijpq,  noccEOS, noccEOS, nocvAOS, nocvAOS)

    call mem_alloc(V3ijkl, noccEOS, noccEOS, noccEOS, noccEOS)  
    call mem_alloc(Gijmc,  noccEOS, noccEOS, noccAOS, ncabsMO)
    call mem_alloc(Rijmc,  noccEOS, noccEOS, noccAOS, ncabsMO)

    call mem_alloc(V4ijkl, noccEOS, noccEOS, noccEOS, noccEOS) 

    ! ***********************************************************
    !   Allocating memory for X matrix
    ! ***********************************************************
    call mem_alloc(X1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(X2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(X3ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(X4ijkl, noccEOS, noccEOS, noccEOS, noccEOS)

    ! ***********************************************************
    !   Allocating memory for B matrix
    ! ***********************************************************
    call mem_alloc(B1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)

    call mem_alloc(B2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(R2rlij, ncabsAO, noccEOS, noccEOS, noccEOS) 

    call mem_alloc(B3ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(R2ijkr,  noccEOS, noccEOS, noccEOS, ncabsAO)     

    call mem_alloc(B4ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijrs,  noccEOS, noccEOS, ncabsAO, ncabsAO)

    call mem_alloc(B5ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijrm,  noccEOS, noccEOS, ncabsAO, noccAOS)

    call mem_alloc(B6ijkl, noccEOS, noccEOS,  noccEOS, noccEOS)
    call mem_alloc(Rijpa,  noccEOS, noccEOS,  nocvAOS, nvirtAOS)

    call mem_alloc(B7ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijcm, noccEOS, noccEOS,  ncabsMO, noccAOS)

    call mem_alloc(B8ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijcr,  noccEOS, noccEOS, ncabsMO, ncabsAO) 

    call mem_alloc(B9ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijca, noccEOS, noccEOS, ncabsMO, nvirtAOS)

    ! Creating a CoccEOS matrix 
    call mem_alloc(CoccEOS, MyFragment%nbasis, noccEOS)
    do i=1, MyFragment%noccEOS
       ix = MyFragment%idxo(i)
       CoccEOS(:,i) = MyFragment%Co(:,ix)
    end do

    ! Creating a CoccAOS matrix 
    call mem_alloc(CoccAOS, MyFragment%nbasis, noccAOS)
    do i=1, MyFragment%noccAOS
       CoccAOS(:,i) = MyFragment%Co(:,i)
    end do

    ! Creating a CvirtAOS matrix 
    call mem_alloc(CvirtAOS, MyFragment%nbasis, nvirtAOS)
    do i=1, MyFragment%nunoccAOS
       CvirtAOS(:,i) = MyFragment%Cv(:,i)
    end do

    ! Creating a CocvAOS matrix 
    call mem_alloc(CocvAOS, MyFragment%nbasis, nocvAOS)
    do i=1, MyFragment%noccAOS
       CocvAOS(:,i) = MyFragment%Co(:,i)
    end do
    do i=1, MyFragment%nunoccAOS
       CocvAOS(:,i+MyFragment%noccAOS) = MyFragment%Cv(:,i)
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

    ! ***********************************************************
    ! Creating the C matrix 
    ! ***********************************************************
    call mem_alloc(Cijab, noccEOS, noccEOS, noccAOS, noccAOS)  

    ! ***********************************************************
    ! Creating the V matrix 
    ! ***********************************************************
    !> Get integrals <ij|f12*r^-1|kl> stored as (i,j,k,l)  (Note INTSPEC is always stored as (2,4,1,3) )      
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)

    !> V1ijrkl
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iiii','RRRRF',V1ijkl)

    !> Gijrp <ij|f12|pq> stored as (i,j,p,q)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iipp','RRRRC',Gijpq)

    !> Rijpq <ij|r^-1|pq> stored as (i,j,p,q)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iipp','RRRRG',Rijpq)

    !> Gijmc <ij|r^-1|mc> stored as (i,j,m,c) where c = cabs
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iimc','RCRRC',Gijmc)    

    !> Rijmc <ij|f12|mc> stored as (i,j,m,c) where c = cabs 
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iimc','RCRRG',Rijmc)

    m = noccEOS*noccEOS  ! <ij G pq> <pq R kl> = <m V2 n> 
    k = nocvAOS*nocvAOS  
    n = noccEOS*noccEOS

    !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','T',m,n,k,1.0E0_realk,Rijpq,m,Gijpq,n,0.0E0_realk,V2ijkl,m)

    m = noccEOS*noccEOS  ! <ij G mc> <mc R kl> = <m V3 n> 
    k = noccAOS*ncabsMO    ! m x k * k x n = m x n
    n = noccEOS*noccEOS

    !> dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','T',m,n,k,1.0E0_realk,Gijmc,m,Rijmc,n,0.0E0_realk,V3ijkl,m)

    !> Creating the V4ijkl = V3jilk !
    call array_reorder_4d(1.0E0_realk,V3ijkl,noccEOS,noccEOS,noccEOS,noccEOS,[2,1,4,3],0.0E0_realk,V4ijkl)

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, '            V matrix - Terms            '
       print *, '----------------------------------------'
       print *, 'norm4D(V1ijkl):', norm4D(V1ijkl)
       print *, '----------------------------------------'   
       print *, '(V2 Term):'
       print *, '----------------------------------------'   
       print *, 'norm4D(V2ijkl):', norm4D(V2ijkl)
       print *, 'norm4D(Gijpq):', norm4D(Gijpq)
       print *, 'norm4D(Rijpq):', norm4D(Rijpq)
       print *, '----------------------------------------'   
       print *, '(V3 Term):'
       print *, '----------------------------------------'   
       print *, 'norm4D(V3ijkl):', norm4D(V3ijkl)  
       print *, 'norm4D(Gijmc):', norm4D(Gijmc)
       print *, 'norm4D(Rijmc):', norm4D(Rijmc)
       print *, '----------------------------------------'   
       print *, '(V4 Term):'
       print *, '----------------------------------------'
       print *, 'norm4D(V4ijkl):', norm4D(V4ijkl)  
    end if

    ! ***********************************************************
    ! Creating the F matrix 
    ! ***********************************************************
    ! Creating a Fij MO matrix occ EOS
    call mem_alloc(Fij, noccEOS, noccEOS)
    Fij = 0E0_realk 
    do i=1, noccEOS
       do j=1, noccEOS      
          ix = MyFragment%idxo(i)
          iy = MyFragment%idxo(j)
          Fij(i,j) = MyFragment%ppfock(ix,iy)
       end do
    end do

    if(DECinfo%F12debug) then
       print *, "size(MyFragment%ppfock,1)", size(MyFragment%ppfock,1)
       print *, "size(MyFragment%ppfock,2)", size(MyFragment%ppfock,2)
       print *, "size(MyFragment%qqfock,1)", size(MyFragment%qqfock,1)
       print *, "size(MyFragment%qqfock,2)", size(MyFragment%qqfock,2)
    endif

    ! Creating a Fmn MO matrix occ AOS
    call mem_alloc(Fmn, noccAOS, noccAOS)
    Fmn = 0E0_realk 

    !> Double Storage! This need to be changed, just for conceptual reasons
    Fmn = MyFragment%ppfock

    ! Creating a Fab MO matrix virt AOS
    call mem_alloc(Fab, nvirtAOS, nvirtAOS)
    Fab = 0E0_realk 

    !> Double storage! This need to be changed, just for conceptual reasons
    Fab = MyFragment%qqfock

    ! ***********************************************************
    ! Creating the X matrix 
    ! ***********************************************************
    ! (Note INTSPEC is always stored as (2,4,1,3) )
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)

    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iiii','RRRR2',X1ijkl)

    m = noccEOS*noccEOS   ! <ij G pq> <pq R kl> = <m V3  n>
    k = nocvAOS*nocvAOS  
    n = noccEOS*noccEOS

    !> Creating the X2ijkl
    !>    dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','T',m,n,k,1.0E0_realk,Rijpq,m,Rijpq,n,0.0E0_realk,X2ijkl,m)

    m = noccEOS*noccEOS   ! <ij G mc> <mc R kl> = <m V3  n>
    k = noccAOS*ncabsMO
    n = noccEOS*noccEOS

    !> Creating the X3ijkl
    !>    dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','T',m,n,k,1.0E0_realk,Rijmc,m,Rijmc,n,0.0E0_realk,X3ijkl,m)

    !> Creating the X4ijkl = X4jilk 
    call array_reorder_4d(1.0E0_realk,X3ijkl,noccEOS,noccEOS,noccEOS,noccEOS,[2,1,4,3],0.0E0_realk,X4ijkl)

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, '          X matrix - Terms              '   
       print *, '----------------------------------------'
       print *, 'norm2D(Fij):'    , norm2D(Fij)
       print *, 'norm4D(X1ijkl):' , norm4D(X1ijkl)
       print *, 'norm4D(X2ijkl):' , norm4D(X2ijkl)
       print *, 'norm4D(X3ijkl):' , norm4D(X3ijkl)
       print *, 'norm4D(X4ijkl):' , norm4D(X4ijkl)
    end if

    ! ***********************************************************
    ! Creating the B matrix 
    ! ***********************************************************
    !> Get integral <ij|[[T,f12],f12]|kl> stored as (i,j,k,l) (Note INTSPEC is always stored as (2,4,1,3) )
    ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)

    !> B1-term
    !> B1ijkl <ij|[[T,f12],f12]|kl> stored as (i,j,k,l)
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iiii','RRRRD',B1ijkl)

    !> B2-term
    !> R2ijrk <ij|f12^2|rk> stored as (i,j,r,k)    r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'riii','RRCR2',R2rlij)

    !> B3-term
    !> R2ijkr <ij|f12^2|kr> stored as (i,j,k,r)    r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iiir','RCRR2',R2ijkr)

    !> B4-term
    !> R2ijrs <ij|f12|rs> stored as (i,j,r,s)      r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iirr','RCRCG',Rijrs)

    !> B5-term
    !> Rijrm <ij|f12|rm> stored as (i,j,r,m)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iirm','RRRCG',Rijrm)

    !> B6-term
    !> Rijpa <ij|f12|pa> stored as (i,j,p,a)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iipa','RRRRG',Rijpa)

    !> Rijcm <ij|f12|cm> stored as (i,j,c,m)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iicm','RRRCG',Rijcm)

    !> Rijcr <ij|f12|cr> stored as (i,j,c,r)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iicr','RCRCG',Rijcr)

    !> Rijca <ij|f12|ca> stored as (i,j,c,a)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iica','RRRCG',Rijca)


    ! ***********************************************************
    !                      B matrix 
    ! ***********************************************************

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, '      Individual B matrix - Terms       '   
       print *, '----------------------------------------'
       print *, '(B1 Term):'
       print *, '----------------------------------------'
       print *, 'norm4D(B1ijkl):', norm4D(B1ijkl)
       print *, '----------------------------------------'   
       print *, '(B2 Term):'
       print *, '----------------------------------------'
       print *, 'norm2D(hJir):',  norm2D(Myfragment%hJir)
       print *, 'norm4D(R2rlij):', norm4D(R2rlij)
       print *, '----------------------------------------'   
       print *, '(B3 Term):'
       print *, '----------------------------------------'
       print *, 'norm2D(hJir):',  norm2D(Myfragment%hJir)
       print *, 'norm4D(R2ijkr):', norm4D(R2ijkr)
       print *, '----------------------------------------'   
       print *, '(B4 Term):'
       print *, '----------------------------------------'
       print *, 'norm2D(Krs):',  norm2D(Myfragment%Krs)
       print *, 'norm4D(Rijrs):', norm4D(Rijrs)
       print *, '----------------------------------------'   
       print *, '(B5 Term):'
       print *, '----------------------------------------'
       print *, 'norm2D(Frs):',  norm2D(Myfragment%Frs)
       print *, 'norm4D(Rijrm):', norm4D(Rijrm)
       print *, '----------------------------------------'   
       print *, '(B6 Term):'
       print *, '----------------------------------------'
       print *, 'norm2D(Fij):',   norm2D(Fij)
       print *, 'norm2D(Fab):',   norm2D(Fab)
       print *, 'norm4D(Rijpa):', norm4D(Rijpa)
       print *, '----------------------------------------'   
       print *, '(B7 Term):'
       print *, '----------------------------------------'
       !print *, 'norm2D(Fnm):',  norm2D(Myfragment%Fnm)
       print *, 'norm4D(Rijcm):', norm4D(Rijcm)
       print *, '----------------------------------------'   
       print *, '(B8 Term):'
       print *, '----------------------------------------'
       print *, 'norm2D(Frm):',  norm2D(Myfragment%Frm)
       !print *, 'norm4D(Rijcm):', norm4D(Rijcm)
       print *, 'norm4D(Rijcr):', norm4D(Rijcr)
       print *, '----------------------------------------'   
       print *, '(B9 Term):'
       print *, '----------------------------------------'
       print *, 'norm2D(Fcp):',  norm2D(Myfragment%Fcp)
       !print *, 'norm4D(Rijpa):', norm4D(Rijpa)
       print *, 'norm4D(Rijca):', norm4D(Rijca)

    end if

    ! *************************************************
    !       Dgemms to get the different B terms       
    ! ************************************************

    !> term2
    !> B2ijkl
    B2ijkl = 0.0E0_realk
    m = noccEOS   ! <k h r> <rl R2 ij> = <kl B2  ij>    m k k n
    k = ncabsAO  
    n = noccEOS*noccEOS*noccEOS  
    !> dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    call dgemm('N','N',m,n,k,1.0E0_realk,Myfragment%hJir,m,R2rlij,k,0.0E0_realk,B2ijkl,m)

    !> term3
    !> B3ijkl Brute force
    B3ijkl = 0.0E0_realk
    do i=1, noccEOS
       do j=1, noccEOS
          tmp  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             tmp =   tmp + R2ijkr(i,j,i,r)*Myfragment%hJir(j,r)
             tmp2 = tmp2 + R2ijkr(i,i,j,r)*Myfragment%hJir(j,r) 
          enddo
          B3ijkl(i,j,i,j) = tmp
          B3ijkl(i,j,j,i) = tmp2
       enddo
    enddo

    !> term4
    !> B4ijkl Brute force
    B4ijkl = 0.0E0_realk
    do i=1, noccEOS
       do j=1, noccEOS

          tmp  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             do s=1, ncabsAO
                do t=1, ncabsAO
                   tmp = tmp + Rijrs(i,j,r,s)*Myfragment%Krs(t,s)*Rijrs(i,j,r,t) + &
                        & Rijrs(i,j,s,r)*Myfragment%Krs(t,s)*Rijrs(i,j,t,r) 

                   tmp2 = tmp2 + Rijrs(j,i,r,s)*Myfragment%Krs(t,s)*Rijrs(i,j,r,t) + &
                        & Rijrs(j,i,s,r)*Myfragment%Krs(t,s)*Rijrs(i,j,t,r) 
                enddo
             enddo
          enddo
          B4ijkl(i,j,i,j) = tmp
          B4ijkl(i,j,j,i) = tmp2
       enddo
    enddo

    !> term5
    !> B5ijkl Brute force with memory savings
    B5ijkl = 0.0E0_realk
    do i=1, noccEOS
       do j=1, noccEOS

          tmp  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             do s=1, ncabsAO
                do m=1, noccAOS
                   tmp = tmp + Rijrm(i,j,r,m)*Myfragment%Frs(s,r)*Rijrm(i,j,s,m) + &
                        & Rijrm(j,i,r,m)*Myfragment%Frs(s,r)*Rijrm(j,i,s,m) 

                   tmp2 = tmp2 + Rijrm(j,i,r,m)*Myfragment%Frs(s,r)*Rijrm(i,j,s,m) + &
                        & Rijrm(i,j,r,m)*Myfragment%Frs(s,r)*Rijrm(j,i,s,m)
                enddo
             enddo
          enddo
          B5ijkl(i,j,i,j) = tmp
          B5ijkl(i,j,j,i) = tmp2
       enddo
    enddo

    !> term6
    !> Need to change this and separate this into two parts one for the Fij and one for the Fab
    !> B6ijkl Brute force with memory savings
    B6ijkl = 0.0E0_realk
    do i=1, noccEOS
       do j=1, noccEOS    
          tmp  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do a=1, nvirtAOS

             do p=1, noccAOS
                do q=1, noccAOS
                   tmp = tmp + Rijpa(i,j,p,a)*Fmn(q,p)*Rijpa(i,j,q,a) + &
                        & Rijpa(j,i,p,a)*Fmn(q,p)*Rijpa(j,i,q,a) 

                   tmp2 = tmp2 + Rijpa(j,i,p,a)*Fmn(q,p)*Rijpa(i,j,q,a) + &
                        & Rijpa(i,j,q,a)*Fmn(q,p)*Rijpa(j,i,q,a) 
                enddo
             enddo

             do p=1, nvirtAOS
                do q=1, nvirtAOS
                   tmp = tmp + Rijpa(i,j,p+noccAOS,a)*Fab(q,p)*Rijpa(i,j,q+noccAOS,a) + &
                        & Rijpa(j,i,p+noccAOS,a)*Fab(q,p)*Rijpa(j,i,q+noccAOS,a) 

                   tmp2 = tmp2 + Rijpa(j,i,p+noccAOS,a)*Fab(q,p)*Rijpa(i,j,q+noccAOS,a) + &
                        & Rijpa(i,j,p+noccAOS,a)*Fab(q,p)*Rijpa(j,i,q+noccAOS,a) 
                enddo
             enddo

          enddo
          B6ijkl(i,j,i,j) = tmp
          B6ijkl(i,j,j,i) = tmp2
       enddo
    enddo

  
    !> B7ijkl Brute force with memory savings
    B7ijkl = 0.0E0_realk
    do i=1, noccEOS
       do j=1, noccEOS    
          tmp  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do c=1, ncabsMO
             do m=1, noccAOS
                do n=1, noccAOS
                   tmp = tmp + Rijcm(i,j,c,m)*Fmn(m,n)*Rijcm(i,j,c,n) + &
                        & Rijcm(j,i,c,m)*Fmn(m,n)*Rijcm(j,i,c,n) 

                   tmp2 = tmp2 + Rijcm(j,i,c,m)*Fmn(m,n)*Rijcm(i,j,c,n) + &
                        & Rijcm(i,j,c,m)*Fmn(m,n)*Rijcm(j,i,c,n) 
                enddo
             enddo
          enddo
          B7ijkl(i,j,i,j) = tmp
          B7ijkl(i,j,j,i) = tmp2
       enddo
    enddo

    !> B8ijkl Brute force with memory savings
    B8ijkl = 0.0E0_realk
    do i=1, noccEOS
       do j=1, noccEOS    
          tmp  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             do c=1, ncabsMO
                do m=1, noccAOS
                   tmp = tmp + Rijcm(i,j,c,m)*MyFragment%Frm(r,m)*Rijcr(i,j,c,r) + &
                        & Rijcm(j,i,c,m)*MyFragment%Frm(r,m)*Rijcr(j,i,c,r) 

                   tmp2 = tmp2 + Rijcm(j,i,c,m)*MyFragment%Frm(r,m)*Rijcr(i,j,c,r) + &
                        & Rijcm(i,j,c,m)*MyFragment%Frm(r,m)*Rijcr(j,i,c,r) 
                enddo
             enddo
          enddo
          B8ijkl(i,j,i,j) = tmp
          B8ijkl(i,j,j,i) = tmp2
       enddo
    enddo

    !> B9ijkl Brute force with memory savings
    B9ijkl = 0.0E0_realk

    do i=1, noccEOS
       do j=1, noccEOS    
          tmp  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do c=1, ncabsMO
             do a=1, nvirtAOS
                do p=1, nocvAOS
                   tmp = tmp + Rijpa(i,j,p,a)*MyFragment%Fcp(c,p)*Rijca(i,j,c,a) + &
                        & Rijpa(j,i,p,a)*MyFragment%Fcp(c,p)*Rijca(j,i,c,a)

                   tmp2 = tmp2 + Rijpa(j,i,p,a)*MyFragment%Fcp(c,p)*Rijca(i,j,c,a) + &
                        & Rijpa(i,j,p,a)*MyFragment%Fcp(c,p)*Rijca(j,i,c,a)
                enddo
             enddo
          enddo
          B9ijkl(i,j,i,j) = tmp
          B9ijkl(i,j,j,i) = tmp2
       !   print *, tmp, tmp2, B9ijkl(i,j,i,j), B9ijkl(i,j,j,i), i,j
       enddo
    enddo

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, '    B matrix - Terms for testing        '   
       print *, '----------------------------------------'
       print *, 'norm4D(B1ijkl):', norm4D(B1ijkl)
       print *, 'norm4D(B2ijkl):', norm4D(B2ijkl)
       print *, 'norm4D(B3ijkl):', norm4D(B3ijkl)
       print *, 'norm4D(B4ijij):', norm4D(B4ijkl)
       print *, 'norm4D(B5ijij):', norm4D(B5ijkl)
       print *, 'norm4D(B6ijij):', norm4D(B6ijkl)
       print *, 'norm4D(B7ijij):', norm4D(B7ijkl)
       print *, 'norm4D(B8ijij):', norm4D(B8ijkl)
       print *, 'norm4D(B9ijij):', norm4D(B9ijkl)
    end if

    if(dopair) then
       call get_mp2f12_pf_E21(V1ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V1energy, 1.0E0_realk)
       call get_mp2f12_pf_E21(V2ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V2energy, -1.0E0_realk)
       call get_mp2f12_pf_E21(V3ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V3energy, -1.0E0_realk)
       call get_mp2f12_pf_E21(V4ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V4energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V1ijkl, noccEOS, V1energy,  1.0E0_realk)
       call get_mp2f12_sf_E21(V2ijkl, noccEOS, V2energy, -1.0E0_realk)
       call get_mp2f12_sf_E21(V3ijkl, noccEOS, V3energy, -1.0E0_realk)
       call get_mp2f12_sf_E21(V4ijkl, noccEOS, V4energy, -1.0E0_realk)
    endif

    E_21 = V1energy + V2energy + V3energy + V4energy

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, '(Single Fragment Energies for V-matrix):'
       print *, '----------------------------------------'
       print *, "E_21_V_term1:", V1energy
       print *, "E_21_V_term2:", V2energy
       print *, "E_21_V_term3:", V3energy
       print *, "E_21_V_term4:", V4energy
       print *, '----------------------------------------'
       print *, "E_21_Vsum:", E_21
    end if

    if(dopair) then
       call get_mp2f12_pf_E22(Fij, X1ijkl, Fragment1, Fragment2, MyFragment, noccEOS, X1energy,  1.0E0_realk)
       call get_mp2f12_pf_E22(Fij, X2ijkl, Fragment1, Fragment2, MyFragment, noccEOS, X2energy, -1.0E0_realk)
       call get_mp2f12_pf_E22(Fij, X3ijkl, Fragment1, Fragment2, MyFragment, noccEOS, X3energy, -1.0E0_realk)
       call get_mp2f12_pf_E22(Fij, X4ijkl, Fragment1, Fragment2, MyFragment, noccEOS, X4energy, -1.0E0_realk)
    else    
       call get_mp2f12_sf_E22(Fij, X1ijkl, noccEOS, X1energy,  1.0E0_realk)
       call get_mp2f12_sf_E22(Fij, X2ijkl, noccEOS, X2energy, -1.0E0_realk)
       call get_mp2f12_sf_E22(Fij, X3ijkl, noccEOS, X3energy, -1.0E0_realk)
       call get_mp2f12_sf_E22(Fij, X4ijkl, noccEOS, X4energy, -1.0E0_realk)
    endif

    E_22 = X1energy + X2energy + X3energy + X4energy 
    
    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, '  E_22 X term (Single Fragment)         '
       print *, '----------------------------------------'
       print *, "E_22_X_term1:", X1energy
       print *, "E_22_X_term2:", X2energy
       print *, "E_22_X_term3:", X3energy
       print *, "E_22_X_term4:", X4energy
       print *, '----------------------------------------'
       print *, "E_22_Xsum:", E_22
    end if

    if(dopair) then
       call get_mp2f12_pf_E23(B1ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B1energy,  1.0E0_realk)
       call get_mp2f12_pf_E23(B2ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B2energy,  1.0E0_realk)
       call get_mp2f12_pf_E23(B3ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B3energy,  1.0E0_realk)
       call get_mp2f12_pf_E23(B4ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B4energy,  1.0E0_realk)
       call get_mp2f12_pf_E23(B5ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B5energy,  -1.0E0_realk)
       call get_mp2f12_pf_E23(B6ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B6energy,  -1.0E0_realk)
       call get_mp2f12_pf_E23(B7ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B7energy,  1.0E0_realk)
       call get_mp2f12_pf_E23(B8ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B8energy,  -2.0E0_realk)
       call get_mp2f12_pf_E23(B9ijkl, Fragment1, Fragment2, MyFragment, noccEOS, B9energy,  -2.0E0_realk)     
    else
       call get_mp2f12_sf_E23(B1ijkl, noccEOS, B1energy,  1.0E0_realk)
       call get_mp2f12_sf_E23(B2ijkl, noccEOS, B2energy,  1.0E0_realk)
       call get_mp2f12_sf_E23(B3ijkl, noccEOS, B3energy,  1.0E0_realk)
       call get_mp2f12_sf_E23(B4ijkl, noccEOS, B4energy,  1.0E0_realk)
       call get_mp2f12_sf_E23(B5ijkl, noccEOS, B5energy,  -1.0E0_realk)
       call get_mp2f12_sf_E23(B6ijkl, noccEOS, B6energy,  -1.0E0_realk)
       call get_mp2f12_sf_E23(B7ijkl, noccEOS, B7energy,  1.0E0_realk)
       call get_mp2f12_sf_E23(B8ijkl, noccEOS, B8energy,  -2.0E0_realk)
       call get_mp2f12_sf_E23(B9ijkl, noccEOS, B9energy,  -2.0E0_realk)     
    endif

    E_23 = B1energy + B2energy + B3energy + B4energy + B5energy + B6energy + B7energy &
         & + B8energy + B9energy

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, '  E_22 B term (Single Fragment)         '
       print *, '----------------------------------------'
       print *, "E_23_B_term1:", B1energy
       print *, "E_23_B_term2:", B2energy   
       print *, "E_23_B_term3:", B3energy   
       print *, "E_23_B_term4:", B4energy   
       print *, "E_23_B_term5:", B5energy   
       print *, "E_23_B_term6:", B6energy   
       print *, "E_23_B_term7:", B7energy   
       print *, "E_23_B_term8:", B8energy   
       print *, "E_23_B_term9:", B9energy   
       print *, '----------------------------------------'
       print *, "E_23_B_sum:", E_23
    end if

    E_F12 = E_21 + E_22 + E_23

    MP2energy = Myfragment%energies(FRAGMODEL_OCCMP2)

    if(.not. DECinfo%onlyoccpart) then
    end if
    
    if(DECinfo%F12debug) then
       print *,   '----------------- DEC-MP2F12 CALCULATION ----------------'
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: MP2 CORRELATION ENERGY = ', MP2energy
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY = ', E_21
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY = ', E_22
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY = ', E_23
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 CORRECTION TO ENERGY = ', E_F12
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: MP2-F12 CORRELATION ENERGY = ', MP2energy+E_F12
    end if

    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: MP2 CORRELATION ENERGY = ', MP2energy
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY = ', E_21
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY = ', E_22
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY = ', E_23
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 CORRECTION TO ENERGY = ', E_F12
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: MP2-F12 CORRELATION ENERGY = ', MP2energy+E_F12

    !> Setting the MP2-F12 correction
    Myfragment%energies(FRAGMODEL_MP2f12) = E_F12

    ! ***********************************************************
    ! Free Memory
    ! ***********************************************************

    !> Need to be free to avoid memory leak for the type(matrix) CMO_RI in CABS.F90
    ! call free_cabs()

    !> C-term
    call mem_dealloc(Fij)
    call mem_dealloc(Fmn)
    call mem_dealloc(Fab)
    call mem_dealloc(Cijab)

    !> Coeff
    call mem_dealloc(CoccEOS)
    call mem_dealloc(CoccAOS)
    call mem_dealloc(CocvAOS)
    call mem_dealloc(Ccabs)
    call mem_dealloc(Cri)
    call mem_dealloc(CvirtAOS)

    !> V-terms
    call mem_dealloc(V1ijkl)

    call mem_dealloc(V2ijkl)
    call mem_dealloc(Gijpq) 
    call mem_dealloc(Rijpq)

    call mem_dealloc(V3ijkl)
    call mem_dealloc(Gijmc)
    call mem_dealloc(Rijmc)

    call mem_dealloc(V4ijkl)

    !> X-terms
    call mem_dealloc(X1ijkl)
    call mem_dealloc(X2ijkl)
    call mem_dealloc(X3ijkl)
    call mem_dealloc(X4ijkl)

    !> B-terms
    call mem_dealloc(B1ijkl)

    call mem_dealloc(B2ijkl)
    call mem_dealloc(R2rlij)  

    call mem_dealloc(B3ijkl)
    call mem_dealloc(R2ijkr)     

    call mem_dealloc(B4ijkl)
    call mem_dealloc(Rijrs)

    call mem_dealloc(B5ijkl)
    call mem_dealloc(Rijrm)

    call mem_dealloc(B6ijkl)
    call mem_dealloc(Rijpa)

    call mem_dealloc(B7ijkl)
    call mem_dealloc(Rijcm)

    call mem_dealloc(B8ijkl)
    call mem_dealloc(Rijcr) 

    call mem_dealloc(B9ijkl)
    call mem_dealloc(Rijca)

  end subroutine get_f12_fragment_energy

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
    tmp = 0E0_realk         ! NB Important reset

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

    do j=1, nocc
       do i=j+1, nocc 
          tmp = tmp +  7.0E0_realk * Bijkl(i,j,i,j) + Bijkl(i,j,j,i)
       enddo
    enddo
    energy = energy + 0.0625E0_realk*tmp
    energy = energy*scalar

    call mem_dealloc(Bijkl)

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
    tmp = 0E0_realk

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
    type(decfrag),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(inout) :: Fragment2
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
    tmp = 0.0_realk

    call mem_alloc(dopair_occ,nocc,nocc)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)

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
  subroutine get_mp2f12_pf_E22(Fij, Xijkl, Fragment1, Fragment2, PairFragment, nocc, energy, scalar)
    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(inout) :: Fragment2
    !> Pair fragment formed from fragment 1 and 2
    type(decfrag), intent(inout) :: PairFragment

    real(realk),intent(in)  :: scalar
    real(realk),intent(out) :: energy
    real(realk),intent(in)  :: Fij(nocc,nocc)
    real(realk),intent(in)  :: Xijkl(nocc,nocc,nocc,nocc)
    integer,intent(in)      :: nocc
    !
    integer     :: i,j
    real(realk) :: tmp,tmp2
    logical,pointer :: dopair_occ(:,:)

    real(realk), pointer :: Bijkl(:,:,:,:)
    tmp = 0E0_realk

    call mem_alloc(Bijkl,nocc,nocc,nocc,nocc)
    call mem_alloc(dopair_occ,nocc,nocc)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)

    do j=1,nocc
       do i=1,nocc
          if(dopair_occ(i,j)) then !Do Pair 1 and 2   
             tmp2 = Fij(i,i) + Fij(j,j)
             Bijkl(i,j,i,j) = -1.0E0_realk*tmp2*Xijkl(i,j,i,j)
             Bijkl(i,j,j,i) = -1.0E0_realk*tmp2*Xijkl(i,j,j,i)
          endif
       enddo
    enddo

    energy = 0.25E0_realk*tmp
    tmp = 0E0_realk         ! NB Important reset

    do j=1, nocc
       do i=j+1, nocc 
          if(dopair_occ(i,j)) then !Do Pair 1 and 2   
             tmp = tmp +  7.0E0_realk * Bijkl(i,j,i,j) + Bijkl(i,j,j,i)
          endif
       enddo
    enddo
    energy = energy + 0.0625E0_realk*tmp
    energy = energy*scalar

    call mem_dealloc(Bijkl)
    call mem_dealloc(dopair_occ)

  end subroutine get_mp2f12_pf_E22

  subroutine get_mp2f12_pf_E23(ijkl, Fragment1, Fragment2, PairFragment, nocc, energy, scalar)
    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(inout) :: Fragment2
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
    tmp = 0E0_realk

    energy = 0.25E0_realk*tmp
    tmp = 0E0_realk         ! NB Important reset

    call mem_alloc(dopair_occ,nocc,nocc)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)

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

  !> Brief: Get <1,2|INTSPEC|3,4> MO integrals wrapper.
  !> Author: Yang M. Wang
  !> Data: Nov 2013
  subroutine get_mp2f12_MO(MyFragment,MySetting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,INTTYPE,INTSPEC,transformed_mo)
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
    !> Number of unoccupied (virtual) orbitals MO in EOS space
    integer :: nunoccEOS
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
    nunoccEOS = MyFragment%nunoccEOS
    nvirtAOS = MyFragment%nunoccAOS
    nocvAOS =   MyFragment%noccAOS + MyFragment%nunoccAOS
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
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
    type(batchtoorb), pointer :: batch2orbAlpha(:),batch2orbGamma(:)
#ifdef VAR_OMP
    integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
    TYPE(DECscreenITEM)   :: DecScreen

    integer :: MinAObatchSize, MaxAObatchSize, GammaBatchSize, AlphaBatchSize

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

    !> Minimum AO batch size
    call determine_maxBatchOrbitalsize(DECinfo%output,MySetting,MinAObatchSize,'R')

    !> Maximum AO batch size (all basis functions)
    MaxAObatchSize = n31
    !> Setting MinAO to AO batch size for debug purposes
    MinAObatchSize = n11

    !> Set alpha and gamma batch size as written above
    GammaBatchSize = n31 ! Needs to be changed, For DEBUG purposes MaxAObatchSize
    AlphaBatchSize = n11 ! Needs to be changes, For DEBUG purposes MinAObatchSize

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

    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
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

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************
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

    ! Setting to FALSE for DEBUG purposes
    Mysetting%scheme%cs_screen = .FALSE.
    Mysetting%scheme%ps_screen = .FALSE.

    ! Integral screening stuff
    doscreen = Mysetting%scheme%cs_screen .or. Mysetting%scheme%ps_screen
    call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mysetting,&
         & nbatchesAlpha,nbatchesGamma,INTSPEC)
    IF(doscreen)then
       call II_getBatchOrbitalScreen(DecScreen,mysetting,&
            & n31,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
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

          ! Get tmp1(beta(n21),delta(n41)|INTSPEC|alphaB(n11),gammaB(n31)) 
          ! ************************************************************************************
          dim1 = i8*n21*n41*dimAlpha*dimGamma   ! dimension for integral array tmp1
          call mem_alloc(tmp1,dim1)
          tmp1 = 0.0_realk

          IF(doscreen) mysetting%LST_GAB_RHS => DECSCREEN%masterGabRHS
          IF(doscreen) mysetting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p

          call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
               & mysetting, tmp1, batchindexAlpha(alphaB), batchindexGamma(gammaB), &
               & batchsizeAlpha(alphaB), batchsizeGamma(gammaB), n21, n41, dimAlpha, dimGamma, FullRHS,&
               & INTSPEC)
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
          
          do i=1, dim1
             if (abs(tmp1(i)) > 10.0E-10) then
               ! write(DECinfo%output,'(1X,a,1i4,f20.10)') 'tmp3:', i, tmp1(i)
             endif
          enddo
          
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
             if (abs(tmp2(i)) > 10.0E-10) then
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
             if (abs(tmp1(i)) > 10.0E-10) then
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
    nullify(mysetting%LST_GAB_LHS)
    nullify(mysetting%LST_GAB_RHS)
    call free_decscreen(DECSCREEN)

    call free_batch(orb2batchGamma, batchdimGamma, batchsizeGamma, batchindexGamma, batch2orbGamma, &
         & orb2batchAlpha, batchdimAlpha, batchsizeAlpha, batchindexAlpha, batch2orbAlpha, nbatchesGamma, nbatchesAlpha)

    ! Free F12 related pointers
    call mem_dealloc(C4T)
    call mem_dealloc(C3T)
    call mem_dealloc(C1T)
    call mem_dealloc(kjli)

  end subroutine get_mp2f12_AO_transform_MO

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

end module f12_integrals_module

