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
  use f12_routines_module!, only: MO_transform_AOMatrix, matrix_print, norm4D, norm2D, get_mp2f12_MO

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

  public :: get_f12_fragment_energy, matrix_print_4d, matrix_print_2d, get_mp2f12_sf_E21

  private

contains
  !> Brief: Gives the single fragment energy for MP2F12
  !> Author: Yang M. Wang
  !> Date: April 2013
  subroutine get_f12_fragment_energy(MyFragment, Taibj, Fragment1,Fragment2,natoms)
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

    !> The introduction of the C-terms as an external term
    !> F12 integrals for the V5_term sum_ab <ij|r^-1|ab> * <ab|t_2|kl>
    real(realk), pointer :: V5ijkl(:,:,:,:)  

    ! ***********************************************************
    !   Allocating for C matrix
    ! ***********************************************************
    !> Fock Fkj 
    real(realk), pointer :: Fkj(:,:)
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
    !> Fock Fac
    real(realk), pointer :: Rijac(:,:,:,:)

    !> t2EOS amplitudes stored in the order T(a,i,b,j)
    real(realk), intent(in), pointer :: Taibj(:,:,:,:) 
    real(realk), pointer :: Tabij(:,:,:,:) 
    
    ! ***********************************************************
    !   Allocating for X matrix Canonical
    ! ***********************************************************
    !> F12 integral X1ijmk
    real(realk), pointer :: X1ijkl(:,:,:,:)
    !> F12 integrals for the X2_term sum_pq <ij|g|pq> * <pq|g|kn>
    real(realk), pointer :: X2ijkl(:,:,:,:) 
    !> F12 integrals for the X3_term sum_pq <ij|g|mc> * <mc|g|kn>
    real(realk), pointer :: X3ijkl(:,:,:,:) 
    !> F12 integrals for the X4_term sum_pq <ij|g|cm> * <cm|g|kn>
    real(realk), pointer :: X4ijkl(:,:,:,:)    
    
    ! ***********************************************************
    !   Allocating for X matrix Non-Canonical
    ! ***********************************************************
    !> F12 integral X1ijkn
    real(realk), pointer :: X1ijkn(:,:,:,:)
    !> F12 integrals for the X2_term sum_pq <ij|f12|pq> * <pq|f12|kn>
    real(realk), pointer :: X2ijkn(:,:,:,:) 
    !real(realk), pointer :: Rijpq(:,:,:,:) <ij|f12|pq> * <pq|f12|kn>
    real(realk), pointer :: Rpqkn(:,:,:,:)  
    !> F12 integrals for the X3_term sum_pq <ij|f12|mc> * <mc|f12|kn>
    real(realk), pointer :: X3ijkn(:,:,:,:) 
    !> F12 integrals for the X4_term sum_pq <ij|f12|cm> * <cm|f12|kn>
    real(realk), pointer :: X4ijkn(:,:,:,:)  
    real(realk), pointer :: Rijcm(:,:,:,:)  
    real(realk), pointer :: Rcmkn(:,:,:,:)  
    
    real(realk), pointer :: X4ijnk(:,:,:,:)    
    real(realk), pointer :: Rmckn(:,:,:,:)  
    
    ! ***********************************************************
    !   Allocating for B matrix
    ! ***********************************************************
    !> F12 integrals for the B1_term <ij|[[T,f12],f12]|kl> 
    real(realk), pointer :: B1ijkl(:,:,:,:)   
    !> F12 integrals for the B2_term <ij|f12^2|rk>  r = RI MO   
    real(realk), pointer :: B2ijkl(:,:,:,:)     
    real(realk), pointer :: R2rlij(:,:,:,:)
    real(realk), pointer :: R2ijrk(:,:,:,:)
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
    !real(realk), pointer :: Rijcm(:,:,:,:)
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
    !> number of occupied MO orbitals in AOS 
    integer :: noccAOS
    !> number of virtual MO orbitals in AOS 
    integer :: nunoccAOS
    !> number of occupied + virtual MO orbitals in EOS 
    integer :: nocvAOS  
    !> number of virtual MO orbitals in AOS 
    integer :: nvirtAOS

    !> number of CABS AO orbitals
    integer :: ncabsAO
    !> number of CABS MO orbitals
    integer :: ncabsMO

    integer :: ix, iy, i, j, m, n, k, l, p, q, c, r, s, t, a, b

    real(realk) :: V1energy, V2energy, V3energy, V4energy, V5energy
    real(realk) :: X1energy, X2energy, X3energy, X4energy, X4energyY
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
       dopair = .TRUE.
    endif

    if(.NOT. present(natoms) .AND. dopair) then
       call lsquit("get_f12_fragment_energy: Missing optional argument natoms")
    endif

    nbasis   = MyFragment%nbasis
    noccEOS  = MyFragment%noccEOS
    nunoccEOS = MyFragment%nunoccEOS
    noccfull = noccEOS

    noccAOS = MyFragment%noccAOS
    nunoccAOS = MyFragment%nunoccAOS
    nocvAOS = MyFragment%noccAOS + MyFragment%nunoccAOS
    nvirtAOS = MyFragment%nunoccAOS

    ncabsAO = size(MyFragment%Ccabs,1)    
    ncabsMO = size(MyFragment%Ccabs,2)


    ! ***********************************************************
    !   Printing Input variables 
    ! ***********************************************************
    
    if(DECinfo%F12debug) then
       print *, "-------------------------------------------------"
       print *, "     F12-integrals.F90 single fragment energy    "
       print *, "-------------------------------------------------"
       print *, "nbasis:    ", nbasis
       print *, "noccEOS:   ", noccEOS
       print *, "nunoccEOS: ", nunoccEOS
       print *, "-------------------------------------------------"
       print *, "noccAOS    ", noccAOS
       print *, "nocvAOS    ", nocvAOS
       print *, "nvirtAOS   ", nvirtAOS
       print *, "ncabsAO    ", ncabsAO
       print *, "ncabsMO    ", ncabsMO
    end if

    ! ***********************************************************
    !   Allocating memory the C matrix 
    ! ***********************************************************
    call mem_alloc(Cijab, noccEOS,  noccEOS,  nvirtAOS,  nvirtAOS) 
    call mem_alloc(Rijac, noccEOS,  noccEOS,  nvirtAOS,  ncabsMO)
    call mem_alloc(Tabij, nvirtAOS, nvirtAOS,  noccEOS,  noccEOS)

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

    call mem_alloc(V5ijkl, noccEOS, noccEOS, noccEOS, noccEOS) 

    ! ***********************************************************
    !   Allocating memory for X matrix Canonical
    ! ***********************************************************
    call mem_alloc(X1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(X2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(X3ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(X4ijkl, noccEOS, noccEOS, noccEOS, noccEOS)

    ! ***********************************************************
    !   Allocating memory for X matrix Non-Canonical
    ! ***********************************************************
    call mem_alloc(X1ijkn, noccEOS, noccEOS, noccEOS, noccAOS)
    
    call mem_alloc(X2ijkn, noccEOS, noccEOS, noccEOS, noccAOS)
    call mem_alloc(Rpqkn,  nocvAOS, nocvAOS, noccEOS, noccAOS)
    
    call mem_alloc(X3ijkn, noccEOS, noccEOS, noccEOS, noccAOS)
    
    call mem_alloc(X4ijnk, noccEOS, noccEOS, noccAOS, noccEOS)
    call mem_alloc(Rmckn,  noccAOS, ncabsMO, noccEOS, noccAOS)
  
    call mem_alloc(X4ijkn, noccEOS, noccEOS, noccEOS, noccAOS)
    call mem_alloc(Rijcm,  noccEOS, noccEOS, ncabsMO, noccAOS)
    call mem_alloc(Rcmkn,  ncabsMO, noccAOS, noccEOS, noccAOS)
  
    ! ***********************************************************
    !   Allocating memory for B matrix
    ! ***********************************************************
    call mem_alloc(B1ijkl, noccEOS, noccEOS, noccEOS, noccEOS)

    call mem_alloc(B2ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(R2rlij, ncabsAO, noccEOS, noccEOS, noccEOS) 
    call mem_alloc(R2ijrk, noccEOS, noccEOS, ncabsAO, noccEOS)     
  
    call mem_alloc(B3ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(R2ijkr, noccEOS, noccEOS, noccEOS, ncabsAO)     

    call mem_alloc(B4ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijrs,  noccEOS, noccEOS, ncabsAO, ncabsAO)

    call mem_alloc(B5ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    call mem_alloc(Rijrm,  noccEOS, noccEOS, ncabsAO, noccAOS)

    call mem_alloc(B6ijkl, noccEOS, noccEOS,  noccEOS, noccEOS)
    call mem_alloc(Rijpa,  noccEOS, noccEOS,  nocvAOS, nvirtAOS)

    call mem_alloc(B7ijkl, noccEOS, noccEOS, noccEOS, noccEOS)
    !call mem_alloc(Rijcm, noccEOS, noccEOS,  ncabsMO, noccAOS)

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
    ! **********************************************************
    
    !> Rijca <ij|f12|ca> stored as (i,j,c,a)       r = RI MO
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iiac','RCRRG',Rijac)
    
    Cijab = 0.0E0_realk
    do i=1, noccEOS
       do j=1, noccEOS      
          do a=1, nvirtAOS
             do b=1, nvirtAOS      
    !            tmp2 = 0.0E0_realk
                tmp  = 0.0E0_realk
                do c=1, ncabsMO
                    tmp =   tmp + Rijac(i,j,a,c)*(Myfragment%Fcp(c,b+noccAOS)) + Rijac(j,i,b,c)*(Myfragment%Fcp(c,a+noccAOS)) 
     !              tmp2 = tmp2 + Rijac(j,i,a,c)*(Myfragment%Fcp(c,b+noccAOS)) + Rijac(i,j,b,c)*(Myfragment%Fcp(c,a+noccAOS)) 
                enddo
                Cijab(i,j,a,b) = tmp 
      !          Cijab(j,i,a,b) = tmp2 
             end do
          end do
       end do
    end do

    m = noccEOS*noccEOS    ! <ij C ab> <ab T kl> = <m V2 n> 
    k = nvirtAOS*nvirtAOS
    n = noccEOS*noccEOS
   ! call array_reorder_4d(1.0E0_realk,Taibj,nvirtAOS,nvirtAOS,noccEOS,noccEOS,[1,3,2,4],0.0E0_realk,Tabij)

   !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
   ! call dgemm('N','N',m,n,k,1.0E0_realk,Cijab,m,Tabij,n,0.0E0_realk,V5ijkl,m)

    !> Brute Force need to do this with dgemm at some point
    do i=1, noccEOS
       do j=1, noccEOS      
          do k=1, noccEOS
             do l=1, noccEOS      
                tmp  =  0.0E0_realk
!                tmp2  = 0.0E0_realk
                do a=1, nvirtAOS
                   do b=1, nvirtAOS      
                      tmp =  tmp +  Cijab(i,j,a,b)*Taibj(a,k,b,l) 
 !                     tmp2 = tmp2 + Cijab(j,i,a,b)*Taibj(a,k,b,l) 
                   enddo
                end do
                V5ijkl(i,j,k,l) = tmp 
           !     V5ijkl(j,i,k,l) = tmp2 
             end do
          end do
       end do
    end do

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

    !> Creating the V4ijkl = V3jilk 
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
       print *, 'norm4D(Gijpq):' , norm4D(Gijpq)
       print *, 'norm4D(Rijpq):' , norm4D(Rijpq)
       print *, '----------------------------------------'   
       print *, '(V3 Term):'
       print *, '----------------------------------------'   
       print *, 'norm4D(V3ijkl):', norm4D(V3ijkl)  
       print *, 'norm4D(Gijmc):' , norm4D(Gijmc)
       print *, 'norm4D(Rijmc):' , norm4D(Rijmc)
       print *, '----------------------------------------'   
       print *, '(V4 Term):'
       print *, '----------------------------------------'
       print *, 'norm4D(V4ijkl):', norm4D(V4ijkl)  
       print *, '----------------------------------------'
       print *, '(V5 Term):'
       print *, '----------------------------------------'
       print *, 'norm4D(Cijab):', norm4D(Cijab)
       print *, 'norm4D(Rijac):', norm4D(Rijac)
       print *, 'norm4D(Taibj):', norm4D(Taibj)
    end if

    if(dopair) then
       call get_mp2f12_pf_E21(V1ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V1energy, 1.0E0_realk)
       call get_mp2f12_pf_E21(V2ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V2energy, -1.0E0_realk)
       call get_mp2f12_pf_E21(V3ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V3energy, -1.0E0_realk)
       call get_mp2f12_pf_E21(V4ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V4energy, -1.0E0_realk)
       call get_mp2f12_pf_E21(V5ijkl, Fragment1, Fragment2, MyFragment, noccEOS, V5energy, -1.0E0_realk)
    else 
       call get_mp2f12_sf_E21(V1ijkl, noccEOS, V1energy,  1.0E0_realk)
       call get_mp2f12_sf_E21(V2ijkl, noccEOS, V2energy, -1.0E0_realk)
       call get_mp2f12_sf_E21(V3ijkl, noccEOS, V3energy, -1.0E0_realk)
       call get_mp2f12_sf_E21(V4ijkl, noccEOS, V4energy, -1.0E0_realk)
       call get_mp2f12_sf_E21(V5ijkl, noccEOS, V5energy,  1.0E0_realk)
    endif

    E_21 = 0.0E0_realk
    E_21 = V1energy + V2energy + V3energy + V4energy + V5energy

    if(DECinfo%F12debug) then
       print *, '----------------------------------------'
       print *, '(Single Fragment Energies for V-matrix):'
       print *, '----------------------------------------'
       print *, "E_21_V_term1:", V1energy
       print *, "E_21_V_term2:", V2energy
       print *, "E_21_V_term3:", V3energy
       print *, "E_21_V_term4:", V4energy
       print *, "E_21_V_term5:", V5energy
       print *, '----------------------------------------'
       print *, "E_21_Vsum:", E_21
    end if
    

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

    ! Creating a Fij MO matrix occ EOS
    call mem_alloc(Fij, noccEOS, noccEOS)
    do i=1, noccEOS
       do j=1, noccEOS      
          ix = MyFragment%idxo(i)
          iy = MyFragment%idxo(j)
          Fij(i,j) = MyFragment%ppfock(ix,iy)
       end do
    end do

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

    if(DECinfo%use_canonical) then

       ! ***********************************************************
       ! Creating the X matrix Canonical
       ! ***********************************************************
       ! (Note INTSPEC is always stored as (2,4,1,3) )
       ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)

       call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iiim','RRRR2',X1ijkl)
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
          print *, '     X matrix - Terms Canonical         '   
          print *, '----------------------------------------'
          print *, 'norm2D(Fij):'    , norm2D(Fij)
          print *, 'norm4D(X1ijkl):' , norm4D(X1ijkl)
          print *, 'norm4D(X2ijkl):' , norm4D(X2ijkl)
          print *, 'norm4D(X3ijkl):' , norm4D(X3ijkl)
          print *, 'norm4D(X4ijkl):' , norm4D(X4ijkl)
       endif

    else !non-canonical

       ! ***********************************************************
       ! Creating the X matrix Non-Canonical
       ! ***********************************************************
       ! (Note INTSPEC is always stored as (2,4,1,3))
       ! (beta,delta,alpha,gamma) (n2,n4,n1,n3)

       call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iiim','RRRR2',X1ijkn)
       call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'ppim','RRRRG',Rpqkn)
       call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'mcim','CRRRG',Rmckn)
       call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iicm','RRRCG',Rijcm)
       call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'cmim','RRCRG',Rcmkn)

       m = noccEOS*noccEOS   ! <ij R pq> <pq R kn> = <m X2  n>
       k = nocvAOS*nocvAOS  
       n = noccEOS*noccAOS

       !> Creating the X2ijkn
       !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       call dgemm('N','N',m,n,k,1.0E0_realk,Rijpq,m,Rpqkn,k,0.0E0_realk,X2ijkn,m)

       m = noccEOS*noccEOS   ! <ij R mc> <mc R kn> = <m X3  n>
       k = noccAOS*ncabsMO
       n = noccEOS*noccAOS

       !> Creating the X3ijkn
       !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       call dgemm('N','N',m,n,k,1.0E0_realk,Rijmc,m,Rmckn,k,0.0E0_realk,X3ijkn,m)

       !> Creating the X3ijkn = X3jink = X4ijnk
       call array_reorder_4d(1.0E0_realk,X3ijkn,noccEOS,noccEOS,noccEOS,noccAOS,[2,1,4,3],0.0E0_realk,X4ijnk)

       !> Creating the X4ijkn
       !>  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       m = noccEOS*noccEOS   ! <ij R mc> <mc R nk> = <m X4  n> instead of <ij R cm> <cm R kn> = <m X4  n>
       k = noccAOS*ncabsMO
       n = noccEOS*noccAOS

       call dgemm('N','N',m,n,k,1.0E0_realk,Rijcm,m,Rcmkn,k,0.0E0_realk,X4ijkn,m)

       if(DECinfo%F12debug) then
          print *, '----------------------------------------'
          print *, '   X matrix - Terms Non-Canonical       '   
          print *, '----------------------------------------'
          print *, 'norm2D(Fkj):   ' , norm2D(Fkj)
          print *, 'norm4D(X1ijkn):' , norm4D(X1ijkn)
          print *, 'norm4D(X2ijkn):' , norm4D(X2ijkn)
          print *, 'norm4D(X3ijkn):' , norm4D(X3ijkn)
          print *, 'norm4D(X4ijkn):' , norm4D(X4ijkn)
          print *, 'norm4D(X4ijnk):' , norm4D(X4ijnk)
       end if

    endif

    X1energy = 0.0E0_realk
    X2energy = 0.0E0_realk
    X3energy = 0.0E0_realk
    X4energy = 0.0E0_realk

    if(DECinfo%use_canonical) then

       if(dopair) then
          call get_mp2f12_pf_E22(Fij, X1ijkl, noccEOS, noccEOS, Fragment1, Fragment2, MyFragment, X1energy,  1.0E0_realk)
          call get_mp2f12_pf_E22(Fij, X2ijkl, noccEOS, noccEOS, Fragment1, Fragment2, MyFragment, X2energy, -1.0E0_realk)
          call get_mp2f12_pf_E22(Fij, X3ijkl, noccEOS, noccEOS, Fragment1, Fragment2, MyFragment, X3energy, -1.0E0_realk)
          call get_mp2f12_pf_E22(Fij, X4ijkl, noccEOS, noccEOS, Fragment1, Fragment2, MyFragment, X4energy, -1.0E0_realk)
       else    
          call get_mp2f12_sf_E22(Fij, X1ijkl, noccEOS, noccAOS, X1energy,  1.0E0_realk)
          call get_mp2f12_sf_E22(Fij, X2ijkl, noccEOS, noccAOS, X2energy, -1.0E0_realk)
          call get_mp2f12_sf_E22(Fij, X3ijkl, noccEOS, noccAOS, X3energy, -1.0E0_realk)
          call get_mp2f12_sf_E22(Fij, X4ijkl, noccEOS, noccAOS, X4energy, -1.0E0_realk)
       endif
   
   else !> Noncanonical

       if(dopair) then
          call get_mp2f12_pf_E22(Fkj, X1ijkn, noccEOS, noccAOS, Fragment1, Fragment2, MyFragment, X1energy,  1.0E0_realk)
          call get_mp2f12_pf_E22(Fkj, X2ijkn, noccEOS, noccAOS, Fragment1, Fragment2, MyFragment, X2energy, -1.0E0_realk)
          call get_mp2f12_pf_E22(Fkj, X3ijkn, noccEOS, noccAOS, Fragment1, Fragment2, MyFragment, X3energy, -1.0E0_realk)
          call get_mp2f12_pf_E22(Fkj, X4ijkn, noccEOS, noccAOS, Fragment1, Fragment2, MyFragment, X4energy, -1.0E0_realk)
       else    
          call get_mp2f12_sf_E22(Fkj, X1ijkn, noccEOS, noccAOS, X1energy,  1.0E0_realk)
          call get_mp2f12_sf_E22(Fkj, X2ijkn, noccEOS, noccAOS, X2energy, -1.0E0_realk)
          call get_mp2f12_sf_E22(Fkj, X3ijkn, noccEOS, noccAOS, X3energy, -1.0E0_realk)
          call get_mp2f12_sf_E22(Fkj, X4ijkn, noccEOS, noccAOS, X4energy, -1.0E0_realk)
       endif
  
    endif

    E_22 = 0.0E0_realk
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
    call get_mp2f12_MO(MyFragment,MyFragment%MyLsitem%Setting,CoccEOS,CoccAOS,CocvAOS,Ccabs,Cri,CvirtAOS,'iiri','RRRC2',R2ijrk)

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
       print *, 'norm4D(R2ijrk):', norm4D(R2ijrk)
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
       print *, 'norm2D(Fmn):',   norm2D(Fmn)
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

!!$    !> term2
!!$    !> B2ijkl
!!$    B2ijkl = 0.0E0_realk
!!$    m = noccEOS   ! <k h r> <rl R2 ij> = <kl B2  ij>    m k k n
!!$    k = ncabsAO  
!!$    n = noccEOS*noccEOS*noccEOS  
!!$    !> dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!$    call dgemm('N','N',m,n,k,1.0E0_realk,Myfragment%hJir,m,R2rlij,k,0.0E0_realk,B2ijkl,m)

    !> term2
    !> B2ijkl Brute force
    ! 4 Do loop and OMP setting this to zero (Not necessary in general)
    B2ijkl = 0.0E0_realk 
    B3ijkl = 0.0E0_realk
    B4ijkl = 0.0E0_realk
    B5ijkl = 0.0E0_realk
    B6ijkl = 0.0E0_realk
    B7ijkl = 0.0E0_realk
    B8ijkl = 0.0E0_realk
    B9ijkl = 0.0E0_realk
   !$OMP PARALLEL PRIVATE(i,j,r,s,t,m,c,n,a,p,q, &
     !$OMP tmp,tmp2) DEFAULT(shared)
    
    !$OMP DO COLLAPSE(2) 
    do i=1, noccEOS
       do j=1, noccEOS
          tmp  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             tmp =   tmp +  R2ijrk(i,j,r,j)*Myfragment%hJir(i,r)
             tmp2 =  tmp2 + R2ijrk(j,i,r,j)*Myfragment%hJir(i,r) 
          enddo
          B2ijkl(i,j,i,j) = tmp
          B2ijkl(i,j,j,i) = tmp2
       enddo
    enddo
    !$OMP END DO NOWAIT

    !> term3
    !> B3ijkl Brute force

    !$OMP DO COLLAPSE(2)
    do i=1, noccEOS
       do j=1, noccEOS
          tmp  = 0.0E0_realk
          tmp2 = 0.0E0_realk
          do r=1, ncabsAO
             tmp =   tmp +  R2ijkr(i,j,i,r)*Myfragment%hJir(j,r)
             tmp2 =  tmp2 + R2ijkr(i,i,j,r)*Myfragment%hJir(j,r) 
          enddo
          B3ijkl(i,j,i,j) = tmp
          B3ijkl(i,j,j,i) = tmp2
       enddo
    enddo
    !$OMP END DO NOWAIT


    !> term4
    !> B4ijkl Brute force

    !$OMP DO COLLAPSE(2)
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
    !$OMP END DO NOWAIT

    !> term5
    !> B5ijkl Brute force with memory savings
    !$OMP DO COLLAPSE(2)   
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
    !$OMP END DO NOWAIT

    
    !> term6
    !> Need to change this and separate this into two parts one for the Fij and one for the Fab
    !> B6ijkl Brute force with memory savings
    
    !$OMP DO COLLAPSE(2)
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
                        & Rijpa(i,j,p,a)*Fmn(q,p)*Rijpa(j,i,q,a) 

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
    !$OMP END DO NOWAIT

    !> B7ijkl Brute force with memory savings
    !$OMP DO COLLAPSE(2)
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
    !$OMD END DO NOWAIT

    !> B8ijkl Brute force with memory savings
    !$OMP DO COLLAPSE(2)
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
    !$OMP END DO NOWAIT 

    !> B9ijkl Brute force with memory savings
    !$OMP DO COLLAPSE(2)
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
    !$OMP END DO NOWAIT 
 
    !$OMP END PARALLEL

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

    E_23 = 0.0E0_realk
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
    
    E_F12 = 0.0E0_realk
    E_F12 = E_21 + E_22 + E_23

    MP2energy = Myfragment%energies(FRAGMODEL_OCCMP2)

    if(.not. DECinfo%onlyoccpart) then
    end if
    
    if(DECinfo%F12debug) then
       print *,   '----------------------------------------------------------------'
       print *,   '                    DEC-MP2F12 CALCULATION                      '
       print *,   '----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: MP2 CORRELATION ENERGY =           ', MP2energy
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E_21
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E_22
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E_23
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 CORRECTION TO ENERGY =         ', E_F12
       write(*,'(1X,a,f20.10)') 'WANGY TOYCODE: MP2-F12 CORRELATION ENERGY =       ', MP2energy+E_F12
    end if

    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: MP2 CORRELATION ENERGY =           ', MP2energy
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E_21
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E_22
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E_23
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: F12 CORRECTION TO ENERGY =         ', E_F12
    write(DECinfo%output,'(1X,a,f20.10)') 'WANGY TOYCODE: MP2-F12 CORRELATION ENERGY =       ', MP2energy+E_F12

    !> Setting the MP2-F12 correction
    Myfragment%energies(FRAGMODEL_MP2f12) = E_F12

    !> Need to be set for the single_fragments
    Myfragment%EoccFOP_Corr = E_F12

    ! ***********************************************************
    ! Free Memory
    ! ***********************************************************

    !> Need to be free to avoid memory leak for the type(matrix) CMO_RI in CABS.F90
    ! call free_cabs()

    !> F-term
    call mem_dealloc(Fij)
    call mem_dealloc(Fkj)
    call mem_dealloc(Fmn)
    call mem_dealloc(Fab)
   
    !> C-term
    call mem_dealloc(Cijab)
    call mem_dealloc(Rijac)
    call mem_dealloc(Tabij)

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

    call mem_dealloc(V5ijkl)

    !> X-terms - Canonical
    call mem_dealloc(X1ijkl)
    call mem_dealloc(X2ijkl)
    call mem_dealloc(X3ijkl)
    call mem_dealloc(X4ijkl)

    !> X-terms - Non-Canonical
    call mem_dealloc(X1ijkn)
   
    call mem_dealloc(X2ijkn)
    call mem_dealloc(Rpqkn)
    
    call mem_dealloc(X3ijkn)
     
    call mem_dealloc(X4ijnk)
    call mem_dealloc(Rmckn)
     
    call mem_dealloc(X4ijkn)
    call mem_dealloc(Rcmkn)
   
    !> B-terms
    call mem_dealloc(B1ijkl)

    call mem_dealloc(B2ijkl)
    call mem_dealloc(R2rlij)  
    call mem_dealloc(R2ijrk)     
    
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
  subroutine get_mp2f12_sf_E22(Fij, Xijkl, n1, n2, energy, scalar)
    implicit none

    real(realk),intent(in)  :: scalar
    real(realk),intent(out) :: energy
    integer,intent(in)      :: n1
    integer,intent(in)      :: n2
    real(realk),intent(in)  :: Fij(n2,n1)
    real(realk),intent(in)  :: Xijkl(n1,n1,n1,n2)
    !
    integer     :: i,j,k
    real(realk) :: tmp,tmp2,tmp3

    real(realk), pointer :: Bijkl(:,:,:,:)
    tmp = 0E0_realk
    tmp2 = 0E0_realk
    tmp3 = 0E0_realk

    call mem_alloc(Bijkl,n1,n1,n1,n1)

    Bijkl = 0.0E0_realk

    if(DECinfo%use_canonical) then
       do j=1,n1
          do i=1,n1
             tmp2 = Fij(i,i) + Fij(j,j)
             Bijkl(i,j,i,j) = -1.0E0_realk*tmp2*Xijkl(i,j,i,j)
             Bijkl(i,j,j,i) = -1.0E0_realk*tmp2*Xijkl(i,j,j,i)
          enddo
       enddo

    else !> Non-canonical
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k,tmp2,tmp3) DEFAULT(shared)
         do i=1,n1
          do j=1,n1
             tmp2 = 0E0_realk
             tmp3 = 0E0_realk
             do k=1,n2
                tmp2 = tmp2 + Xijkl(i,j,i,k)*Fij(k,j) + Xijkl(j,i,j,k)*Fij(k,i)  
                tmp3 = tmp3 + Xijkl(j,i,i,k)*Fij(k,j) + Xijkl(i,j,j,k)*Fij(k,i)  
             enddo
             Bijkl(i,j,i,j) = -1.0E0_realk*tmp2
             Bijkl(i,j,j,i) = -1.0E0_realk*tmp3

          enddo
       enddo
      !$OMP END PARALLEL DO
    endif

    do i=1, n1
       tmp = tmp + Bijkl(i,i,i,i)
    enddo
    energy = 0.25E0_realk*tmp


    tmp2 = 0E0_realk

    !Spawning of threads making unique copies of i,j and tmp on different
    !memory spaces.
    !$OMP PARALLEL PRIVATE(i,j,tmp) DEFAULT(shared)
    tmp  = 0E0_realk         ! NB Important reset
    !$OMP DO 
    do j=1, n1
       do i=j+1, n1
          tmp = tmp +  7.0E0_realk * Bijkl(i,j,i,j) + Bijkl(i,j,j,i)
       enddo
    enddo
    !The first thread finished can go on
    !$OMP END DO NOWAIT

    !This needs to be done in serial, they cannot read and write simultaneously
    !$OMP CRITICAL
    tmp2 = tmp2 + tmp
    !$OMP END CRITICAL
    !$OMP END PARALLEL  
    energy = energy + 0.0625E0_realk*tmp2
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
  subroutine get_mp2f12_pf_E22(Fkj, Xijkl, n1, n2, Fragment1, Fragment2, PairFragment, energy, scalar)
    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(inout) :: Fragment2
    !> Pair fragment formed from fragment 1 and 2
    type(decfrag), intent(inout) :: PairFragment

    real(realk),intent(in)  :: scalar
    real(realk),intent(out) :: energy
    integer,intent(in)  :: n1
    integer,intent(in)  :: n2
    real(realk),intent(in)  :: Fkj(n2,n1)
    real(realk),intent(in)  :: Xijkl(n1,n1,n1,n2)
    !
    integer     :: i,j,k
    real(realk) :: tmp,tmp2,tmp3
    logical,pointer :: dopair_occ(:,:)

    real(realk), pointer :: Bijkl(:,:,:,:)

    tmp = 0E0_realk
    tmp2 = 0E0_realk
    tmp3 = 0E0_realk
    
    call mem_alloc(Bijkl,n1,n1,n1,n1)
    call mem_alloc(dopair_occ,n1,n1)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)

    Bijkl = 0.0E0_realk

    if(DECinfo%use_canonical) then
       do j=1,n1
          do i=1,n1
             if(dopair_occ(i,j)) then !Do Pair 1 and 2   
                tmp2 = Fkj(i,i) + Fkj(j,j)
                Bijkl(i,j,i,j) = -1.0E0_realk*tmp2*Xijkl(i,j,i,j)
                Bijkl(i,j,j,i) = -1.0E0_realk*tmp2*Xijkl(i,j,j,i)
             endif
          enddo
       enddo

    else 
       do i=1,n1
          do j=1,n1
             tmp2 = 0E0_realk
             tmp3 = 0E0_realk
             do k=1,n2
                tmp2 = tmp2 + Xijkl(i,j,i,k)*Fkj(k,j) + Xijkl(j,i,j,k)*Fkj(k,i)  
                tmp3 = tmp3 + Xijkl(j,i,i,k)*Fkj(k,j) + Xijkl(i,j,j,k)*Fkj(k,i)  
             enddo
             Bijkl(i,j,i,j) = -1.0E0_realk*tmp2
             Bijkl(i,j,j,i) = -1.0E0_realk*tmp3
          enddo
       enddo
    endif

    energy = 0.25E0_realk*tmp
    tmp = 0E0_realk         ! NB Important reset

    do j=1, n1
       do i=j+1, n1 
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


end module f12_integrals_module

