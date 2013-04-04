!GIVE MODULE NAME
MODULE pbc_setup
  USE precision
  USE fundamental
  USE TYPEDEF
  USE matrix_module
  USE lattice_vectors
  USE lattice_type
  use pbc_msc
  USE memory_handling
  USE multipole_pbc
  USE harmonics_pbc
  USE pbc_matrix_operations
  USE pbc_scfdiis
  USE pbc_interactions
  USE pbc_ff_contrib
contains

!Subroutine to copy atoms and interface it with the integration code developed
!by Simen Reine.
SUBROUTINE set_pbc_molecules(INPUT,SETTING,lupri,luerr,nbast,Dmat,ll)
!  use pbc_matrix_module
! PROPERTIES SECTION
!  use lsdalton_rsp_mod
  implicit none
  !INPUT VARIABLES
  INTEGER, INTENT(IN) :: lupri, luerr, nbast
  TYPE(lvec_list_t), INTENT(INOUT) :: ll
  TYPE(DALTONINPUT),INTENT(INOUT) :: INPUT
  TYPE(LSSETTING):: SETTING
  TYPE(MATRIX),intent(INOUT)  ::  Dmat
  !Local variables
  TYPE(MOLECULEINFO) :: refcell
  TYPE(MOLECULEINFO),pointer :: latt_cell(:)
  TYPE(MATRIX)  ::  D(1),D1(1)!,F(1)
  TYPE(matrix),pointer :: nfdensity(:),f_1(:),ovl(:)
  TYPE(matrix),pointer :: g_2(:)
  !TYPE(lattice_cell_info_t),pointer :: sphermom(:)
  REAL(realk) :: Dfull(nbast,nbast),kvec(3,3),origin(3)
  complex(complexk), pointer :: fck(:), S_ab(:),k_fock(:,:),K_Sab(:,:)
  real(realk), pointer :: Tlat(:,:)
  REAL(realk)::  latt_vec_std(3),focknorm
  REAL(realk)::  E_cell,E_kin,E_ff,E_XC,E_K,E_J,E_en,E_nuc
  REAL(realk)::  E_1,E_nn
  real(realk) :: TS,TE
!  REAL(realk) :: PI=3.14159265358979323846D0
  INTEGER :: sze,num_latvectors
  INTEGER :: maxmultmom,n1,n2,n3,nfsze,fdim(3),Tlmax
  INTEGER :: i,j,k,scfit,iunit,l1,l2,l3,nbasterik,ierror
  character*(20) :: mattxt,string1,numtostring1,numtostring2,numtostring3
  TYPE(BZgrid_t) :: BZ


!!!!!!!!!!!!Information for me
!!!!! input%Basis%regular%atomtype(:)%shell(:)%segment(1)%exponents(:)
write(lupri,*) 'shells ',input%Basis%regular%atomtype(1)%nangmom
write(lupri,*) 'Number of atom species ',input%Basis%regular%natomtypes
write(lupri,*) 'Exponents ',(input%Basis%regular%atomtype(1)%shell(1)%segment(1)%exponents)
!  call screening_ovl(input%Basis%regular)
  call build_lvec_list(ll,nbast) 

  pbc_control%ldef%is_active=ll%ldef%is_active
  nearfield=ll%nf
  write(*,*) '1 active not active:',ll%ldef%is_active(1)
  write(*,*) '2 active not active:',ll%ldef%is_active(2)
  write(*,*) '3 active not active:',ll%ldef%is_active(3)

  call set_refcell(refcell,input%molecule)
  
 ! write(lupri,*) ll%ldef%avec
  origin=1.0d0
 !origin of cell should be 0d0, that means that I have to change the origin for
 !the next cell to what it should be and the atomic position in
 !each cell should be according to the origin of the new cell, a copy of
 !the reference cell.

  
  write(*,*) 'number of layers', ll%max_layer
  call latt_2_std_coord(origin,latt_vec_std,ll%ldef%avec)
  
  num_latvectors=size(ll%lvec)
  write(lupri,*) 'Number of vectors ', num_latvectors

  allocate(latt_cell(num_latvectors))

  call set_lattice_cells(latt_cell,num_latvectors,input%molecule,ll,lupri)


  write(lupri,*) 'setting%p%atom%center', latt_cell(2)%atom(1)%center(1)


  write(lupri,*) 'cutoff' 
  call find_cutoff_onep(lupri,luerr,setting,nbast,ll,&
              latt_cell, refcell,num_latvectors)

  call build_nflvec_list(ll,nbast) 
  n_neighbour=ll%nneighbour

 
  SELECT CASE(ll%wannier_direct)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BLOCH FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CASE('indirectly')

    write(lupri,*) 'We are now in the ',ll%wannier_direct, ' method'
    write(lupri,*) 'CS screen:', setting%scheme%CS_SCREEN 
    write(lupri,*) 'PS screen:', setting%scheme%PS_SCREEN 
    write(*,*) 'CS screen:', setting%scheme%CS_SCREEN 
    write(*,*) 'PS screen:', setting%scheme%PS_SCREEN 

   ! allocate(S_ab(nbast*nbast))
   ! S_ab=0.0000
  
   ! allocate(fck(nbast*nbast))
   ! fck=0.0000E-100
  
!
    call pbc_init_recvec(ll%ldef%avec,kvec)


    write(lupri,*) 'kvec', kvec!(3,3)
    write(lupri,*) 'a1 dot b1 = ',dot_product(ll%ldef%avec(1,:),kvec(1,:))
    write(lupri,*) 'lat vec', ll%ldef%avec!(3,3)
    write(lupri,*) 'module test: ', mod(3,3)
    !STOP

    call pbc_init_Bzgrid(kvec,ll%nk1,ll%nk2,ll%nk3,'nosym',bZ,nbast,nbast,lupri)
    lat_data%num_k1=ll%nk1
    lat_data%num_k2=ll%nk2
    lat_data%num_k3=ll%nk3
    lat_data%num_kpoints=bz%nk 
    lat_data%reclatvec(:,:)=kvec(:,:)

  kvec=0.0d0

  write(lupri,*) 'kvec(3,3)', kvec(3,3)

  write(*,*) 'll%nneighbour', ll%nneighbour
  write(*,*) 'll%nf', ll%nf

  call pbc_get_nfsize(n1,n2,n3,ll%nneighbour,lupri)
  nfsze=(2*n1+1)*(2*n2+1)*(2*n3+1)
  write(*,*) 'nfsize: ',nfsze
  write(lupri,*) 'nfsize: ',nfsze
  call mem_alloc(nfdensity,num_latvectors)
  call mem_alloc(f_1,num_latvectors)
  call mem_alloc(Ovl,num_latvectors)
  call mem_alloc(g_2,num_latvectors)
    k=0
    l1=0
    l2=0
    l3=0
    !do l1=-3,3

    if(ll%testcase) THEN !THIS IS MAINLY FOR DEBUGGING
      iunit = 345
      scfit=1
    !  write(string1,'(I5)')  scfit
    !  string1=adjustl(string1)
    !  write(numtostring1,'(I5)')  l1
    !  write(numtostring2,'(I5)')  l2
    !  write(numtostring3,'(I5)')  l3
    !  numtostring1=adjustl(numtostring1)
    !  numtostring2=adjustl(numtostring2)
    !  numtostring3=adjustl(numtostring3)

    !  write(*,*) string1,numtostring1,numtostring2,numtostring3

    !  write(mattxt,'(A20)') 'PBCDMAT'//trim(string1)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)


      call find_latt_index(k,0,0,0,fdim,ll,ll%max_layer)
      !do i=1,num_latvectors
      !  call init_lvec_data(ll%lvec(i),nbast)
      !enddo
      if(.not.ll%read_file) then
        ll%lvec(k)%d_mat(1,1)=0.18197943668877323D0
        ll%lvec(k)%d_mat(1,2)=0.28409190914049431D0
        ll%lvec(k)%d_mat(1,3)=0.18197943668877323D0
        ll%lvec(k)%d_mat(1,4)=0.28409190914049431D0
        ll%lvec(k)%d_mat(2,1)=0.28409190914049431D0
        ll%lvec(k)%d_mat(2,2)=0.44350182805060839D0
        ll%lvec(k)%d_mat(2,3)=0.28409190914049431D0
        ll%lvec(k)%d_mat(2,4)=0.44350182805060839D0
        ll%lvec(k)%d_mat(3,1)=0.18197943668877323D0
        ll%lvec(k)%d_mat(3,2)=0.28409190914049431D0
        ll%lvec(k)%d_mat(3,3)=0.18197943668877323D0
        ll%lvec(k)%d_mat(3,4)=0.28409190914049431D0
        ll%lvec(k)%d_mat(4,1)=0.28409190914049431D0
        ll%lvec(k)%d_mat(4,2)=0.44350182805060839D0
        ll%lvec(k)%d_mat(4,3)=0.28409190914049431D0
        ll%lvec(k)%d_mat(4,4)=0.44350182805060839D0
      else
        mattxt=adjustl(ll%debugdensfile)
        mattxt=trim(mattxt)
        write(*,*) mattxt
        CALL lsOPEN(IUNIT,mattxt,'old','UNFORMATTED')
        !OPEN(UNIT=iunit,FILE=trim(mattxt),STATUS='OLD',IOSTAT=ierror)
        read(iunit) nbasterik
        if(nbasterik .ne. nbast) then
          write(*,*) 'Not the right dimensions for the density matrix'
          write(*,*) 'are you sure you have the same basis or molecule?'
          write(lupri,*) 'Not the right dimensions for the density matrix'
          write(lupri,*) 'are you sure you have the same basis or molecule?'
          call LSquit('Not correct dimension in density matrix',lupri)
        endif
        DO j=1,nbasterik
         read(iunit) (ll%lvec(k)%d_mat(i,j),i=1,nbasterik)
        ENDDO
        CALL lsCLOSE(IUNIT,'KEEP')
      endif
      do n1=1,num_latvectors
         call mat_init(nfdensity(n1),nbast,nbast)
         call mat_zero(nfdensity(n1))
      enddo

      call find_latt_index(n1,0,0,0,fdim,ll,ll%max_layer)
      Call mat_set_from_full(ll%lvec(n1)%d_mat,1.D0,nfdensity(n1))

      write(lupri,*) 'density used'
      do j=1,nbast
         write(lupri,*) (ll%lvec(n1)%d_mat(i,j),i=1,nbast)
      enddo
      write(*,*) 'density used'
      !call write_matrix(ll%lvec(n1)%d_mat,nbast,nbast)
    else !TESTCASE

  
      do n1=1,num_latvectors
         call mat_init(nfdensity(n1),nbast,nbast)
         call mat_zero(nfdensity(n1))
      enddo
      call find_latt_index(n1,0,0,0,fdim,ll,ll%max_layer)
      call mat_copy(1.0_realk,Dmat,nfdensity(n1)) 
    endif!END ELSEIF TENSTCASE


  CALL LSTIMER('START ',TS,TE,LUPRI)
  call pbc_overlap_k(lupri,luerr,setting,input%molecule,nbast,&
      ll,latt_cell,refcell,num_latvectors,ovl)
  CALL LSTIMER('pbc_overlap_k',TS,TE,LUPRI)
 
  !CALCULATES kinetic energy of electrons
  CALL LSTIMER('START ',TS,TE,LUPRI)
  call pbc_kinetic_k(lupri,luerr,setting,input%molecule,nbast,&
   ll,latt_cell,refcell,num_latvectors,nfdensity,f_1,E_kin)
  CALL LSTIMER('pbc_kinetic_k',TS,TE,LUPRI)

  !CALCULATES electron nuclei attraction
  CALL LSTIMER('START ',TS,TE,LUPRI)
  call pbc_nucattrc_k(lupri,luerr,setting,input%molecule,nbast,&
     ll,latt_cell,refcell,num_latvectors,nfdensity,f_1,E_en)
  CALL LSTIMER('pbc_nucattrc_k',TS,TE,LUPRI)

  !CALCULATES nuclear repulsion
  CALL LSTIMER('START ',TS,TE,LUPRI)
  CALL pbc_nucpot(lupri,luerr,setting,input%molecule,ll,&
                  latt_cell,refcell,num_latvectors,E_nn)
  CALL LSTIMER('pbc_nucpot',TS,TE,LUPRI)

  
  maxmultmom=ll%lmax
  Tlmax=ll%Tlmax
  

  !CALL LSTIMER('START ',TS,TE,LUPRI)
  !call find_cutoff_twop(lupri,luerr,setting,nbast,ll,&
  !            latt_cell, refcell,num_latvectors,nfdensity)
  !CALL LSTIMER('pbc_find_twop',TS,TE,LUPRI)
 

  CALL LSTIMER('START ',TS,TE,LUPRI)
  call pbc_electron_rep_k(lupri,luerr,setting,input%molecule,nbast,&
     ll,latt_cell,refcell,num_latvectors,nfdensity,g_2,E_J)
  CALL LSTIMER('pbc Coul',TS,TE,LUPRI)

  CALL LSTIMER('START ',TS,TE,LUPRI)
  call pbc_multipole_expan_k(lupri,luerr,setting,nbast,ll,&
    &latt_cell,refcell,num_latvectors,maxmultmom)
  CALL LSTIMER('pbc_multipole',TS,TE,LUPRI)
  
  call mem_alloc(Tlat,(Tlmax+1)**2,(Tlmax+1)**2)

  call pbc_controlmm(20,Tlat,Tlmax,maxmultmom,.false.,ll%ldef%avec,&
     nbast,lupri,nfdensity,num_latvectors,nfsze,ll,g_2,E_ff,E_nn,refcell)

  CALL LSTIMER('START ',TS,TE,LUPRI)
  call pbc_exact_xc_k(lupri,luerr,setting,input%molecule,nbast,&
     ll,latt_cell,refcell,num_latvectors,nfdensity,g_2,E_K)
  CALL LSTIMER('pbc Xchange',TS,TE,LUPRI)

  write(lupri,*) 'nlayers exch',ll%kx1,ll%kx2,ll%kx3
  write(*,*) 'nlayers exch',ll%kx1,ll%kx2,ll%kx3

  ll%fc1=max(ll%oneop1,ll%col1)
  ll%fc1=max(ll%fc1,ll%Kx1)
  ll%fc2=max(ll%oneop2,ll%col2)
  ll%fc2=max(ll%fc2,ll%Kx2)
  ll%fc3=max(ll%oneop3,ll%col3)
  ll%fc3=max(ll%fc3,ll%Kx3)


!  DO n1=1,Bz%nk
!  call pbc_zdevectorize_mat(kdep(n1)%kfockmat,nbast,nbast,kdep(n1)%kfockvec)
!  call pbc_zdevectorize_mat(kdep(n1)%koverlapmat,nbast,nbast,kdep(n1)%koverlapvec)
!  ENDDO

  do i=1,num_latvectors
     call mat_free(nfdensity(i))
  enddo
  call mem_dealloc(nfdensity)

    focknorm=0.0d0
    write(lupri,*) num_latvectors
     do i=1,num_latvectors
      do j=1,nbast*nbast
         focknorm=focknorm + ll%lvec(i)%fck_vec(j)**2
      enddo
     enddo
    !CALL lsOPEN(IUNIT,'pbch2_t.dat','UNKNOWN','FORMATTED')
    E_cell=E_kin+E_en+E_J+E_K+E_ff+E_nn
    write(lupri,'(A,I4)') 'numbers of lattice vectors', num_latvectors
    write(lupri,'(A,I8)') 'number of basis', nbast
    write(lupri,'(A,F16.6)')  'Norm of fock matrix', focknorm
    write(lupri,'(A,F16.6)')  'Cell energy', E_cell
    write(lupri,'(A,F16.6)')  'Electronic energy', E_cell-E_nn-E_ff-E_K
    write(lupri,'(A,F16.6)')  'Nuclear repulsion energy', E_nn
    write(lupri,'(A,F16.6)')  'N. F. Coulomb energy', E_J
    write(lupri,'(A,F16.6)')  'Exact xch energy', E_K
    write(lupri,'(A,F16.6)')  'Far field', E_ff
    write(*,'(A,F16.6)')  'N. F. Coulomb energy', E_J
    write(*,'(A,F16.6)')  'Exact xch energy', E_K
    write(*,'(A,F16.6)')  'Far field', E_ff
    write(*,'(A,F16.6)')  'Cell energy', E_cell
    write(*,'(A,F16.6,X,F16.6,X,F16.6)')  'one part energy', E_kin+E_en,E_kin,E_en
    write(*,'(A,F16.6)')  'Electronic energy', E_cell-E_nn-E_ff-E_K
    write(*,'(A,F16.6)')  'Nuclear repulsion energy', E_nn
    !CALL lsCLOSE(IUNIT,'KEEP')
    

  if(ll%compare_elmnts) then

    allocate(k_fock(nbast,nbast))
    allocate(k_Sab(nbast,nbast))

    do n1=-3,3
    write(numtostring1,*) n1
    numtostring1=adjustl(numtostring1)
    mattxt='minFmat1'//trim(numtostring1)//'00.dat'
      !call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
    CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
    call find_latt_index(k,n1,0,0,fdim,ll,ll%max_layer)
    write(iunit,*) k
    DO j=1,nbast
       write(iunit,*) (ll%lvec(k)%fck_vec(i+(j-1)*nbast),i=1,nbast)
    ENDDO
    call lsclose(iunit,'KEEP')
    enddo !n1
    write(*,*) 'Fock matrix written to disk'
    write(lupri,*) 'Fock matrix written to disk'

    if(ll%store_mats) then
      call pbc_read_fock_matrix(ll,nbast,nbast,'            ')
      call pbc_fockmat_write(ll,nbast,nbast,7,2,'            ',lupri)
    endif

    !call COMPARE_MATRICES(lupri,nbast,num_latvectors,nfsze,maxmultmom,ll)
    call readerikmats(input%molecule,setting,k_fock,k_Sab,nbast,ll,&
    num_latvectors,latt_cell,nfsze,maxmultmom,bz,tlat,lupri,luerr)

    deallocate(k_fock)
    deallocate(k_Sab)
    call mem_dealloc(Tlat)

  else

    call mem_dealloc(Tlat)
    !call init_pbc_elstr(bz%fck,nbast,nbast)
    !call init_pbc_elstr(bz%smat,nbast,nbast)
	do i=1,num_latvectors
	 call free_Moleculeinfo(latt_cell(i))
	enddo
    deallocate(latt_cell)

    E_1=E_kin+E_en+E_nn
    if(ll%num_its .gt. 0) call pbc_startzdiis(input%molecule,setting,nbast,ll,&
    num_latvectors,nfsze,maxmultmom,bz,ovl,f_1,g_2,E_1,lupri,luerr)

  endif

	do i=1,num_latvectors
     if(f_1(i)%init_magic_tag.EQ.mat_init_magic_value) then
       call mat_free(f_1(i))
     endif
     if(ovl(i)%init_magic_tag.EQ.mat_init_magic_value) then
       call mat_free(ovl(i))
     endif
	enddo
    call mem_dealloc(f_1)
    call mem_dealloc(g_2)
    call mem_dealloc(ovl)
	call free_Moleculeinfo(refcell)







!!!Wannier orbitals, may most probably delete it.

  CASE('directly')

	sze=nbast*size(ll%lvec)
!   ALLOCATE(S_ab(sze,sze))
!   S_ab=0.0000
!    write(lupri,*) 'Overlap'
  
    write(lupri,*) 'We are now in the ',ll%wannier_direct, ' method'

    allocate(fck(sze*sze))
    fck=0.0000E-100
  
    call pbc_overlap_int(lupri,luerr,setting,input%molecule,&
    nbast,ll,latt_cell,refcell,num_latvectors)
  
    call pbc_kinetic_int(lupri,luerr,setting,input%molecule,nbast,fck,&
     sze,ll,latt_cell,refcell,num_latvectors)
  
    call pbc_nucattrc_int(lupri,luerr,setting,input%molecule,nbast,fck,&
     sze,ll,latt_cell,refcell,num_latvectors)
     
    call pbc_get_nfsize(n1,n2,n3,ll%nneighbour,lupri)
    nfsze=(2*n1+1)*(2*n2+1)*(2*n3+1)
    write(lupri,*) 'nfsize: ',nfsze
    call mem_alloc(nfdensity,nfsze)

    k=0
    l1=0
    l2=0
    l3=0
    !do l1=-3,3
    iunit = -1
    scfit=1
    write(string1,'(I5)') scfit
    string1=adjustl(string1)
    k=1
    write(numtostring1,'(I5)')  l1
    write(numtostring2,'(I5)')  l2
    write(numtostring3,'(I5)')  l3
    numtostring1=adjustl(numtostring1)
    numtostring2=adjustl(numtostring2)
    numtostring3=adjustl(numtostring3)
    
    write(mattxt,'(A20)') 'PBCDMAT'//trim(string1)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)

    mattxt=adjustl(mattxt)
    mattxt=trim(mattxt)
    write(*,*) mattxt

    CALL lsOPEN(IUNIT,mattxt,'old','UNFORMATTED')
    read(iunit) nbasterik
    if(nbasterik .eq. nbast) THEN
      DO j=1,nbasterik
       read(iunit) (ll%lvec(k)%d_mat(i,j),i=1,nbasterik)
      ENDDO
      do n1=1,nfsze
         call mat_init(nfdensity(n1),nbast,nbast)
         call mat_zero(nfdensity(n1))
      enddo
      call find_latt_index(n1,0,0,0,fdim,ll,3)
      Call mat_set_from_full(ll%lvec(k)%d_mat,1.D0,nfdensity(n1))
      CALL lsCLOSE(IUNIT,'KEEP')
    else
    !mat_to_full(a, alpha, afull)
     CALL lsCLOSE(IUNIT,'KEEP')
    !enddo



    call mat_init(D1(1),nbast,nbast)
    call mat_init(D(1),nbast,nbast)
    call mat_zero(D(1))
    Dfull=0.0D0
    call make_uptriag1_mat(Dfull,nbast,nbast)
  
    Call mat_set_from_full(Dfull,1.D0,D(1))
    Call mat_set_from_full(Dfull,1.D0,D1(1))
  
    do n1=1,nfsze
       call mat_init(nfdensity(n1),nbast,nbast)
       call mat_zero(nfdensity(n1))
       Call mat_set_from_full(Dfull,1.D0,nfdensity(n1))
    enddo
  endif
 
    call pbc_electron_rep(lupri,luerr,setting,input%molecule,nbast,&
     ll,latt_cell,refcell,num_latvectors,nfdensity,nfsze)
  
!  call pbc_complete_Fock_mtx(lupri,nbast,S_ab,sze,cutoff,ll)

    do n1=1,nfsze
       call mat_free(nfdensity(n1))
    enddo

    call mem_dealloc(nfdensity)

    call mat_free(D1(1))
    call mat_free(D(1))


  END SELECT


  call free_lvec_list(ll)
  write(lupri,*) 'Program ended successfully !'
  write(*,*) 'Program ended successfully !'

END SUBROUTINE set_pbc_molecules

SUBROUTINE pbc_init_recvec(realspace,recvec)
IMPLICIT NONE
REAL(realk),INTENT(IN) :: realspace(3,3)
REAL(realk),INTENT(INOUT) :: recvec(3,3)
!REAL(realk) :: PI=3.14159265358979323846D0

    recvec(1,1)=2.*pi*(realspace(2,2)*realspace(3,3)-&
              realspace(3,2)*realspace(2,3))
    recvec(2,1)=2.*pi*(realspace(3,2)*realspace(1,3)-&
              realspace(1,2)*realspace(3,3))
    recvec(3,1)=2.*pi*(realspace(1,2)*realspace(2,3)-&
              realspace(2,2)*realspace(1,3))

    recvec(1,2)=2.*pi*(realspace(2,3)*realspace(3,1)-&
              realspace(3,3)*realspace(2,1))
    recvec(2,2)=2.*pi*(realspace(3,3)*realspace(1,1)-&
              realspace(1,3)*realspace(3,1))
    recvec(3,2)=2.*pi*(realspace(1,3)*realspace(2,1)-&
              realspace(2,3)*realspace(1,1))

    recvec(1,3)=2.*pi*(realspace(2,1)*realspace(3,3)-&
              realspace(3,1)*realspace(2,3))
    recvec(2,3)=2.*pi*(realspace(3,1)*realspace(1,2)-&
              realspace(1,1)*realspace(3,2))
    recvec(3,3)=2.*pi*(realspace(1,1)*realspace(2,2)-&
    	      &realspace(2,1)*realspace(1,2))

    recvec=2.*pi*recvec/(recvec(1,1)*realspace(1,1)+recvec(2,1)*realspace(2,1)+&
         recvec(3,1)*realspace(3,1))!This gives NaN for the moment

END SUBROUTINE pbc_init_recvec

!SUBROUTINE pbc_scf(SETTING,lupri,luerr,nbast)
!  USE TYPEDEF
!  USE matrix_module
!  use lattice_vectors
!  IMPLICIT NONE
!  !INPUT VARIABLES
!  TYPE(LSSETTING):: SETTING
!  INTEGER, INTENT(IN) :: lupri, luerr, nbast
!  !LOCAL VARIABLES
!  TYPE(MOLECULEINFO) :: refcell
!  TYPE(MATRIX)  :: S
!  REAL(realk), ALLOCATABLE ::  S_ab(:,:)
!  LOGICAL :: threshold
!
!  threshold =.true.
!  !Do while(threshold)
!  !call make_S_matrix
!  !call make_kin_matrix !set_pbc_molecules !then we get fock matrix
!  !call make_nucl_matrix
!  !call make_coulomb
!  !call make_fock
!  !call complete fock_matrix
!  !call solvfock
!  !call make_density
!  !call find_energy
!  !ENDDO
!
!END SUBROUTINE pbc_scf
END MODULE pbc_setup
