!GIVE MODULE NAME
MODULE pbc_setup
use files
!#ifdef MOD_UNRELEASED
  USE precision
  USE fundamental
  USE TYPEDEF
  USE matrix_module
  USE lattice_vectors
  USE lattice_type
  use pbc_msc
  USE memory_handling
!  USE multipole_pbc
!  USE harmonics_pbc
  USE pbc_matrix_operations
  USE pbc_scfdiis
  USE pbc_interactions
  USE pbc_ff_contrib
  private
  public :: set_pbc_molecules
contains

!Subroutine to copy atoms and interface it with the integration code developed
!by Simen Reine.
SUBROUTINE set_pbc_molecules(INPUT,SETTING,lupri,luerr,nbast,Dmat,lattice)
!  use pbc_matrix_module
! PROPERTIES SECTION
!  use lsdalton_rsp_mod
  implicit none
  !INPUT VARIABLES
  INTEGER, INTENT(IN) :: lupri, luerr, nbast
  TYPE(lvec_list_t), INTENT(INOUT) :: lattice
  TYPE(DALTONINPUT),INTENT(INOUT) :: INPUT
  TYPE(LSSETTING):: SETTING
  TYPE(MATRIX),intent(INOUT)  ::  Dmat
  !Local variables
  TYPE(MOLECULEINFO) :: refcell
  TYPE(MATRIX)  ::  D(1),D1(1)!,F(1)
  TYPE(matrix),pointer :: nfdensity(:),f_1(:),ovl(:)
  TYPE(matrix),pointer :: g_2(:)
  !TYPE(lattice_cell_info_t),pointer :: sphermom(:)
  REAL(realk) :: Dfull(nbast,nbast),kvec(3,3),origin(3)
  complex(complexk), pointer :: fck(:), S_ab(:),k_fock(:,:),K_Sab(:,:)
  real(realk), pointer :: Tlat(:,:)
  REAL(realk)::  latt_vec_std(3),focknorm
  REAL(realk)::  E_cell,E_kin,E_ff,E_XC,E_K,E_J,E_en,E_nuc
  REAL(realk)::  E_1,E_nnff,maxdens
  real(realk) :: TS,TE
!  REAL(realk) :: PI=3.14159265358979323846D0
  INTEGER :: sze,num_latvectors
  INTEGER :: maxmultmom,n1,n2,n3,nfsze,Tlmax
  INTEGER :: i,j,k,scfit,iunit,l1,l2,l3,nbasterik,ierror
  character*(20) :: mattxt,string1,numtostring1,numtostring2,numtostring3
  character(len=3) ::nline
  TYPE(BZgrid_t) :: BZ


!!!!!!!!!!!!Information for me
!!!!! input%Basis%binfo(regbasparam)%atomtype(:)%shell(:)%segment(1)%exponents(:)
write(lupri,*) 'shells ',input%Basis%binfo(regbasparam)%atomtype(1)%nangmom
write(lupri,*) 'Number of atom species ',input%Basis%binfo(regbasparam)%natomtypes
write(lupri,*) 'Exponents ',(input%Basis%binfo(regbasparam)%atomtype(1)%shell(1)%segment(1)%exponents)
!  call screening_ovl(input%Basis%binfo(regbasparam))
  call build_lvec_list(lattice,nbast) 

  write(*,*) 'nearest neighbour in pbcmain', lattice%nneighbour
  pbc_control%ldef%is_active=lattice%ldef%is_active
  nearfield=lattice%nf
  write(*,*) '1 active not active:',lattice%ldef%is_active(1)
  write(*,*) '2 active not active:',lattice%ldef%is_active(2)
  write(*,*) '3 active not active:',lattice%ldef%is_active(3)

  call set_refcell(refcell,input%molecule)
  
 ! write(lupri,*) lattice%ldef%avec
  origin=1.0d0
 !origin of cell should be 0d0, that means that I have to change the origin for
 !the next cell to what it should be and the atomic position in
 !each cell should be according to the origin of the new cell, a copy of
 !the reference cell.

  
  write(*,*) 'number of layers', lattice%max_layer
  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)
  
  num_latvectors=size(lattice%lvec)
  write(lupri,*) 'Number of vectors ', num_latvectors


  call set_lattice_cells(num_latvectors,input%molecule,lattice,lupri)


  write(lupri,*) 'setting%p%atom%center', lattice%lvec(2)%molecule%atom(1)%center(1)

  call reset_integral_computed(lattice,num_latvectors,'all')

  write(lupri,*) 'cutoff' 
  call find_cutoff_onep(lupri,luerr,setting,nbast,refcell,lattice)

!  write(*,*) 'nearest neighbour in pbcmain 2 ', lattice%nneighbour
!  call build_nflvec_list(lattice,nbast) 
!  write(*,*) 'nearest neighbour in pbcmain 3 ', lattice%nneighbour
!  n_neighbour=lattice%nneighbour

 
  SELECT CASE(lattice%wannier_direct)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BLOCH FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CASE('indirectly')

    write(lupri,*) 'We are now in the ',lattice%wannier_direct, ' method'
    write(lupri,*) 'CS screen:', setting%scheme%CS_SCREEN 
    write(lupri,*) 'PS screen:', setting%scheme%PS_SCREEN 
    write(*,*) 'CS screen:', setting%scheme%CS_SCREEN 
    write(*,*) 'PS screen:', setting%scheme%PS_SCREEN 

   ! allocate(S_ab(nbast*nbast))
   ! S_ab=0.0000
  
   ! allocate(fck(nbast*nbast))
   ! fck=0.0000E-100
  
!
    call pbc_init_recvec(lattice%ldef%avec,kvec,lattice%ldef%is_active,lupri)


    write(lupri,*) 'kvec', kvec!(3,3)
    write(lupri,*) 'a1 dot b1 = ',dot_product(lattice%ldef%avec(1,:),kvec(1,:))
    write(lupri,*) 'lat vec', lattice%ldef%avec!(3,3)
    write(lupri,*) 'module test: ', mod(3,3)

    call pbc_init_Bzgrid(kvec,lattice%nk1,lattice%nk2,lattice%nk3,'nosym',bZ,nbast,nbast,lupri)
    lat_data%num_k1=lattice%nk1
    lat_data%num_k2=lattice%nk2
    lat_data%num_k3=lattice%nk3
    lat_data%num_kpoints=bz%nk 
    lat_data%reclatvec(:,:)=kvec(:,:)
  !  write(*,*) 'Number kpoints', Bz%Nk
  !  write(*,*) 'value of kpoint', BZ%kpnt(1)%n(1:3)
  !  call pbc_get_kpoint(1,kvec(1,:))
  !  write(*,*) 'kpoint realvalue',kvec(1,:)
  !  STOP

  kvec=0.0d0

  write(lupri,*) 'kvec(3,3)', kvec(3,3)

  write(*,*) 'lattice%nneighbour', lattice%nneighbour
  write(*,*) 'lattice%nf', lattice%nf
  write(*,*) 'Number of k points', Bz%nk !TODO Should be nk_nosym ????
  write(lupri,*) 'Number of k points', Bz%nk

    k=0
    l1=0
    l2=0
    l3=0

#ifdef DEBUGPBC

  if(lattice%compare_elmnts) then
  mattxt=adjustl('PBCDMAT000')
  mattxt=trim(mattxt)
  write(*,*) mattxt
  call find_latt_index(n1,0,0,0,lattice,lattice%max_layer)
  CALL lsOPEN(IUNIT,mattxt,'old','FORMATTED')
  !OPEN(UNIT=iunit,FILE=trim(mattxt),STATUS='OLD',IOSTAT=ierror)
  read(iunit,*) nbasterik
  if(nbasterik .ne. nbast) then
    write(*,*) 'Not the right dimensions for the density matrix'
    write(*,*) 'are you sure you have the same basis or molecule?'
    write(lupri,*) 'Not the right dimensions for the density matrix'
    write(lupri,*) 'are you sure you have the same basis or molecule?'
    call LSquit('Not correct dimension in density matrix',lupri)
  endif
  !k=0
  !DO i=1,nbast
   DO j=1,nbast
   !k=k+1
   !nline='no'
   !if(k .eq. nbast) THEN
   !  k=0
   !  nline='yes'
   !ENDIF
   !read(iunit,*,advance=nline) lattice%lvec(k)%d_mat(i,j)
   read(iunit,*) (lattice%lvec(n1)%d_mat(i,j),i=1,nbasterik)
  ! ENDDO
  ENDDO
  write(*,*) 'dmat'
  write(*,*) lattice%lvec(n1)%d_mat
  CALL lsCLOSE(IUNIT,'KEEP')

  !do n1=1,num_latvectors
  !write(*,*) 'debugpbc num_latvectors',num_latvectors,n1
  !   call mat_init(nfdensity(n1),nbast,nbast)
  !write(*,*) 'debugpbc 4'
  !   call mat_zero(nfdensity(n1))
  !write(*,*) 'debugpbc 4'
  !enddo

  !call find_latt_index(n1,0,0,0,lattice,lattice%max_layer)
  Call mat_set_from_full(lattice%lvec(n1)%d_mat,1.D0,dmat)
  write(*,*) 'density'
  write(lupri,*) 'density'
  call mat_print(dmat,1,nbast,1,nbast,lupri)
  write(*,*) 'density written to LSDALTON.OUT'


  elseif(lattice%read_file) then

        call find_latt_index(n1,0,0,0,lattice,lattice%max_layer)
        mattxt=adjustl(lattice%debugdensfile)
        mattxt=trim(mattxt)
        write(*,*) mattxt
        CALL lsOPEN(IUNIT,mattxt,'old','FORMATTED')
        !OPEN(UNIT=iunit,FILE=trim(mattxt),STATUS='OLD',IOSTAT=ierror)
        read(iunit,*) nbasterik
        if(nbasterik .ne. nbast) then
          write(*,*) 'Not the right dimensions for the density matrix'
          write(*,*) 'are you sure you have the same basis or molecule?'
          write(*,*) 'your dim = ',nbast
          write(*,*) 'Read dim = ',nbasterik
          write(lupri,*) 'Not the right dimensions for the density matrix'
          write(lupri,*) 'are you sure you have the same basis or molecule?'
          write(lupri,*) 'your dim = ',nbast
          write(lupri,*) 'Read dim = ',nbasterik
          call LSquit('Not correct dimension in density matrix',lupri)
        endif
        DO j=1,nbasterik
         read(iunit,*) (lattice%lvec(n1)%d_mat(i,j),i=1,nbasterik)
        ENDDO
        CALL lsCLOSE(IUNIT,'KEEP')
  
        Call mat_set_from_full(lattice%lvec(n1)%d_mat,1.D0,dmat)

    elseif(lattice%testcase) THEN !THIS IS FOR DEBUGGING
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


      call find_latt_index(k,0,0,0,lattice,lattice%max_layer)
      !do i=1,num_latvectors
      !  call init_lvec_data(lattice%lvec(i),nbast)
      !enddo
      call mem_alloc(lattice%lvec(k)%d_mat,nbast,nbast)
      if(.not.lattice%read_file) then
        write(lupri,*) 'READING density mat'
        lattice%lvec(k)%d_mat(1,1)=0.18197943668877323D0
        lattice%lvec(k)%d_mat(1,2)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(1,3)=0.18197943668877323D0
        lattice%lvec(k)%d_mat(1,4)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(2,1)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(2,2)=0.44350182805060839D0
        lattice%lvec(k)%d_mat(2,3)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(2,4)=0.44350182805060839D0
        lattice%lvec(k)%d_mat(3,1)=0.18197943668877323D0
        lattice%lvec(k)%d_mat(3,2)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(3,3)=0.18197943668877323D0
        lattice%lvec(k)%d_mat(3,4)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(4,1)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(4,2)=0.44350182805060839D0
        lattice%lvec(k)%d_mat(4,3)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(4,4)=0.44350182805060839D0
        Call mat_set_from_full(lattice%lvec(k)%d_mat,1.D0,dmat)
        call mem_dealloc(lattice%lvec(k)%d_mat)
      else
        mattxt=adjustl(lattice%debugdensfile)
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
         read(iunit) (lattice%lvec(k)%d_mat(i,j),i=1,nbasterik)
        ENDDO
        CALL lsCLOSE(IUNIT,'KEEP')

      call find_latt_index(n1,0,0,0,lattice,lattice%max_layer)
      Call mat_set_from_full(lattice%lvec(n1)%d_mat,1.D0,dmat)
      !call mat_copy(1.0_realk,Dmat,nfdensity(n1)) 
      call mem_dealloc(lattice%lvec(k)%d_mat)
    endif
        

else


      !call find_latt_index(n1,0,0,0,lattice,lattice%max_layer)
      !call mat_copy(1.0_realk,Dmat,nfdensity(n1)) 
    endif



#else

    write(*,*) 'before lattice%testcase'
    if(lattice%testcase) THEN !THIS IS FOR DEBUGGING
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


      call find_latt_index(k,0,0,0,lattice,lattice%max_layer)
      !do i=1,num_latvectors
      !  call init_lvec_data(lattice%lvec(i),nbast)
      !enddo
      call mem_alloc(lattice%lvec(k)%d_mat,nbast,nbast)
      if(.not.lattice%read_file) then
        lattice%lvec(k)%d_mat(1,1)=0.18197943668877323D0
        lattice%lvec(k)%d_mat(1,2)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(1,3)=0.18197943668877323D0
        lattice%lvec(k)%d_mat(1,4)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(2,1)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(2,2)=0.44350182805060839D0
        lattice%lvec(k)%d_mat(2,3)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(2,4)=0.44350182805060839D0
        lattice%lvec(k)%d_mat(3,1)=0.18197943668877323D0
        lattice%lvec(k)%d_mat(3,2)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(3,3)=0.18197943668877323D0
        lattice%lvec(k)%d_mat(3,4)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(4,1)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(4,2)=0.44350182805060839D0
        lattice%lvec(k)%d_mat(4,3)=0.28409190914049431D0
        lattice%lvec(k)%d_mat(4,4)=0.44350182805060839D0
      else
        mattxt=adjustl(lattice%debugdensfile)
        mattxt=trim(mattxt)
        write(*,*) mattxt
        CALL lsOPEN(IUNIT,mattxt,'old','FORMATTED')
        !OPEN(UNIT=iunit,FILE=trim(mattxt),STATUS='OLD',IOSTAT=ierror)
        read(iunit,*) nbasterik
        if(nbasterik .ne. nbast) then
          write(*,*) 'Not the right dimensions for the density matrix'
          write(*,*) 'are you sure you have the same basis or molecule?'
          write(*,*) 'Your dimension is',nbast,'read dimension is',nbasterik
          write(lupri,*) 'Not the right dimensions for the density matrix'
          write(lupri,*) 'are you sure you have the same basis or molecule?'
          call LSquit('Not correct dimension in density matrix',lupri)
        endif
        DO j=1,nbasterik
         read(iunit,*) (lattice%lvec(k)%d_mat(i,j),i=1,nbasterik)
        ENDDO
        CALL lsCLOSE(IUNIT,'KEEP')
      endif
      !do n1=1,num_latvectors
      !   call mat_init(nfdensity(n1),nbast,nbast)
      !   call mat_zero(nfdensity(n1))
      !enddo

      call find_latt_index(n1,0,0,0,lattice,lattice%max_layer)
      Call mat_set_from_full(lattice%lvec(n1)%d_mat,1.D0,Dmat)

      write(lupri,*) 'density used'
      do j=1,nbast
         write(lupri,*) (lattice%lvec(n1)%d_mat(i,j),i=1,nbast)
      enddo
      write(*,*) 'density used'
      !call write_matrix(lattice%lvec(n1)%d_mat,nbast,nbast)
      call mem_dealloc(lattice%lvec(n1)%d_mat)

    else !TESTCASE

  
      !do n1=1,num_latvectors
      !   call mat_init(nfdensity(n1),nbast,nbast)
      !   call mat_zero(nfdensity(n1))
      !enddo
      !call find_latt_index(n1,0,0,0,lattice,lattice%max_layer)
      !call mat_copy(1.0_realk,Dmat,nfdensity(n1)) 


    endif!END ELSEIF TESTCASE DEBUGGING
#endif

#ifdef DEBUGPBC
      !write(lupri,*) 'density used'
      !call mat_to_full(nfdensity(n1), 1.0_realk,lattice%lvec(n1)%d_mat)
      !call write_matrix(lattice%lvec(n1)%d_mat,nbast,nbast,lupri)
      !write(*,*) 'density used'
      write(lupri,*) 'Density first'
      call mat_print(dmat,1,nbast,1,nbast,lupri)
#endif

  !Preparation for the far field contribution, setup multipole
  !moments and the Tlattice
  maxmultmom=lattice%lmax
  Tlmax=lattice%Tlmax

 ! call mem_alloc(Tlat,(Tlmax+1)**2,(Tlmax+1)**2)

 ! write(*,*) 'density used ',num_latvectors
 ! call pbc_controlmm(20,Tlat,Tlmax,maxmultmom,.false.,lattice%ldef%avec,&
 !    nbast,lupri,dmat,num_latvectors,lattice,E_ff,E_nnff,refcell)

#ifdef DEBUGPBC
  if(lattice%compare_elmnts) then
    !write(*,*) 'hei'
    call readerikmats(input%molecule,setting,k_fock,k_Sab,nbast,lattice,&
    & num_latvectors,nfsze,maxmultmom,bz,tlat,lupri,luerr)
  endif

#endif
    write(lupri,'(A,I4)') 'numbers of lattice vectors', num_latvectors
    write(lupri,'(A,I8)') 'number of basis', nbast
    
  if(lattice%compare_elmnts) then

    allocate(k_fock(nbast,nbast))
    allocate(k_Sab(nbast,nbast))

    do n1=-3,3
    write(numtostring1,*) n1
    numtostring1=adjustl(numtostring1)
    mattxt='minFmat1'//trim(numtostring1)//'00.dat'

    !call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
    !CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
    !call find_latt_index(k,n1,0,0,lattice,lattice%max_layer)
    !write(iunit,*) k
    !DO j=1,nbast
    !   write(iunit,*) (lattice%lvec(k)%fck_vec(i+(j-1)*nbast),i=1,nbast)
    !ENDDO
    !call lsclose(iunit,'KEEP')
    enddo !n1
    write(*,*) 'Fock matrix written to disk'
    write(lupri,*) 'Fock matrix written to disk'

    if(lattice%store_mats) then
      call pbc_read_fock_matrix(lattice,nbast,nbast,'            ')
      call pbc_fockmat_write(lattice,nbast,nbast,7,2,'            ',lupri)
    endif

    !call COMPARE_MATRICES(lupri,nbast,num_latvectors,nfsze,maxmultmom,lattice)
    call readerikmats(input%molecule,setting,k_fock,k_Sab,nbast,lattice,&
    & num_latvectors,nfsze,maxmultmom,bz,tlat,lupri,luerr)

    deallocate(k_fock)
    deallocate(k_Sab)
    call mem_dealloc(Tlat)

  else
!    call mem_dealloc(Tlat)
    !call init_pbc_elstr(bz%fck,nbast,nbast)
    !call init_pbc_elstr(bz%smat,nbast,nbast)
	do i=1,num_latvectors
	 call free_Moleculeinfo(lattice%lvec(i)%molecule)
	enddo
    write(*,*) 'before call to scf loops'

   call mat_abs_max_elm(dmat,maxdens)
   write(lupri,*) 'max element in initial dmat',maxdens
   write(*,*) 'max element in initial dmat',maxdens

  call pbc_startzdiis(input%molecule,setting,nbast,lattice,&
  & num_latvectors,maxmultmom,bz,dmat,lupri,luerr)

  call pbc_end_Bzgrid(bZ)
!  call mat_free(nfdensity(n1))
!  call mem_dealloc(nfdensity)

 !   E_1=E_kin+E_en+E_nuc
 !   if(lattice%num_its .gt. 0) call pbc_startzdiis(input%molecule,setting,nbast,lattice,&
 !   num_latvectors,nfsze,maxmultmom,bz,ovl,f_1,g_2,E_nuc,lupri,luerr)

  endif

!	do i=1,num_latvectors
    ! if(f_1(i)%init_magic_tag.EQ.mat_init_magic_value) then
    !   call mat_free(f_1(i))
    ! endif
    ! if(ovl(i)%init_magic_tag.EQ.mat_init_magic_value) then
    !   call mat_free(ovl(i))
    ! endif
    !    enddo
    !call mem_dealloc(f_1)
    !call mem_dealloc(g_2)
    !call mem_dealloc(ovl)
        call free_Moleculeinfo(refcell)







!!!Wannier orbitals, may most probably delete it.

  CASE('directly')

	sze=nbast*size(lattice%lvec)
!   ALLOCATE(S_ab(sze,sze))
!   S_ab=0.0000
!    write(lupri,*) 'Overlap'
  
    write(lupri,*) 'We are now in the ',lattice%wannier_direct, ' method'

    allocate(fck(sze*sze))
    fck=0.0000E-100_realk
  
    call pbc_overlap_int(lupri,luerr,setting,nbast,lattice, &
		 & refcell,num_latvectors)
  
    call pbc_kinetic_int(lupri,luerr,setting,input%molecule%natoms,nbast,fck,&
     & sze,lattice,refcell,num_latvectors)
  
    call pbc_nucattrc_int(lupri,luerr,setting,input%molecule%natoms,nbast,fck,&
     & sze,lattice,refcell,num_latvectors)
     
    call pbc_get_nfsize(n1,n2,n3,lattice%nneighbour,lupri)
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
       read(iunit) (lattice%lvec(k)%d_mat(i,j),i=1,nbasterik)
      ENDDO
      do n1=1,nfsze
         call mat_init(nfdensity(n1),nbast,nbast)
         call mat_zero(nfdensity(n1))
      enddo
      call find_latt_index(n1,0,0,0,lattice,3)
      Call mat_set_from_full(lattice%lvec(k)%d_mat,1.D0,nfdensity(n1))
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
 
    call pbc_electron_rep(lupri,luerr,setting,input%molecule%natoms,nbast,&
     & lattice,refcell,num_latvectors,nfdensity,nfsze)
  
!  call pbc_complete_Fock_mtx(lupri,nbast,S_ab,sze,cutoff,ll)

    do n1=1,nfsze
       call mat_free(nfdensity(n1))
    enddo

    call mem_dealloc(nfdensity)

    call mat_free(D1(1))
    call mat_free(D(1))


  END SELECT


  !call free_lvec_list(lattice)
  call mem_dealloc(lattice%lvec)
  write(lupri,*) 'Program ended successfully !'
  write(*,*) 'Program ended successfully !'

END SUBROUTINE set_pbc_molecules


!TODO move this s.r. to pbc-msc ?? 

!> @brief Calculate basis vectors of the reciprocal space.
!> @param realspace Primitive vectors of the lattice.
!> @param recvec Reciprocal vectors.
!> @param is_active Active dimensions in the calculation.
!> @param lu 
SUBROUTINE pbc_init_recvec(realspace,recvec,is_active,lu)
	IMPLICIT NONE
	! input
	REAL(realk),INTENT(IN) :: realspace(3,3)
	REAL(realk),INTENT(INOUT) :: recvec(3,3)
	LOGICAL,INTENT(IN)     :: is_active(3)
	INTEGER,INTENT(IN)     :: lu
	! local
	REAL(realk)            :: t2ct3(3), norm_a(2), tmp(3)
	REAL(realk)            :: a(3, 2), b(3, 2), e(3, 2), rotmat(3, 3)
	REAL(realk)            :: vol
	INTEGER                :: active_dims, i, j
    
	! count the number of active dimensions
	active_dims = 0
	do i = 1, 3
      if(is_active(i)) active_dims = active_dims + 1
	enddo
   recvec(:,:) = 0.0_realk

	SELECT CASE(active_dims)

	CASE(3)

		recvec(1,1)=2.*pi*(realspace(2,2)*realspace(3,3)-&
			& realspace(3,2)*realspace(2,3))
		recvec(2,1)=2.*pi*(realspace(3,2)*realspace(1,3)-&
			& realspace(1,2)*realspace(3,3))
		recvec(3,1)=2.*pi*(realspace(1,2)*realspace(2,3)-&
			& realspace(2,2)*realspace(1,3))

		recvec(1,2)=2.*pi*(realspace(2,3)*realspace(3,1)-&
			& realspace(3,3)*realspace(2,1))
		recvec(2,2)=2.*pi*(realspace(3,3)*realspace(1,1)-&
			& realspace(1,3)*realspace(3,1))
		recvec(3,2)=2.*pi*(realspace(1,3)*realspace(2,1)-&
			& realspace(2,3)*realspace(1,1))

		recvec(1,3)=2.*pi*(realspace(2,1)*realspace(3,3)-&
			& realspace(3,1)*realspace(2,3))
		recvec(2,3)=2.*pi*(realspace(3,1)*realspace(1,2)-&
			& realspace(1,1)*realspace(3,2))
		recvec(3,3)=2.*pi*(realspace(1,1)*realspace(2,2)-&
			&realspace(2,1)*realspace(1,2))

		t2ct3(1)=realspace(2,2)*realspace(3,3)-&
			& realspace(3,2)*realspace(2,3)

		t2ct3(2)=realspace(3,2)*realspace(1,3)-&
			& realspace(1,2)*realspace(3,3)

		t2ct3(3)=realspace(1,2)*realspace(2,3)-&
			& realspace(2,2)*realspace(1,3)

		vol= dot_product(t2ct3,realspace(:,1))

		recvec(:,:)=recvec(:,:)/vol

	CASE(2)

		! the formula for the 2D reciprocal vectors is
	   ! $b_i = 2\pi \frac{ E a_j }{ a_i \cdot E a_j }$	
		! where $E = (e_i \otimes \e_j - e_j \otimes e_i)$, $e_i = a_i/|a_i|$ 
		! and $i\neq j, \, i,j\in\{1,2\}$
		
		! pick the 'active' primitive vectors
		do i = 1, 3
			if ( is_active(i) ) then
				a(:, i) = realspace(:, i)
			endif
		enddo
		
		! calculate reciprocal vectors
		do i = 1, 2
			norm_a(i) = sqrt( dot_product(a(:, i), a(:, i)) )
			e(:, i) = a(:, i) / norm_a(i)
		enddo
		do i = 1, 3
			rotmat(:, i) = e(i, 1)*e(:, 2) - e(i, 2)*e(:, 1)
		enddo
		do i = 1, 3
			tmp(i) =  dot_product(rotmat(:, i), e(:, 2)) 
		enddo
		b(:, 1) = 2.0_realk * pi * tmp(:) &
			/ ( norm_a(1) * dot_product(e(:, 1), tmp(:)) )
		do i = 1, 3
			tmp(i) =  dot_product(rotmat(:, i), e(:, 1)) 
		enddo
		b(:, 2) = 2.0_realk * pi *tmp(:) &
			/ ( norm_a(2) * dot_product(e(:, 2), tmp(:)) )

		! set reciprocal vectors
		j = 1
		do i = 1, 3
			if ( is_active(i) ) then
				recvec(:, i) = b(:, j)
				j = j + 1
			else
				recvec(:, i) = 0.0_realk
			endif
		enddo
		
!		write (*, *) 'test'
!		write (*, *) 'a1b1', dot_product(a(:,1), b(:,1)), 'should be 2 pi'
!		write (*, *) 'a1b2', dot_product(a(:,1), b(:,2)), 'should be 0'
!		write (*, *) 'a2b1', dot_product(a(:,2), b(:,1)), 'should be 0'
!		write (*, *) 'a2b2', dot_product(a(:,2), b(:,2)), 'should be 2 pi'
!		
!	 	write (*, *) recvec

	CASE(1)
		if(realspace(1,1) .ne. 0) then
			recvec(1,1)=2._realk*pi/realspace(1,1)
			if(realspace(2,1) .ne. 0._realk .or. realspace(3,1) .ne. 0._realk) then
				write(*,*) 'Use just one axis it is a one dimensional system'
				call LSquit('Not correct usage of lattice vectors',lu)
			endif
		elseif(realspace(2,1) .ne. 0) then
			recvec(2,1)=2._realk*pi/realspace(2,1)
			if(realspace(1,1) .ne. 0._realk .or. realspace(3,1) .ne. 0._realk) then
				write(*,*) 'Use just one axis it is a one dimensional system'
				call LSquit('Not correct usage of lattice vectors',lu)
			endif
		elseif(realspace(3,1) .ne. 0) then
			recvec(3,1)=2._realk*pi/realspace(3,1)
			if(realspace(2,1) .ne. 0._realk .or. realspace(1,1) .ne. 0._realk) then
				write(*,*) 'Use just one axis it is a one dimensional system'
				call LSquit('Not correct usage of lattice vectors',lu)
			endif
		endif

	END SELECT



	!recvec=2.*pi*recvec/(recvec(1,1)*realspace(1,1)+recvec(2,1)*realspace(2,1)+&
	!     recvec(3,1)*realspace(3,1))!This gives NaN for the moment

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
!#else
!contains
!  subroutine pbc_setup_empty()
!  end subroutine pbc_setup_empty
!#endif

END MODULE pbc_setup
