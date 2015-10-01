SUBROUTINE readerikmats(molecule,setting,fock,Sabk,ndim,ll,numrealvec,&
           & nfsze,lmax,bz,tlat,lupri,luerr)
  USE TYPEDEF
  USE precision
  USE matrix_module
  USE lattice_vectors
  USE lattice_type
  USE memory_handling
!  USE multipole_pbc
!  USE harmonics_pbc
  USE pbcffdata
  USE PBC_MSC
  USE pbc_interactions
  USE pbc_ff_contrib
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,lupri,luerr,numrealvec,nfsze,lmax
  COMPLEX(complexk), INTENT(INOUT) :: fock(ndim,ndim),Sabk(ndim,ndim)
  TYPE(lvec_list_t),INTENT(INOUT) :: ll
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule
  TYPE(LSSETTING) :: setting
  TYPE(BZgrid_t),intent(in) :: bz
  REAL(realk) ,intent(In) :: tlat((1+lmax)**2,(1+lmax)**2)
  !LOCAL VARIABLES
  !COMPLEX(complexk) :: fock_old(7,ndim,ndim),fock_vec(7,ndim*ndim)
  !COMPLEX(complexk) :: fockMO(7,ndim,ndim)
  REAL(realk) :: eigv(ndim),phase1,phase2,phase3
  COMPLEX(complexk) :: phase
  REAL(realk),allocatable :: error(:,:)
  !COMPLEX(COMPLEXK) :: U(ndim,ndim),C_tmp(7,ndim,ndim)!,fockMO(ndim,ndim)
  INTEGER :: i,j,errlm,tol!tol should be a real but now I just have it for test
  INTEGER :: k,kpt,nmax
  real(realk) :: kvec(3,3)
  REAL(realk),allocatable :: weight(:)
  real(realk),allocatable :: nucmom(:)
  TYPE(matrix), pointer :: nfdensity(:),F_1(:),ovl(:)
  TYPE(matrix), pointer :: g_2(:)
  TYPE(moleculeinfo) :: refcell
  TYPE(pbc_scfiterations_t) :: pbc_it(7)
  TYPE(pbc_elstr_t) :: kdep_tmp(bz%nk)
  TYPE(lvec_data_t),allocatable :: C_tmp_t(:)
  INTEGER :: scfit
  CHARACTER(LEN=10) :: numtostring1,numtostring2,numtostring3,iter
  character(len=15) :: mattxt,string1,filename
  integer :: l1,l2,l3,nil,fdim(3),iunit,n1,iunit2
  real(realk) :: fcao(ndim*ndim),E_kin,E_en,E_J,E_K,E_ff,E_nn,E_nuc
  TYPE(lvec_data_t) :: erikll(7)

  write(*,*) 'hei'
  fdim=0
  call mem_alloc(nfdensity,numrealvec)
  allocate(f_1(numrealvec))
  allocate(g_2(numrealvec))
  allocate(ovl(numrealvec))
  allocate(nucmom((lmax+1)**2))
  write(*,*) lmax


  call set_refcell(refcell,molecule)

  if (ll%ldef%is_active(1)) fdim(1) = 1
  if (ll%ldef%is_active(2)) fdim(2) = 1
  if (ll%ldef%is_active(3)) fdim(3) = 1
  
  iunit=1
  iunit2=2
  nmax=3
  do n1=1,numrealvec
     call mat_init(nfdensity(n1),ndim,ndim)
     call mat_zero(nfdensity(n1))
  enddo
    write(*,*) 'hei'

  k=5
  !DO k=2,2   !!!!!!!!!!12
   write(string1,'(I5)') k
   string1=adjustl(string1)
   DO l3=-ll%ndmat*fdim(3),ll%ndmat*fdim(3)
    DO l2=-ll%ndmat*fdim(2),ll%ndmat*fdim(2)
     DO l1=-ll%ndmat*fdim(1),ll%ndmat*fdim(1)

       write(numtostring1,'(I5)')  abs(l1)
       write(numtostring2,'(I5)')  abs(l2)
       write(numtostring3,'(I5)')  abs(l3)
       write(iter,'(I5)')  k
       iter=adjustl(iter)
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       if(l1 .ge. 0) then
         mattxt='dmt__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       else
         mattxt='dmt__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       endif

       call find_latt_index(n1,l1,l2,l3,ll,ll%max_layer)
       call pbc_readopmat2(mattxt,ll%lvec(n1)%d_mat,ndim,.true.,.false.)

       write(numtostring1,'(I5)')  l1
       write(numtostring2,'(I5)')  l2
       write(numtostring3,'(I5)')  l3
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       iter=adjustl(iter)
 !     numstring
       mattxt='erikdmt'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
       write(iunit,*) ndim
       DO j=1,ndim
         write(iunit,*) (ll%lvec(n1)%d_mat(i,j),i=1,ndim)
       ENDDO
       CALL lsCLOSE(IUNIT,'KEEP')
       Call mat_set_from_full(ll%lvec(n1)%d_mat,1.D0,nfdensity(n1))
   
!
     ENDDO !l1
    ENDDO  !l2
   ENDDO   !l3

!   DO n1=1,nfsze
!         write(lupri,*) 'check n1', n1,nfsze
!         CALL mat_print(nfdensity(n1),1,ndim,1,ndim,lupri)
!   enddo
#ifdef DEBUGPBC

#else

   k=5
   DO l3=-ll%max_layer*fdim(3),ll%max_layer*fdim(3)
    DO l2=-ll%max_layer*fdim(2),ll%max_layer*fdim(2)
     DO l1=-ll%max_layer*fdim(1),ll%max_layer*fdim(1)

       write(numtostring1,'(I5)')  abs(l1)
       write(numtostring2,'(I5)')  abs(l2)
       write(numtostring3,'(I5)')  abs(l3)
       write(iter,'(I5)')  k
       iter=adjustl(iter)
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       if (abs(l1) .le. 3) then
       if(l1 .ge. 0) then
         mattxt='xch__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
         !mattxt='diisfck0'//trim(iter)//'_p0'//trim(numtostring1)//'_p0'//trim(numtostring2)//'_p0'//trim(numtostring3)
       else
         mattxt='xch__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
         !mattxt='diisfck0'//trim(iter)//'_n0'//trim(numtostring1)//'_p0'//trim(numtostring2)//'_p0'//trim(numtostring3)
       endif

       call pbc_readopmat2(mattxt,ll%lvec(n1)%d_mat,ndim,.true.,.false.)

       write(numtostring1,'(I5)')  l1
       write(numtostring2,'(I5)')  l2
       write(numtostring3,'(I5)')  l3
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       iter=adjustl(iter)
 !     numstring
       mattxt='erikxch'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       !mattxt='erikfck'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
       !write(iunit,*) ndim
       DO j=1,ndim
         write(iunit,*) (ll%lvec(n1)%d_mat(i,j),i=1,ndim)
       ENDDO
       CALL lsCLOSE(IUNIT,'KEEP')
     endif
   
!
     ENDDO !l1
    ENDDO  !l2
   ENDDO   !l3
#endif

!   do l1=-3,3
!   if (abs(l1) .le. 3) then
!       write(numtostring1,'(I5)')  l1
!       numtostring1=adjustl(numtostring1)
!     mattxt='minham'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
!     CALL lsOPEN(IUNIT,mattxt,'OLD','FORMATTED')
!    call find_latt_index(k,l1,0,0,ll,ll%max_layer)
!     DO j=1,ndim
!       read(iunit,*) (ll%lvec(k)%fck_vec(i+(j-1)*ndim),i=1,ndim)
!     ENDDO
!     CALL lsCLOSE(IUNIT,'KEEP')
!   endif
! enddo

 do kpt=1,bz%nk
   call init_pbc_elstr(kdep_tmp(kpt),ndim,ndim)
 enddo
  
   call pbc_overlap_k(lupri,luerr,setting,molecule%natoms,ndim,ll,&
                    & refcell,numrealvec,ovl)

   call pbc_kinetic_k(lupri,luerr,setting,molecule%natoms,ndim,&
   & ll,refcell,numrealvec,nfdensity,f_1)

   !ll%nf=5

   call pbc_nucattrc_k(lupri,luerr,setting,molecule%natoms,ndim,&
     & ll,refcell,numrealvec,nfdensity,f_1)

   call pbc_exact_xc_k(lupri,luerr,setting,molecule%natoms,ndim,&
     & ll,refcell,numrealvec,nfdensity,g_2)


  call pbc_electron_rep_k(lupri,luerr,setting,molecule%natoms,ndim,&
     & ll,refcell,numrealvec,nfdensity,g_2)

   !This is needed to form fck
  call pbc_comp_nucmom(refcell,nucmom,lmax,lupri)

  !fixme nucmom 
!   call lsquit('fixme nucmom no value assigned to this variable',-1)
   call pbc_ff_fck(ll%tlmax,tlat,ll%lmax,ndim,ll,nfdensity,nucmom,&
                   & g_2,E_ff,E_nn,lupri)

  CALL pbc_nucpot(lupri,luerr,setting,molecule%natoms,ll,&
                  & refcell,numrealvec,E_nuc)
#ifdef DEBUGPBC
   do n1=1,numrealvec
       call mat_free(nfdensity(n1))
   enddo
   call mem_dealloc(nfdensity)
   write(*,*) 'Nuclear N.F. repulsion', E_nuc
   call LSQUIT('Finished computing matrices with erik',lupri)
#endif



   do n1=1,nfsze
      call mat_zero(nfdensity(n1))
   enddo


  !ENDDO !k

   deallocate(nucmom)
do n1=-3,3
    write(numtostring1,'(I5)') n1
    numtostring1=adjustl(numtostring1)
    mattxt='minEFmat3'//trim(numtostring1)//'00.dat'
      !call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
    !CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
    !call find_latt_index(k,n1,0,0,ll,ll%max_layer)
    !!write(iunit,*) ndim
    !DO j=1,ndim
    !   write(iunit,*) (ll%lvec(k)%fck_vec(i+(j-1)*ndim),i=1,ndim)
    !ENDDO
    !call lsclose(iunit,'KEEP')
 enddo
    write(*,*) 'Fock matrix written to disk'
    write(lupri,*) 'Fock matrix written to disk'

 do kpt=1,bz%nk
!    kdep_tmp(kpt)%kfockvec=0.0_realk
! do n1=-9,9
!
!    call pbc_get_kpoint(m,kvec)
!    call find_latt_index(k,n1,0,0,ll,ll%max_layer)
!    phase1=kvec(1)*ll%lvec(k)%std_coord(1)
!    phase2=kvec(2)*ll%lvec(k)%std_coord(2)
!    phase3=kvec(3)*ll%lvec(k)%std_coord(3)
!    phase=CMPLX(0.,(phase1+phase2+phase3),complexk)
!
!    kdep_tmp(kpt)%kfcockvec=kdep_tmp(kpt)%kfcockvec+ll%lvec(k)%fck_vec*exp(phase)
! 
! enddo
 
!    call pbc_zdevectorize_mat(kdep_tmp(kpt)%kfockmat,ndim,ndim,&
!         kdep_tmp(kpt)%kfockvec)
!
!    call solve_kfcsc_mat(kpt,ndim,kdep_tmp(kpt)%kfockmat,&
!    kdep(kpt)%koverlapmat,kdep_tmp(kpt)%kcdensitymat,eigv,lupri)
!!    kdep_tmp(kpt)%keigv(:)=eigv(:)
!    write(*,*) 'eig',eigv

!    
 enddo
! write(*,*) kdep_tmp(1)%kfockvec
! write(*,*) 'mat'
! write(*,*) kdep_tmp(1)%kfockmat
  !call COMPARE_MATRICES(lupri,ndim,numrealvec,nfsze,lmax,ll)
!  DO k=1,3
!  DO l1=-3,3
!   DO l2=0,0
!    DO l3=0,0
!       write(numtostring1,'(I5)')  abs(l1)
!       write(numtostring2,'(I5)')  abs(l2)
!       write(numtostring3,'(I5)')  abs(l3)
!       write(iter,'(I5)')  k
!       iter=adjustl(iter)
!       numtostring1=adjustl(numtostring1)
!       numtostring2=adjustl(numtostring2)
!       numtostring3=adjustl(numtostring3)
!       if(l1 .ge. 0) then
!         mattxt='kin__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
!       else
!         mattxt='kin__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
!       endif
!
!       !call pbc_readopmat2(mattxt,erikll(k)%d_mat,ndim,.true.,.false.)
!
!       write(numtostring1,'(I5)')  l1
!       write(numtostring2,'(I5)')  l2
!       write(numtostring3,'(I5)')  l3
!       numtostring1=adjustl(numtostring1)
!       numtostring2=adjustl(numtostring2)
!       numtostring3=adjustl(numtostring3)
!       iter=adjustl(iter)
! !     numstring
!       mattxt='erikin'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
!       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
!       write(iunit,*) ndim
!       !DO j=1,ndim
!       !  write(iunit,*) (erikll(k)%d_mat(i,j),i=1,ndim)
!       !ENDDO
!       CALL lsCLOSE(IUNIT,'KEEP')
!
!       write(numtostring1,'(I5)')  abs(l1)
!       write(numtostring2,'(I5)')  abs(l2)
!       write(numtostring3,'(I5)')  abs(l3)
!       numtostring1=adjustl(numtostring1)
!       numtostring2=adjustl(numtostring2)
!       numtostring3=adjustl(numtostring3)
!       if(l1 .ge. 0) then
!         mattxt='nnf'//trim(iter)//'__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
!       else
!         mattxt='nnf'//trim(iter)//'__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
!       endif
!
!       call pbc_readopmat2(mattxt,erikll(k)%d_mat,ndim,.true.,.false.)
!
!       write(numtostring1,'(I5)')  l1
!       write(numtostring2,'(I5)')  l2
!       write(numtostring3,'(I5)')  l3
!       numtostring1=adjustl(numtostring1)
!       numtostring2=adjustl(numtostring2)
!       numtostring3=adjustl(numtostring3)
! !     numstring
!       mattxt='eriknnf'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
!       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
!       write(iunit,*) ndim
!       DO j=1,ndim
!         write(iunit,*) (erikll(k)%d_mat(i,j),i=1,ndim)
!       ENDDO
!       CALL lsCLOSE(IUNIT,'KEEP')
!
!       write(numtostring1,'(I5)')  abs(l1)
!       write(numtostring2,'(I5)')  abs(l2)
!       write(numtostring3,'(I5)')  abs(l3)
!       numtostring1=adjustl(numtostring1)
!       numtostring2=adjustl(numtostring2)
!       numtostring3=adjustl(numtostring3)
!       if(l1 .ge. 0) then
!         mattxt='jop'//trim(iter)//'__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
!       else
!         mattxt='jop'//trim(iter)//'__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
!       endif
!
!       call pbc_readopmat2(mattxt,erikll(k)%d_mat,ndim,.true.,.false.)
!
!       write(numtostring1,'(I5)')  l1
!       write(numtostring2,'(I5)')  l2
!       write(numtostring3,'(I5)')  l3
!       numtostring1=adjustl(numtostring1)
!       numtostring2=adjustl(numtostring2)
!       numtostring3=adjustl(numtostring3)
! !     numstring
!       mattxt='erikjop'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
!       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
!       write(iunit,*) ndim
!       DO j=1,ndim
!         write(iunit,*) (erikll(k)%d_mat(i,j),i=1,ndim)
!       ENDDO
!       CALL lsCLOSE(IUNIT,'KEEP')
!
!       write(numtostring1,'(I5)')  abs(l1)
!       write(numtostring2,'(I5)')  abs(l2)
!       write(numtostring3,'(I5)')  abs(l3)
!       write(iter,'(I5)')  scfit
!       numtostring1=adjustl(numtostring1)
!       numtostring2=adjustl(numtostring2)
!       numtostring3=adjustl(numtostring3)
!       if(l1 .ge. 0) then
!         mattxt='xch'//trim(iter)//'__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
!       else
!         mattxt='xch'//trim(iter)//'__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
!       endif
!
!       call pbc_readopmat2(mattxt,erikll(k)%d_mat,ndim,.true.,.false.)
!
!       write(numtostring1,'(I5)')  l1
!       write(numtostring2,'(I5)')  l2
!       write(numtostring3,'(I5)')  l3
!       numtostring1=adjustl(numtostring1)
!       numtostring2=adjustl(numtostring2)
!       numtostring3=adjustl(numtostring3)
! !     numstring
!       mattxt='erikxch'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
!       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
!       write(iunit,*) ndim
!       DO j=1,ndim
!         write(iunit,*) (erikll(k)%d_mat(i,j),i=1,ndim)
!       ENDDO
!       CALL lsCLOSE(IUNIT,'KEEP')
!
!    enddo
!   enddo
!  enddo
!  enddo

  
  

  
END SUBROUTINE readerikmats 

SUBROUTINE COMPARE_MATRICES(lupri,ndim,numrealvec,nfsze,lmax,ll)
  USE TYPEDEF
  USE precision
  USE matrix_module
  USE lattice_vectors
  USE lattice_type
!  USE multipole_pbc
!  USE harmonics_pbc
  USE pbcffdata
  USE PBC_MSC
  USE pbc_interactions
  USE pbc_ff_contrib
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,lupri,numrealvec,nfsze,lmax
  TYPE(lvec_list_t),INTENT(IN) :: ll
  !LOCAL VARIABLES
  INTEGER :: i,j !,toltol should be a real but now I just have it for test
  !TYPE(moleculeinfo) :: refcell
  INTEGER,save :: scfit=0
  CHARACTER(LEN=10) :: numtostring1,numtostring2,numtostring3,iter
  character(len=20) :: mattxt,string1
  integer     ::   iunit,k,fdim(3)!,ldum
  integer :: l1,l2,l3,nbast
  TYPE(lvec_data_t) :: erikll(7)

  !PE: nbast was not set and the compiler kept complaining, I set it to Zero
  !here, the responsible person should make sure the correct value is set. Even
  !though there were several complaints through the mailing list this was not
  !fixed by the responsible person. KEEP THE CODE CLEAN
  nbast = 0

  ! write fmat to file
  scfit=scfit+1
  l2=0
  l3=0
  k=0
  write(string1,'(I5)') scfit
  string1=adjustl(string1)
  do l1=-3,3
  do l2=0,0
  do l3=0,0
  k=k+1
  iunit = -1
  write(numtostring1,'(I5)') l1
  write(numtostring2,'(I5)') l2
  write(numtostring3,'(I5)') l3
  numtostring1=adjustl(numtostring1)
  numtostring2=adjustl(numtostring2)
  numtostring3=adjustl(numtostring3)

  write(mattxt,'(A20)') 'PBCFMAT'//trim(string1)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)

  mattxt=adjustl(mattxt)
  mattxt=trim(mattxt)
  !write(*,*) mattxt

  !CALL lsOPEN(IUNIT,mattxt,'old','UNFORMATTED')
  !read(iunit) nbast
  !allocate(erikll(k)%fck_vec(nbast*nbast))
  !allocate(erikll(k)%d_mat(nbast,nbast))
  !DO j=1,nbast
  !  read(iunit) (erikll(k)%fck_vec(i+(j-1)*nbast),i=1,nbast)
  !ENDDO
  !CALL lsCLOSE(IUNIT,'KEEP')
  enddo
  enddo
  enddo

  !printing erik's fock to file in readable format
!  k=0
!  do l1=-3,3
!  do l2=0,0
!  do l3=0,0
!  iunit = -1
!  k=k+1
!  write(numtostring1,'(I5)')  l1
!  write(numtostring2,'(I5)')  l2
!  write(numtostring3,'(I5)')  l3
!  write(iter,'(I5)')  scfit
!  numtostring1=adjustl(numtostring1)
!  numtostring2=adjustl(numtostring2)
!  numtostring3=adjustl(numtostring3)
!  iter=adjustl(iter)
! !numstring
!  mattxt='erikFMAT'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
!  CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
!  write(iunit,*) nbast
!  DO j=1,nbast
!    write(iunit,*) (erikll(k)%fck_vec(i+(j-1)*nbast),i=1,nbast)
!  ENDDO
!  CALL lsCLOSE(IUNIT,'KEEP')
!  enddo
!  enddo
!  enddo

   !printing may fock to file in readable format
  k=0
  do l1=-3,3
  do l2=0,0
  do l3=0,0
  iunit = -1
  k=k+1
  write(numtostring1,'(I5)')  l1
  write(numtostring2,'(I5)')  l2
  write(numtostring3,'(I5)')  l3
  write(iter,'(I5)')  scfit
  numtostring1=adjustl(numtostring1)
  numtostring2=adjustl(numtostring2)
  numtostring3=adjustl(numtostring3)
  iter=adjustl(iter)
 !numstring
  mattxt='minFMAT'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
  !CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
  !write(iunit,*) nbast
  !DO j=1,nbast
  !  write(iunit,*) (ll%lvec(k)%fck_vec(i+(j-1)*nbast),i=1,nbast)
  !ENDDO
  !CALL lsCLOSE(IUNIT,'KEEP')
  enddo
  enddo
  enddo
 
  do l1=0,0
  do l2=0,0
  do l3=0,0
  iunit = -1
  call find_latt_index(k,l1,l2,l3,ll,3)
  write(numtostring1,'(I5)')  l1
  write(numtostring2,'(I5)')  l2
  write(numtostring3,'(I5)')  l3
  write(iter,'(I5)')  scfit
  numtostring1=adjustl(numtostring1)
  numtostring2=adjustl(numtostring2)
  numtostring3=adjustl(numtostring3)
  iter=adjustl(iter)
 !numstring
  mattxt='minDMAT'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
  CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
  write(iunit,*) nbast
  DO j=1,nbast
    write(iunit,*) (ll%lvec(k)%d_mat(i,j),i=1,nbast)
  ENDDO
  CALL lsCLOSE(IUNIT,'KEEP')
  enddo
  enddo
  enddo

  l1=0
  l2=0
  l3=0
  !do l1=-3,3
  iunit = -1
!  k=1
  write(numtostring1,'(I5)')  l1
  write(numtostring2,'(I5)')  l2
  write(numtostring3,'(I5)')  l3
  numtostring1=adjustl(numtostring1)
  numtostring2=adjustl(numtostring2)
  numtostring3=adjustl(numtostring3)
  write(mattxt,'(A20)') 'PBCDMAT'//trim(string1)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)

  mattxt=adjustl(mattxt)
  mattxt=trim(mattxt)
  !write(*,*) mattxt

   CALL lsOPEN(IUNIT,mattxt,'old','UNFORMATTED')
   read(iunit) nbast
   DO j=1,nbast
     read(iunit) (erikll(k)%d_mat(i,j),i=1,nbast)
   ENDDO
      CALL lsCLOSE(IUNIT,'KEEP')
  !enddo
  
  write(numtostring1,'(I5)')  l1
  write(numtostring2,'(I5)')  l2
  write(numtostring3,'(I5)')  l3
  write(iter,'(I5)')  scfit
  numtostring1=adjustl(numtostring1)
  numtostring2=adjustl(numtostring2)
  numtostring3=adjustl(numtostring3)
  iter=adjustl(iter)
 !numstring
  mattxt='erikDMAT'//trim(iter)//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
  CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
  write(iunit,*) nbast
  DO j=1,nbast
    write(iunit,*) (erikll(k)%d_mat(i,j),i=1,nbast)
  ENDDO
  CALL lsCLOSE(IUNIT,'KEEP')

  k=0
  DO l1=-3,3
   DO l2=0,0
    DO l3=0,0
       k=k+1
       write(numtostring1,'(I5)')  abs(l1)
       write(numtostring2,'(I5)')  abs(l2)
       write(numtostring3,'(I5)')  abs(l3)
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       if(l1 .ge. 0) then
         mattxt='kin__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       else
         mattxt='kin__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       endif

       call pbc_readopmat2(mattxt,erikll(k)%d_mat,nbast,.true.,.false.)

       write(numtostring1,'(I5)')  l1
       write(numtostring2,'(I5)')  l2
       write(numtostring3,'(I5)')  l3
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       iter=adjustl(iter)
 !     numstring
       mattxt='erikin3'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
       write(iunit,*) nbast
       DO j=1,nbast
         write(iunit,*) (erikll(k)%d_mat(i,j),i=1,nbast)
       ENDDO
       CALL lsCLOSE(IUNIT,'KEEP')

       write(numtostring1,'(I5)')  abs(l1)
       write(numtostring2,'(I5)')  abs(l2)
       write(numtostring3,'(I5)')  abs(l3)
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       if(l1 .ge. 0) then
         mattxt='ham__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       else
         mattxt='ham__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       endif

       call pbc_readopmat2(mattxt,erikll(k)%d_mat,nbast,.true.,.false.)

       write(numtostring1,'(I5)')  l1
       write(numtostring2,'(I5)')  l2
       write(numtostring3,'(I5)')  l3
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       iter=adjustl(iter)
 !     numstring
       mattxt='erikham3'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
       write(iunit,*) nbast
       DO j=1,nbast
         write(iunit,*) (erikll(k)%d_mat(i,j),i=1,nbast)
       ENDDO
       CALL lsCLOSE(IUNIT,'KEEP')

       write(numtostring1,'(I5)')  abs(l1)
       write(numtostring2,'(I5)')  abs(l2)
       write(numtostring3,'(I5)')  abs(l3)
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       if(l1 .ge. 0) then
         mattxt='jop__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       else
         mattxt='jop__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       endif

       call pbc_readopmat2(mattxt,erikll(k)%d_mat,nbast,.true.,.false.)

       write(numtostring1,'(I5)')  l1
       write(numtostring2,'(I5)')  l2
       write(numtostring3,'(I5)')  l3
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       iter=adjustl(iter)
 !     numstring
       mattxt='erikjop3'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
       write(iunit,*) nbast
       DO j=1,nbast
         write(iunit,*) (erikll(k)%d_mat(i,j),i=1,nbast)
       ENDDO
       CALL lsCLOSE(IUNIT,'KEEP')

       write(numtostring1,'(I5)')  abs(l1)
       write(numtostring2,'(I5)')  abs(l2)
       write(numtostring3,'(I5)')  abs(l3)
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       if(l1 .ge. 0) then
         mattxt='xch__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       else
         mattxt='xch__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       endif

       call pbc_readopmat2(mattxt,erikll(k)%d_mat,nbast,.true.,.false.)

       write(numtostring1,'(I5)')  l1
       write(numtostring2,'(I5)')  l2
       write(numtostring3,'(I5)')  l3
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       iter=adjustl(iter)
 !     numstring
       mattxt='erikxch3'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
       write(iunit,*) nbast
       DO j=1,nbast
         write(iunit,*) (erikll(k)%d_mat(i,j),i=1,nbast)
       ENDDO
       CALL lsCLOSE(IUNIT,'KEEP')
       write(numtostring1,'(I5)')  abs(l1)
       write(numtostring2,'(I5)')  abs(l2)
       write(numtostring3,'(I5)')  abs(l3)
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)


if(l1 .ge. 0) then
         mattxt='nnf__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       else
         mattxt='nnf__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       endif

       call pbc_readopmat2(mattxt,erikll(k)%d_mat,nbast,.true.,.false.)

       write(numtostring1,'(I5)')  l1
       write(numtostring2,'(I5)')  l2
       write(numtostring3,'(I5)')  l3
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       iter=adjustl(iter)
 !     numstring
       mattxt='eriknnf3'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
       write(iunit,*) nbast
       DO j=1,nbast
         write(iunit,*) (erikll(k)%d_mat(i,j),i=1,nbast)
       ENDDO
       CALL lsCLOSE(IUNIT,'KEEP')

       write(numtostring1,'(I5)')  abs(l1)
       write(numtostring2,'(I5)')  abs(l2)
       write(numtostring3,'(I5)')  abs(l3)
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)


if(l1 .ge. 0) then
         mattxt='fck__p'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       else
         mattxt='fck__n'//trim(numtostring1)//'__p'//trim(numtostring2)//'__p'//trim(numtostring3)
       endif

       call pbc_readopmat2(mattxt,erikll(k)%d_mat,nbast,.true.,.false.)

       write(numtostring1,'(I5)')  l1
       write(numtostring2,'(I5)')  l2
       write(numtostring3,'(I5)')  l3
       write(iter,'(I5)')  scfit
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       iter=adjustl(iter)
 !     numstring
       mattxt='erikfck1'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
       write(iunit,*) nbast
       DO j=1,nbast
         write(iunit,*) (erikll(k)%d_mat(i,j),i=1,nbast)
       ENDDO
       CALL lsCLOSE(IUNIT,'KEEP')



    enddo
   enddo
  enddo

!  call pbc_readopmat2('ham__p1__p0__p0',erikll(k)%d_mat,nbast,.true.,.false.)
!  call pbc_readopmat2('ham__p2__p0__p0',erikll(k)%d_mat,nbast,.true.,.false.)
!  call pbc_readopmat2('ham__p3__p0__p0',erikll(k)%d_mat,nbast,.true.,.false.)
!  call pbc_readopmat2('ham__n1__p0__p0',erikll(k)%d_mat,nbast,.true.,.false.)
!  call pbc_readopmat2('ham__n2__p0__p0',erikll(k)%d_mat,nbast,.true.,.false.)
!  call pbc_readopmat2('ham__n3__p0__p0',erikll(k)%d_mat,nbast,.true.,.false.)

!  iunit = -1
!  k=1
!   write(mattxt,'(A11)') 'erikjop_000'
!   CALL lsOPEN(IUNIT,mattxt,'UNKNOWN','FORMATTED')
!   write(iunit,*) nbast
!   DO j=1,nbast
!     write(iunit,*) (erikll(k)%d_mat(i,j),i=1,nbast)
!   ENDDO
!      CALL lsCLOSE(IUNIT,'KEEP')


!  call pbc_readopmat2('ovl__p0__p0__p0',erikll(k)%d_mat,nbast,.true.,.false.)
!
!  k=1
!   write(mattxt,'(A11)') 'erikovl_000'
!   CALL lsOPEN(IUNIT,mattxt,'UNKNOWN','FORMATTED')
!   write(iunit,*) nbast
!   DO j=1,nbast
!     write(iunit,*) (erikll(k)%d_mat(i,j),i=1,nbast)
!   ENDDO
!      CALL lsCLOSE(IUNIT,'KEEP')
!
!  iunit = -1
!  k=1
!  call find_latt_index(nil,0,0,0,ll,ll%max_layer)
!  write(*,*) nbast,ndim
!   write(mattxt,'(A6,4I2)') 'minjop',scfit,0,l2,l3
!   CALL lsOPEN(IUNIT,mattxt,'UNKNOWN','FORMATTED')
!   write(iunit,*) nbast!,nil
!   DO j=1,nbast
!     write(iunit,*) (ll%lvec(nil)%fck_vec(i+(j-1)*nbast),i=1,nbast)
!   ENDDO
!      CALL lsCLOSE(IUNIT,'KEEP')
! !enddo
!   write(mattxt,'(A7,4I2)') 'PBC1MAT',scfit,l1,l2,l3
!   CALL lsOPEN(IUNIT,mattxt,'old','UNFORMATTED')
!   read(iunit) nbast
!   DO j=1,nbast
!     read(iunit) (erikll(k)%fck_vec(i+(j-1)*nbast),i=1,nbast)
!   ENDDO
!      CALL lsCLOSE(IUNIT,'KEEP')
!
!   write(mattxt,'(A7,4I2)') 'erik1MAT',scfit,l1,l2,l3
!   CALL lsOPEN(IUNIT,mattxt,'UNKNOWN','FORMATTED')
!   write(iunit,*) nbast
!   DO j=1,nbast
!     write(iunit,*) (erikll(k)%fck_vec(i+(j-1)*nbast),i=1,nbast)
!   ENDDO
!      CALL lsCLOSE(IUNIT,'KEEP')







END SUBROUTINE COMPARE_MATRICES
