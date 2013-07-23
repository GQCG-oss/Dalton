#ifdef MOD_UNRELEASED
MODULE pbc_scfdiis
  USE TYPEDEF
  USE precision
  USE matrix_module
  USE lattice_vectors
  USE lattice_type
  USE multipole_pbc
  USE harmonics_pbc
  USE pbc_matrix_operations
  USE pbc_interactions
  USE pbcffdata
  USE PBC_MSC
  USE PBC_kspce_rspc_operations
  USE pbc_ff_contrib

  CONTAINS
SUBROUTINE pbc_zeigsolve(A,B,N,M,eigv,lupri)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,M,lupri
  COMPLEX(COMPLEXK), INTENT(INOUT) :: A(N,M),B(N,M)
  REAL(realk),INTENT(INOUT) :: eigv(N)
  !LOCAL VARIABLES
  COMPLEX(COMPLEXK), pointer :: work(:)
  !real(realk) :: rwork(8*n)
  real(realk),pointer :: rwork(:)
  INTEGER :: info,lwork
  INTEGER,save :: ncalls=0
  ncalls=ncalls+1
  lwork=2*n-1

!(c_tmp,Sabk_tmp,ndim,ndim,eigv,lupri)
!  allocate(work(lwork))
  call mem_alloc(work,lwork)
  call mem_alloc(rwork,3*N-2)
  
  !call zheev('V','U',N,A,N,eig,work,3*N+1,rwork,INFO)

  !if(info .ne. 0) THEN
  !  call LSQUIT('pbc_zeigsolve: INFO not zero, while solving eigenvalue',lupri)
  !endif

!  write(*,*) 'just calling to solve eigenvalue equation'
  call zhegv(1,'V','U',n,A,n,B,n,eigv,work,lwork,rwork,info)

  if(info .ne. 0) THEN
    write(lupri,*) 'ERROR: zhegv problems, info=', info
    call write_zmatrix(B,n,n)
    write(*,*) 'ERROR: zhegv problems, info=', info
    call LSQUIT('pbc_zeigsolve: INFO not zero, while solving eigenvalue',lupri)
  endif

  !deallocate(work)
  call mem_dealloc(work)
  call mem_dealloc(rwork)


END SUBROUTINE pbc_zeigsolve

SUBROUTINE pbc_zggeigsolve(kindex,A,B,smatk,N,M,eigv,lupri)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: kindex,N,M,lupri
  COMPLEX(COMPLEXK), INTENT(INOUT) :: A(N,M),B(N,M),smatk(N,M)
  REAL(realk),INTENT(INOUT) :: eigv(N)
  !LOCAL VARIABLES
  COMPLEX(COMPLEXK), pointer :: work(:)
  real(realk),pointer :: rwork(:)
  COMPLEX(COMPLEXK) :: alphavec(n),betavec(n)
  COMPLEX(COMPLEXK) :: vl(n,m),ccoeff(n,m)
  COMPLEX(COMPLEXK) :: ztemp
  INTEGER :: info,iband,i,j,lwork
  INTEGER,save :: ncalls=0
  ncalls=ncalls+1
  lwork=2*n*n
  

!(c_tmp,Sabk_tmp,ndim,ndim,eigv,lupri)
!  allocate(work(lwork))
  call mem_alloc(work,lwork)
  call mem_alloc(rwork,8*N)
  
  !write(lupri,*) 'Smat before energy'
  !call write_zmatrix(B,n,n,lupri)
  !write(lupri,*) 'Fockmat before energy'
  !call write_zmatrix(A,n,n,lupri)

    !write(*,*) 'debug after read tlatt',kindex
  call zggev('N','V',n,A,n,B,n,alphavec,betavec,vl,1,ccoeff,n,work,2*n,rwork,info)
    !write(*,*) 'debug after read tlatt',kindex

  if(info .ne. 0) THEN
    write(lupri,*) 'ERROR: zggev problems, info=', info
    write(*,*) 'ERROR: zggev problems, info=', info
    call LSQUIT('pbc_zggeigsolve: INFO not zero, while solving eigenvalue',lupri)
  endif

  !write(lupri,*) 'ccoeff before energy'
  !call write_zmatrix(ccoeff,n,n,lupri)

  !fix normalization which zggev does in a strange form
  call pbc_fixzggevnorm(N,ccoeff,smatk,alphavec,betavec,lupri)

  !write(lupri,*) 'ccoeff after gevnorm'
  !call write_zmatrix(ccoeff,n,n,lupri)
!  do i=1,n
!   do j=1,m
!    A(i,j) = ccoeff(i,j)
!   enddo
!  enddo

  ! we now have the crystal orbital coefficients in WORK(icocoeff_k)
  ! and the orbital energies as eps(j) = alpha(j)/beta(j)
  loop_bands: do iband = 1,n
     if (abs(alphavec(iband)) + abs(betavec(iband)) .lt. 1.0D-8) then
        ! because the virtual orbitals are used in the calculation of
        ! the DIIS error vector we need to actually zero the discarded
        ! orbital
        ccoeff(:,iband) = 0.0D0
        ztemp = 9.0D2
        write(LUPRI,*) 'Warning: Orbital at (b=',iband,', k=',kindex, &
             & ') projected out.'
     else if (abs(betavec(iband)) .le. 1.0D-3 * abs(alphavec(iband))) then
        ztemp = 1.0D3
        write(LUPRI,*) 'pbc_focksolver: Orbital energies given by ', &
             & 'alpha(j)/beta(j). Error: beta(j) very close to zero.'
        write(LUPRI,*) 'alpha = ',alphavec(iband)
        write(LUPRI,*) 'beta  = ',betavec(iband)
        !call LSquit('pbc_focksolver: Zero eigenvalue denominator.')
     else
        ztemp = alphavec(iband) / betavec(iband)
     end if

     eigv(iband) = real(ztemp)
#if 1
     if (.false.) then
#else
     if (eigv(iband) .lt. -1.0D3) then
        write(LUPRI,*) 'Warning: Flipping sign of eps(b=',iband,',k=', &
             & kindex,') = ',eigv(iband)
        eigv(iband) = abs(eigv(iband))
     else if (eigv(iband) .gt. 1.0D3) then
        write(LUPRI,*) 'Warning: eps(b=',iband,',k=',kindex,') = ', &
             & eigv(iband)
#endif
     else if (abs(aimag(ztemp)) .gt. 1.0D-5) then
        write(LUPRI,*) 'pbc_focksolver: Error, complex orbital energy.'
!        write(LUPRI,'(A,I5,A,3(F10.6,A))') '  k-point no. ',kindex, &
!             & ': kvec = (',kvec(1),',',kvec(2),',',kvec(3),')'
        write(LUPRI,*) '  Band index: ',iband
        write(LUPRI,*) '  Energy(band,kvec) = ',ztemp
        write(LUPRI,*) '  alpha(band,kvec) = ',alphavec(iband)
        write(LUPRI,*) '  beta(band,kvec)  = ',betavec(iband)
        write(LUPRI,*) '  call             = ', ncalls

        call LSquit('pbc_focksolver: Complex orbital energy.',lupri)
     end if
  end do loop_bands

  !deallocate(work)
  call mem_dealloc(work)
  call mem_dealloc(rwork)


  do i=1,n
   do j=1,m
    A(i,j) = ccoeff(i,j)
   enddo
  enddo

 

END SUBROUTINE pbc_zggeigsolve


SUBROUTINE pbc_deigsolve(A,N,M,eig,lupri)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,M,lupri
  REAL(realk), INTENT(INOUT) :: A(N,M),eig(N)
  !LOCAL VARIABLES
  real(realk) :: work(3*n-1)
  INTEGER :: info
  info=0
  call dsyev('V','U',N,A,N,eig,work,3*N-1,INFO)

  if(info .ne. 0) THEN
    call LSQUIT('pbc_deigsolve: INFO not zero, while solving eigenvalue',lupri)
  endif

END SUBROUTINE pbc_deigsolve


SUBROUTINE pbc_zdiagonalize(sigma,A,ndim,V,U)  
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim
  REAL(realk),INTENT(IN) :: sigma(ndim)
  COMPLEX(COMPLEXK),INTENT(INOUT) :: A(ndim,ndim),V(ndim,ndim)
  COMPLEX(COMPLEXK),INTENT(OUT) :: U(ndim,ndim)
  !LOCAL VARIABLES
  INTEGER :: i,j
  COMPLEX(COMPLEXK) :: D(ndim,ndim)
  COMPLEX(COMPLEXK) :: alpha,beta

  DO i=1,ndim
   DO j=1,ndim
    U(i,j) = V(i,j)/sqrt(sigma(j))
    D(i,j) =cmplx(0d0,0d0,complexk)
   ENDDO
  ENDDO

  alpha=cmplx(1D0,0d0,complexk)
  beta=cmplx(0D0,0d0,complexk)
  call zgemm('C','N',ndim,ndim,ndim,alpha,U,ndim,A,ndim,&
        beta,D,ndim)

  Write(*,*) 'U dagger fock'
  call write_zmatrix(D,ndim,ndim)

  call zgemm('N','N',ndim,ndim,ndim,alpha,D,ndim,U,ndim,&
        beta,A,ndim)

END SUBROUTINE pbc_zdiagonalize


SUBROUTINE ztransformbackC(C,ndim,U)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim
  COMPLEX(COMPLEXK),INTENT(INOUT) :: C(ndim,ndim)
  COMPLEX(COMPLEXK),INTENT(INOUT) :: U(ndim,ndim)
  COMPLEX(COMPLEXK) :: ctmp(ndim,ndim)
  COMPLEX(COMPLEXK) :: alpha,beta

  alpha=cmplx(1D0,0d0,complexk)
  beta =cmplx(0D0,0d0,complexk)
  ctmp(:,:)=C(:,:)
  call zgemm('N','N',ndim,ndim,ndim,alpha,U,ndim,ctmp,ndim,&
        beta,C,ndim)

END SUBROUTINE ztransformbackC

SUBROUTINE transform_toMOfock(C_tmp,fock,fockMO,ndim,lupri)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,lupri
  COMPLEX(COMPLEXK),INTENT(IN) :: C_tmp(ndim,ndim),fock(ndim,ndim)  
  COMPLEX(COMPLEXK),INTENT(OUT) :: fockMO(ndim,ndim)
  !LOCAL VARIABLES
  !INTEGER :: i,j,mu,nu
  COMPLEX(COMPLEXK) :: fockMOtmp(ndim,ndim)
  !LOCAL VARIABLES
  COMPLEX(COMPLEXK) :: alpha,beta

  alpha=CMPLX(1D0,0D0,complexk)
  beta =CMPLX(0D0,0D0,complexk)
  
  !FOR DEBUGGING
 ! write(*,*) 'C(1)' 
 ! call write_zmatrix(C_tmp,ndim,ndim)
 ! write(*,*) 'Fock(1)'
 ! call write_zmatrix(fock,ndim,ndim)
 fockMO=0d0
 fockMOtmp=0d0
 !write(*,*) 'C_tmp'
 !write(lupri,*) 'C_tmp'
 !call write_zmatrix(C_tmp,ndim,ndim,lupri)

! do i=1,ndim
!  do j=1,ndim
!   do mu=1,ndim
!    do nu=1,ndim
!    fockMO=fockMO(i,j)+c_tmp(i,mu)&
!                   *fock(mu,nu)*c_tmp(nu,j)
!    enddo
!   enddo
!  enddo
! enddo
  

  call zgemm('C','N',ndim,ndim,ndim,alpha,C_tmp,ndim,fock,ndim,&
        beta,fockMOtmp,ndim)


  call zgemm('N','N',ndim,ndim,ndim,alpha,fockMOtmp,ndim,C_tmp,ndim,&
        beta,fockMO,ndim)

 !write(*,*) 'fockMO'
 !write(lupri,*) 'fockMO'
 !call write_zmatrix(fockMO,ndim,ndim,lupri)
END SUBROUTINE transform_toMOfock

SUBROUTINE  pbc_geterrorvec(error,fockMO,ndim,errlm,nelectrons)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,errlm,nelectrons
  REAL(realk),INTENT(INOUT) :: error(errlm)
  COMPLEX(COMPLEXK),INTENT(IN) :: fockMO(ndim,ndim)
  !LOCAL variables
  INTEGER :: i,j,k

  k=0

!  write(*,*) 'DEBUG inside geterror'
!  call write_zmatrix(fockMO,ndim,ndim)
!  write(*,*) 'error vector'
  DO i=1,nelectrons/2
   DO j=nelectrons/2+1,ndim
    k=k+1
    error(k)=abs(real(fockMO(i,j),realk))
    !write(*,*) error(k)
   ENDDO
 ENDDO
! write(*,*) 'fockMO'
! call write_zmatrix(fockMO,ndim,ndim)
 !stop
      

END SUBROUTINE  pbc_geterrorvec

SUBROUTINE pbc_kcomputeenergy()
  IMPLICIT NONE
END SUBROUTINE pbc_kcomputeenergy

SUBROUTINE pbc_dcomputeenergy()
  IMPLICIT NONE
END SUBROUTINE pbc_dcomputeenergy


SUBROUTINE solve_kfcsc_mat(kindex,ndim,fock_old,Sabk,C_tmp,eigv,lupri)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,lupri,kindex
  COMPLEX(COMPLEXK),INTENT(IN) :: fock_old(ndim,ndim)
  COMPLEX(COMPLEXK), INTENT(IN) :: Sabk(ndim,ndim)
  COMPLEX(COMPLEXK),intent(INOUT) :: C_tmp(ndim,ndim)
  REAL(realk),intent(INOUT) :: eigv(ndim)
  !LOCAL VARIABLES
  INTEGER :: i,j
  COMPLEX(COMPLEXK) :: Sabk_tmp(ndim,ndim)


  DO i=1,ndim
   eigv(i)=0.0d0
   DO j=1,ndim
   C_tmp(i,j)=fock_old(i,j)
   Sabk_tmp(i,j)=Sabk(i,j)
  ! fock_old(1,i,j)=fock(i,j)
   ENDDO
  ENDDO
         
    !call write_zmatrix(sabk_tmp,ndim,ndim,lupri)
    call pbc_zeigsolve(c_tmp,Sabk_tmp,ndim,ndim,eigv,lupri)
    !call pbc_zggeigsolve(kindex,c_tmp,Sabk_tmp,Sabk,ndim,ndim,eigv,lupri)

 ! write(*,*) 'fock before transformation'
!  call write_zmatrix(C_tmp,ndim,ndim)

  !creates U matrix and U^FU matrix in C_tmp
!  call pbc_zdiagonalize(eigv,C_tmp,ndim,Sabk_tmp,U)

!  write(*,*) 'Transformed fock'
!  call write_zmatrix(C_tmp,ndim,ndim)
!  write(*,*) 'U matrix'
!  call write_zmatrix(U,ndim,ndim)

!  call pbc_zeigsolve(C_tmp,ndim,ndim,eigv,lupri)

!  call ztransformbackC(C_tmp,ndim,U)
 
!  write(*,*) 'k dependent density matrix'
!  call write_zmatrix(C_tmp,ndim,ndim)

END SUBROUTINE solve_kfcsc_mat

! Mady by Johannes 
!SUBROUTINE FOR COMPUTING energy for a k vector
SUBROUTINE pbc_k_energy(Aop,nrows,ncols,Bz,energy_1k,energy_2k,kvec,lupri)
IMPLICIT NONE
TYPE(lvec_list_t) :: Aop
INTEGER,INTENT(IN) :: lupri
INTEGER,INTENT(IN) :: nrows,ncols
real(realk),intent(in) :: kvec(3)
COMPLEX(COMPLEXK),INTENT(INOUT) :: energy_1k(nrows,ncols),energy_2k(nrows,ncols)
TYPE(BZgrid_t),intent(inout) :: bz
character(len=23)  :: diis
INTEGER            :: k,nlats

  !COMPUTE THE energy for k

  nlats=size(Aop%lvec)

     call zero_pbc_elstr(Bz%fck)
     call zero_pbc_elstr(Bz%Smat)

     call pbc_rspc_to_kspc_mat(Aop,Bz,nrows,kvec,2)
     call lsquit('fixme diis no value assigned to this variable ',-1)
!     call pbc_get_onep_matrix(Aop,energy_1k,nrows,ncols,nlats,kvec,diis)

!     call pbc_get_twop_matrix(Aop,energy_2k,nrows,ncols,nlats,kvec,diis)



END SUBROUTINE pbc_k_energy

! Mady by Johannes 
!SUBROUTINE FOR COMPUTING energy for one cell
SUBROUTINE pbc_cell_energy(Aop,Bz,nrows,ncols,Dk,kvec,energy_cell,kpt,nocc,lupri)
IMPLICIT NONE
INTEGER,INTENT(IN) :: lupri,nrows,ncols,kpt,nocc
TYPE(lvec_list_t) :: Aop
TYPE(BZgrid_t),intent(inout) :: Bz
real(realk),intent(in) :: kvec(3)
real(realk),intent(inout) ::energy_cell
COMPLEX(COMPLEXK):: energy1_k(nrows,ncols),energy2_k(nrows,ncols)
COMPLEX(COMPLEXK):: Dk(nrows,ncols),dijkl
!LOCAL VARIABLES
INTEGER :: i,j,k,l
REAL(realk) :: onepenergy,twopenergy,energy_k



call pbc_k_energy(Aop,nrows,ncols,Bz,energy1_k,energy2_k,kvec,lupri)

 onepenergy=0._realk
 twopenergy=0._realk
 DO i=1,nocc
 onepenergy=onepenergy+2.*bz%keigv((kpt-1)*nrows+i)
 enddo
 DO i=1,nrows
  DO j=1,ncols
  !onepenergy=onepenergy+real(energy1_k(i,j)*Dk(j,i),realk)
  twopenergy=twopenergy+0.5*real(energy2_k(i,j)*Dk(j,i),realk)
   !DO k=1,nrows
   ! DO l=1,ncols
   !  dijkl=DK(i,j)*Dk(k,l)-0.5*Dk(i,l)*Dk(k,j)
   !  twopenergy=twopenergy+0.5*

   ! ENDDO
   !ENDDO
  ENDDO
 ENDDO

 energy_k=onepenergy-twopenergy
 write(lupri,*) 'energy_k',onepenergy,twopenergy
 energy_cell=energy_cell+real(energy_k*bz%kpnt(kpt)%weight/BZ%NK_nosym,realk)
 


END SUBROUTINE pbc_cell_energy

!Made by Johannes
!SUBROUTINE for the HF iterations, 
!solves the equation F(k)C(k)=eps(k)C(k) for each k
SUBROUTINE pbc_startzdiis(molecule,setting,ndim,ll,numrealvec,&
           nfsze,maxmultmom,bz,ovl,f_1,g_2,E_1,lupri,luerr)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,lupri,luerr,numrealvec,nfsze,maxmultmom
  TYPE(lvec_list_t),INTENT(INOUT) :: ll
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule
  TYPE(LSSETTING) :: setting
  TYPE(BZgrid_t),intent(inout) :: bz
  REAL(realk),INTENT(INOUT) :: E_1
  TYPE(matrix),target,intent(inout) :: f_1(numrealvec),ovl(numrealvec)
  TYPE(matrix),target,intent(inout) :: g_2(numrealvec)
  !LOCAL VARIABLES
  COMPLEX(COMPLEXK) :: fockMO(7,ndim,ndim),fock(ndim,ndim),smatk(ndim,ndim)
  COMPLEX(COMPLEXK) :: C_k(ndim,ndim),D_k(ndim,ndim),C_0(ndim,ndim)
  REAL(realk) :: eigv(ndim),cellenergies(ndim),errortest,dmat(ndim,ndim)
  REAL(realk),pointer :: error(:,:),nucmom(:)
  INTEGER :: i,j,errlm,tol!tol should be a real but now I just have it for test
  INTEGER :: k,kpt,n1,iunit,fdim(3),indexx,layer
  INTEGER :: il1,il2,il3,realcut(3)
  real(realk) :: kvec(3),Ecell,E_nn
  REAL(realk),pointer :: tlat(:,:),weight(:)
  TYPE(matrix), pointer :: nfdensity(:)
  TYPE(moleculeinfo),pointer :: latt_cell(:)
  TYPE(moleculeinfo) :: refcell
  TYPE(pbc_scfiterations_t) :: pbc_it(7)
  CHARACTER(LEN=10) :: stiter
  CHARACTER(LEN=12) :: diis,diismats
  character(len=20) :: mattxt
  LOGICAL :: diis_exit
  REAL(realk) :: E_J,E_K,E_XC,E_ff,E_cell
  real(realk)         :: TS,TE,TST,TET,TOT,TWT !For finding time usage

  write(lupri,*) 'Entering routine startzdiis'

    write(stiter,'(I5)') 1
    stiter=adjustl(stiter)
    diis='diis_'//trim(stiter)//'_'

  call mem_alloc(tlat,(maxmultmom+1)**2,(maxmultmom+1)**2)

  !read the tlatticetensor to build up the multipole moments
  call read_pbc_tlattice(tlat,maxmultmom,'Tlatticetensor.dat',lupri)
 
  !config%molecule%nelectrons
  !write(*,*) 'number of electrons' molecule%nelecetrons
  



  call mem_alloc(nfdensity,numrealvec)
  !allocate(nfdensity(numrealvec))

  DO i=1,numrealvec  
     call mat_init(nfdensity(i),ndim,ndim)
     call mat_zero(nfdensity(i))
  ENDDO

  !get the fock matrices f^0l
  if(ll%store_mats)then
    call pbc_read_fock_matrix(ll,ndim,ndim,'            ') !here i initalize oper(1),oper(2)
    !get the overlap matrices S^0l
    call pbc_read_matrix(ll,ndim,ndim,1,1,'            ')
  else
    call pbc_add_fock_matrix(f_1,g_2,ll,ndim,ndim,numrealvec)
  endif

   DO k=1,BZ%nk

      call pbc_get_kpoint(k,kvec)

      call zero_pbc_elstr(Bz%fck)
      call zero_pbc_elstr(Bz%Smat)

      call pbc_rspc_to_kspc_mat(ll,Bz,ndim,kvec,2)
      if(ll%store_mats)then
        call pbc_rspc_to_kspc_mat(ll,Bz,ndim,kvec,1)
      else
        realcut(1)=ll%oneop1
        realcut(2)=ll%oneop2
        realcut(3)=ll%oneop3
        call pbc_trans_mat_to_kspc(ovl,numrealvec,ll,Bz,ndim,kvec,realcut)
      endif

      call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)
      call pbc_zdevectorize_mat(smatk,ndim,ndim,bz%smat%zelms)

      call solve_kfcsc_mat(k,ndim,fock,smatk,&
      C_k,bz%keigv((k-1)*ndim+1:k*ndim),lupri)
      
      !write(*,*) bz%keigv((k-1)*ndim+1:k*ndim)

      if(bz%kpnt(k)%is_gamma ) C_0(:,:)=C_k(:,:)
      !if(bz%kpnt(k)%is_gamma ) write(lupri,*) 'C_k'
      !if(bz%kpnt(k)%is_gamma ) call write_zmatrix(C_k,ndim,ndim,lupri)
      !if(bz%kpnt(k)%is_gamma ) write(lupri,*) 'C_0'
      !if(bz%kpnt(k)%is_gamma ) call write_zmatrix(C_0,ndim,ndim,lupri)
      
      call pbc_get_ddensity(D_k,C_k,ndim,molecule%nelectrons/2 &
          ,lupri)
      
      !write(*,*) 'Debug volbz',bz%nk,bz%nk_nosym
      call kspc_2_rspc_loop_k(nfdensity,D_k,ll,kvec,bz%kpnt(k)%weight,BZ%NK_nosym,ndim)
          


      !if(k.eq. 1) write(lupri,*) 'First C coefficient for k = 0'
      !write(lupri,*) 'First C coefficients for k'
      !if(k .eq. 1) call write_zmatrix(kdep_tmp(k)%kcdensitymat,ndim,ndim,lupri)
      !call write_zmatrix(kdep_tmp(k)%kcdensitymat,ndim,ndim,lupri)


    !  write(*,*) 'C(0)'
    !  call write_zmatrix(pbc_it(1)%kdep_it(k)%kcdensitymat,ndim,ndim)
    !  stop


      !write(lupri,*) 'eigenvalues'
      !write(lupri,*) kdep_tmp(k)%keigv
  !  endif

   ENDDO

   !do i=1,numrealvec
   !call find_latt_vectors(i,il1,il2,il3,fdim,ll)
   !if(abs(il1) .gt. ll%ndmat) CYCLE
   !if(abs(il2) .gt. ll%ndmat) CYCLE
   !if(abs(il3) .gt. ll%ndmat) CYCLE
   !write(lupri,*) 'Dmat 0l', il1,il2,il3
   !call mat_to_full(nfdensity(i),1._realk,dmat)
   !call pbc_matlab_print(dmat,ndim,ndim,'dmat',lupri)
   !write(lupri,*) 'Dmat 0l', il1,il2,il3
   !call mat_print(nfdensity(i),1,ndim,1,ndim,lupri)
   !enddo
   !write(*,*) 'density',nfdensity(13)%elms
   !call pbc_densitymat_write(nfdensity,ll,ndim,ndim,8,diis)
   call pbc_densitymat_write(nfdensity,ll,ndim,ndim,8,'            ')

   call pbc_get_onehamenergy(numrealvec,f_1,nfdensity,E_1)

   call print_bands(bz,ndim,'band-energy')
   call pbc_trans_k_energy(ll,cellenergies,nfsze,ndim,molecule%nelectrons,bz)
   write(lupri,*) 'E(HOMO) =', cellenergies(1)
   write(lupri,*) 'E(LUMO) =', cellenergies(2)
   write(*,*) 'E(HOMO) =', cellenergies(1)
   write(*,*) 'E(LUMO) =', cellenergies(2)
   !write(lupri,*) 'not exactly E(cell) =', cellenergies(3)

   call pbc_free_read_matrices(ll)

  allocate(latt_cell(numrealvec))

  call set_refcell(refcell,molecule)
  call set_lattice_cells(latt_cell,numrealvec,molecule,ll,lupri)
  errlm=molecule%nelectrons/2*(ndim-molecule%nelectrons/2.)
  !write(*,*) 'errlm= ',errlm
  call mem_alloc(error,7,errlm)

 !computing the nuclear moments
  call mem_alloc(nucmom,(1+maxmultmom)**2)
  call pbc_comp_nucmom(refcell,nucmom,maxmultmom,nfsze,lupri)
  ! self consistent iterations

  k=0
  i=0 !either 0 or 1, check it, for zero no errors
  tol=0

  Ecell=0.0_realk
  diis_exit = .false. !when diis_exit the iterations are finished
  DO WHILE(tol .lt. ll%num_its)! 20)!should have an input parameter
    CALL LSTIMER('START ',TOT,TWT,LUPRI)
    k=k+1
    i=i+1
    tol=tol+1
    !write(*,*) 'toleration of iterations: ',tol
    !We keep only data of the last 7 iterations

    if(i .ge. ll%num_store) Then !should have an input parameter to decide how man we store
      k=ll%num_store-1
      i=ll%num_store
      !DO j=1,6
      !   C_tmp(j,:,:)=c_tmp(j+1,:,:)
      !ENDDO
    endif

    !to keep order of the iterations, should maybe have one for 
    !the matrix elements
    write(stiter,'(I5)') tol+1
    stiter=adjustl(stiter)
    diis='diis_'//trim(stiter)//'_'
    write(stiter,'(I5)') tol
    stiter=adjustl(stiter)
    diismats='diis_'//trim(stiter)//'_'

    call pbc_get_onehamenergy(numrealvec,f_1,nfdensity,E_1)

    CALL LSTIMER('START ',TST,TET,LUPRI)
    CALL LSTIMER('START ',TS,TE,LUPRI)
    call pbc_electron_rep_k(lupri,luerr,setting,molecule,ndim,&
      ll,latt_cell,refcell,numrealvec,nfdensity,g_2,E_J)
    CALL LSTIMER('pbc Coul',TS,TE,LUPRI)

    !ll%compare_elmnts=.true.
    !CALL LSTIMER('START ',TS,TE,LUPRI)
    !CALL find_cutoff_twop(lupri,luerr,setting,ndim,ll,&
    !          latt_cell, refcell,numrealvec,nfdensity)
    !CALL LSTIMER('pbc_find_twop',TS,TE,LUPRI)


    CALL LSTIMER('START ',TS,TE,LUPRI)
    call pbc_exact_xc_k(lupri,luerr,setting,molecule,ndim,&
     ll,latt_cell,refcell,numrealvec,nfdensity,g_2,E_K)
    CALL LSTIMER('pbc xchange',TS,TE,LUPRI)
    CALL LSTIMER('rep xchange',TST,TET,LUPRI)

    write(lupri,*) 'nlayers exch',ll%kx1,ll%kx2,ll%kx3
    write(*,*) 'nlayers exch',ll%kx1,ll%kx2,ll%kx3

    ll%fc1=max(ll%oneop1,ll%col1)
    ll%fc1=max(ll%fc1,ll%Kx1)
    ll%fc2=max(ll%oneop2,ll%col2)
    ll%fc2=max(ll%fc2,ll%Kx2)
    ll%fc3=max(ll%oneop3,ll%col3)
    ll%fc3=max(ll%fc3,ll%Kx3)
     !ll%compare_elmnts=.false.

    CALL LSTIMER('START ',TS,TE,LUPRI)
    call pbc_fform_fck(maxmultmom,tlat,ll%tlmax,ndim,nfsze,ll,nfdensity,nucmom,&
                   g_2,E_ff,E_nn,lupri)
    CALL LSTIMER('pbc farfield',TS,TE,LUPRI)

!  if(ll%compare_elmnts) then
!    call COMPARE_MATRICES(lupri,ndim,numrealvec,nfsze,maxmultmom,ll)
!  endif

!  if(k.eq. 2) then
!    iunit=-1
!    do n1=-3,3
!      write(numtostring1,'(I5)') n1
!      numtostring1=adjustl(numtostring1)
!      mattxt='minFmatm3'//trim(numtostring1)//'00.dat'
!        !call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
!      CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
!      call find_latt_index(indexx,n1,0,0,fdim,ll,ll%max_layer)
!      write(iunit,*) ndim
!      DO j=1,ndim
!         write(iunit,*) (ll%lvec(indexx)%fck_vec(i+(j-1)*ndim),i=1,ndim)
!      ENDDO
!      call lsclose(iunit,'KEEP')
!    enddo
!    write(*,*) 'Fock matrix written to disk'
!    write(lupri,*) 'Fock matrix written to disk'
!    i=2
!  endif


  !get the fock matrices f^0l
  if(ll%store_mats) then
    CALL LSTIMER('START',TS,TE,LUPRI)
    !call pbc_read_fock_matrix(ll,ndim,ndim,diismats)
    call pbc_read_fock_matrix(ll,ndim,ndim,'            ')
    CALL LSTIMER('Reading fock',TS,TE,LUPRI)
    call pbc_fockmat_write(ll,ndim,ndim,7,2,diismats,lupri)
    !get the overlap matrices S^0l
    call pbc_read_matrix(ll,ndim,ndim,1,1,'            ')
  else
    call pbc_add_fock_matrix(f_1,g_2,ll,ndim,ndim,numrealvec)
    call pbc_fockmat_write(ll,ndim,ndim,7,2,diismats,lupri)
    realcut(1)=ll%oneop1
    realcut(2)=ll%oneop2
    realcut(3)=ll%oneop3
  endif
  !if(tol .eq. 1) call pbc_fockmat_write(ll,ndim,ndim,7,2,diis)

  DO layer=1,numrealvec  
     call mat_zero(nfdensity(layer))
  ENDDO


  CALL LSTIMER('START',TST,TET,LUPRI)
  do kpt=1,bz%nk

      call pbc_get_kpoint(kpt,kvec)

      call zero_pbc_elstr(Bz%fck)
      call zero_pbc_elstr(Bz%Smat)

      !We need k-space fock matrix
      !call pbc_rspc_to_kspc_mat(ll,Bz,ndim,kvec,2)
      !call pbc_rspc_to_kspc_mat(ll,Bz,ndim,kvec,1)

      !call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)
      !call pbc_zdevectorize_mat(smatk,ndim,ndim,bz%smat%zelms)


      if(bz%kpnt(kpt)%is_gamma )then
       
      CALL LSTIMER('START',TS,TE,LUPRI)

      !We need k-space fock matrix in gamma point now
      call pbc_rspc_to_kspc_mat(ll,Bz,ndim,kvec,2)

      call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)
        !write(*,*) pbc_it(i)%kdep_it(1)%kfockvec
        !write(*,*) 'fock(0)'
        !call write_zmatrix(fock,ndim,ndim)
        !write(lupri,*) 'fock(0)'
        !call write_zmatrix(fock,ndim,ndim,lupri)
!
        call transform_toMOfock(C_0,fock,fockMO(i,:,:),ndim,lupri)

        call mem_alloc(weight,i)
        weight=0.D0
        !! k is still the gamma point
        if(tol .gt. ll%num_store) then
          do j=1,ll%num_store-1
           error(j,:)=error(j+1,:)
          enddo
          error(ll%num_store,:)=0d0
        endif
        !get the error vectors
        call pbc_geterrorvec(error(i,:),fockMO(i,:,:),ndim,errlm,molecule%nelectrons)
        errortest=dot_product(error(i,:),error(i,:))
        errortest=sqrt(errortest)
        if(errortest .le. ll%error) diis_exit=.true.

        !Get diis weights
        call pbc_diisweights(errlm,error,weight,i,lupri)

        CALL LSTIMER('diis weights',TS,TE,LUPRI)
        if(tol .eq. 1) weight(1)=1.0d0

        write(*,*) 
        write(*,*) 'Iteration nr. ', tol
!        write(*,*) 'Weights'
!        write(*,*) weight
        write(lupri,*) 
        write(lupri,*) 'Iteration nr. ', tol
        write(lupri,*) 'Weights'
        write(lupri,*) weight
        

        !call subroutine that reads matrices writes them to disk
        !and sum with corresponding weights
        !write(*,*) 'Debug before get_weights'
        !call pbc_get_weighted_fock(i,tol,7,ndim,weight,ll)
        !write(*,*) 'Debug after get_weights'

        !We need to get gamma point fock matrix again
        !call pbc_rspc_to_kspc_mat(ll,Bz,ndim,kvec,2)

        !call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)

        ! We have to reset the gamma point again since we use it again.
        call zero_pbc_elstr(Bz%fck)
        call zero_pbc_elstr(Bz%Smat)
      endif ! is_gamma

      if(.not. diis_exit) then
        !CALL LSTIMER('START',TS,TE,LUPRI)
        call pbc_get_weighted_fock(i,tol,ll%num_store,ndim,weight,ll,lupri)
        !CALL LSTIMER('Weighted fock',TS,TE,LUPRI)
      endif

      if(ll%store_mats)then
        !We need k-space overlap
        call pbc_rspc_to_kspc_mat(ll,Bz,ndim,kvec,1)
      else
        !We need k-space overlap
        call pbc_trans_mat_to_kspc(ovl,numrealvec,ll,Bz,ndim,kvec,realcut)
      endif

      !Put overlap in a matrix form
      call pbc_zdevectorize_mat(smatk,ndim,ndim,bz%smat%zelms)

      !We need to get k point fock matrix
      call pbc_rspc_to_kspc_mat(ll,Bz,ndim,kvec,2)
      !Put it in a matrix form
      call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)

      !solves F(k)C(k)=eps(k)S(k)C(k)
      call solve_kfcsc_mat(kpt,ndim,fock,smatk,&
      C_k,bz%keigv((kpt-1)*ndim+1:kpt*ndim),lupri)

      !write(*,*) bz%keigv((k-1)*ndim+1:k*ndim)
      
      !gets D(k)
      call pbc_get_ddensity(D_k,C_k,ndim,molecule%nelectrons/2 &
          ,lupri)
      
      !Converts D(k) to D^0l
      call kspc_2_rspc_loop_k(nfdensity,D_k,ll,kvec,bz%kpnt(kpt)%weight,BZ%NK_nosym,ndim)

    if(bz%kpnt(kpt)%is_gamma ) C_0(:,:)=C_k(:,:)

    !if(bz%kpnt(kpt)%is_gamma ) then
    !  write(lupri,*) 'Testing for C coefficients in gamma point'
    !  write(lupri,*) 'Computed C coefficients'
    !  call write_zmatrix(C_k,ndim,ndim,lupri)
    !  write(lupri,*) 'Copy of C coefficients'
    !  call write_zmatrix(C_0,ndim,ndim,lupri)
    !endif

    !if(diis_exit .or. tol .eq.ll%num_its) then
    !  call pbc_cell_energy(ll,Bz,ndim,ndim,D_k,kvec,Ecell,kpt,&
    !  molecule%nelectrons/2,lupri)
    !endif

    enddo !kpt
    CALL LSTIMER('k point energy',TST,TET,LUPRI)

    if(associated(weight)) call mem_dealloc(weight)
    
    write(*,*)
    call pbc_densitymat_write(nfdensity,ll,ndim,ndim,8,'            ')
    call pbc_free_read_matrices(ll)
    !call pbc_densitymat_write(nfdensity,ll,ndim,ndim,8)
    !call pbc_trans_to_realspc(ll,nfdensity,numrealvec,ndim,bz,pbc_it(i)%kdep_it)
    call print_bands(bz,ndim,'band-energy')
    call pbc_trans_k_energy(ll,cellenergies,nfsze,ndim,molecule%nelectrons,&
    bz)

    E_cell=E_1+E_j+E_K+E_ff
    write(lupri,*) 'E(HOMO) =', cellenergies(1), tol
    write(lupri,*) 'E(LUMO) =', cellenergies(2), tol
    write(lupri,*) 'Cell Energy electrons=', E_cell
    write(lupri,*) 'K energy', E_k
    write(lupri,*) 'J energy', E_J
    write(lupri,*) 'h_1=',E_1
    write(lupri,*) 'Far field=', E_ff,E_nn
    write(*,*) 'H_1=',E_1
    write(*,*) 'E(HOMO) =', cellenergies(1)
    write(*,*) 'E(LUMO) =', cellenergies(2)
    write(*,*) 'Cell Energy =', E_cell
    write(*,*) 'K energy', E_k
    write(*,*) 'J energy', E_J
    write(*,*) 'Nuclear+h_1=',E_1
    write(*,*) 'Far field=', E_ff,E_nn

    if(diis_exit) exit

    !write(*,*) 'error(1,1)',error(1,1)
    !if(k .gt. 1) write(*,*) 'error(2,1)',error(2,1)!then
      !allocate(weight(k))
      !call pbc_diisweights(errlm,error,weight,k,lupri) 
    !endif
    !if(k.eq. 3)stop

    !do j=1,numrealvec
    !   call mat_zero(nfdensity(j))
    !enddo
    CALL LSTIMER('Diis Iteration',TOT,TWT,LUPRI)

  ENDDO

  do i=1,numrealvec
	call free_Moleculeinfo(latt_cell(i))
	call mat_free(nfdensity(i))
  enddo
  call free_Moleculeinfo(refcell)
  call mem_dealloc(tlat)
  call mem_dealloc(nfdensity)
  call mem_dealloc(nucmom)
  call mem_dealloc(error)
  deallocate(latt_cell)

  if(diis_exit) then
    write(lupri,'(A19)') 'CONVERGENCE REACHED'
    write(lupri,'(A19)') 'CONVERGENCE REACHED'
    write(lupri,'(A19)') 'CONVERGENCE REACHED'
    write(*,'(A19)') 'CONVERGENCE REACHED'
    write(*,'(A19)') 'CONVERGENCE REACHED'
    write(*,'(A19)') 'CONVERGENCE REACHED'
  else
    write(lupri,'(A28)') 'FINISHED WITHOUT CONVERGENCE'
    write(lupri,'(A28)') 'FINISHED WITHOUT CONVERGENCE'
    write(lupri,'(A28)') 'FINISHED WITHOUT CONVERGENCE'
    write(*,'(A28)') 'FINISHED WITHOUT CONVERGENCE'
    write(*,'(A28)') 'FINISHED WITHOUT CONVERGENCE'
    write(*,'(A28)') 'FINISHED WITHOUT CONVERGENCE'
  endif

  write(lupri,*) 'final E(HOMO) =', cellenergies(1)
  write(lupri,*) 'final E(LUMO) =', cellenergies(2)
  write(lupri,*) 'Cell Energy =', E_cell
  write(lupri,*) 'Nuclear+h_1=',E_1
  write(lupri,*) 'Far field=', E_ff,E_nn
  write(lupri,*) 'K energy', E_k
  write(lupri,*) 'J energy', E_J
  write(*,*) 'E(HOMO) =', cellenergies(1)
  write(*,*) 'E(LUMO) =', cellenergies(2)
  write(*,*) 'Cell Energy =', E_cell
  write(*,*) 'K energy', E_k
  write(*,*) 'J energy', E_J
  write(*,*) 'Nuclear+h_1=',E_1
  write(*,*) 'Far field=', E_ff,E_nn
  !write(*,*) 'number of electrons', molecule%nelectrons

END SUBROUTINE pbc_startzdiis

SUBROUTINE pbc_get_initial_density(kdep,bz,nrow,ncol,DMAT,nk,lupri)
IMPLICIT NONE
INTEGER,INTENT(IN) :: nrow,ncol,lupri,nk
TYPE(BZgrid_t),intent(in) :: bz
!REAL(realk),INTENT(IN) :: Smat(nrow,ncol)
REAL(realk),INTENT(INOUT) :: DMAT(nrow*ncol)
TYPE(pbc_elstr_t),INTENT(IN) :: kdep(nk)
!LOCAL variable
TYPE(pbc_elstr_t) :: kdep_tmp(nk)
INTEGER :: info!,LWORK(nrow),IPIV(nrow)
real(realk) :: lattindex(3)
REAL(realk) :: w(nrow,nk),rwork(3*nrow-2)
COMPLEX(COMPLEXK) :: work(3*nrow) 
COMPLEX(COMPLEXK) :: sk_tmp(nrow,ncol),winv(nrow,ncol,nk)
COMPLEX(COMPLEXK) :: alpha,beta
integer :: i,j

!V(:,:) =Smat(:,:)
w=0D0
if(nrow .ne. ncol) then
  call LSQUIT('pbc_get_initial_density: dmat not n X n',lupri)
endif
write(lupri,*) 'number of k points',bz%nk
DO i=1,nk
call init_pbc_elstr(kdep_tmp(i),nrow,ncol)
kdep_tmp(i)%koverlapmat=kdep(i)%koverlapmat

!call write_matrix(DMAT,nrow,ncol)
!write(*,*)
!call write_matrix(SMAT,nrow,ncol)
!call DGETRI(nrow,V,nrow,IPIV,work,3*nrow,info)
!if(info .ne. 0) then
!call LSQUIT('pbc_get_initial_density: info not 0 for inversion',lupri)
!endif
!write(*,*) 'inverse of S 0 0 0'
!call write_matrix(V,nrow,ncol)
!call write_matrix(DMAT,nrow,ncol)
!call ZHETRI(nrow,skmat(:,:,i),nrow,IPIV,work,3*nrow,info)
call ZHEEV('V','U',nrow,kdep_tmp(i)%koverlapmat,ncol,w(:,i),work,3*nrow,rwork,info)

if(info .ne. 0) then
  call LSQUIT('pbc_get_initial_density: info not 0 for eigenvalues',lupri)
endif

!winv(:,:)=0d0
do j=1,nrow
winv(j,j,i)=CMPLX(1./sqrt(w(j,i)),0D0,complexk)
enddo
!write(*,*) w(:,i)
alpha=CMPLX(1D0,0d0,complexk)
beta=CMPLX(0D0,0D0,complexk)
call zgemm('N','N',nrow,ncol,ncol,alpha,winv,nrow,kdep_tmp(i)%koverlapmat,&
                                                     nrow,beta,sk_tmp,nrow)

call zgemm('C','N',nrow,ncol,ncol,alpha,kdep_tmp(i)%koverlapmat,&
                        nrow,sk_tmp,nrow,beta,kdep_tmp(i)%kddensitymat,nrow)

!call zgemm('T','N',nrow,ncol,ncol,alpha,V_tmp,nrow,V_tmp,nrow,beta,DMAT,nrow)
!write(*,*) kdep_tmp(i)%kddensitymat
ENDDO

!write(*,*) 'nfsize in initial',nfsze
!call pbc_trans_to_realspc(ll,nfdensity,nfsze,ndim,bz,kdep_tmp)
lattindex(:)=0
write(*,*)
call transformk_2_realmat(kdep_tmp,bz,&
          DMAT,nrow,lattindex,lupri)

!write(*,*) 'D(00)'
!write(*,*) DMAT
!CMAT(:,:)=sqrt(DMAT(:,:))
!write(*,*) 'inverse of S 0 0 0,info=',info,'w =',w
!
!call write_matrix(V,nrow,ncol)
!write(*,*)'sqrt(inv(S))'
!!call write_matrix(CMAT,nrow,ncol)
!alpha=1D0
!beta=0D0
!call dgemm('N','N',nrow,ncol,ncol,alpha,winv,nrow,V,nrow,beta,V_tmp,nrow)
!call dgemm('T','N',nrow,ncol,ncol,alpha,V,nrow,V_tmp,nrow,beta,V_tmp,nrow)
!call dgemm('T','N',nrow,ncol,ncol,alpha,V_tmp,nrow,V_tmp,nrow,beta,DMAT,nrow)
!call write_matrix(DMAT,nrow,ncol)

END SUBROUTINE pbc_get_initial_density

!SUBROUTINE to the the electronic one-body energy
SUBROUTINE  pbc_get_onehamenergy(numvecs,f_1,nfdensity,E_1)
INTEGER,INTENT(IN) :: numvecs
TYPE(matrix),target,intent(IN) :: f_1(numvecs),nfdensity(numvecs)
REAL(REALK),INTENT(INOUT) ::  E_1
!LOCAL
INTEGER :: celli

 E_1=0
 do celli=1,numvecs
    if(f_1(celli)%init_magic_tag.EQ.mat_init_magic_value) THEN
      E_1=E_1+mat_dotproduct(f_1(celli),nfdensity(celli))
    endif
 enddo
 

END SUBROUTINE  pbc_get_onehamenergy


SUBROUTINE pbc_get_ddensity(ddensity,C_tmp,nbast,nkmobas,lupri)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nbast,nkmobas,lupri
  !TYPE(matrix),INTENT(INOUT) :: nfdensity(nfsize)
  COMPLEX(COMPLEXK),intent(INOUT) :: ddensity(nbast,nbast)
  COMPLEX(COMPLEXK),intent(in) :: C_tmp(nbast,nbast)
  !real(realk),intent(INOUT) :: ddensity(nbast,nbast)
  !real(realk),intent(in) :: C_tmp(nbast,nbast)
  !TYPE(lvec_data_t) :: C_tmp(nfsize)
  !LOCAL VARIABLES
  COMPLEX(COMPLEXK), pointer :: density_tmp(:,:)
  INTEGER :: i,j,mu,nu
  COMPLEX(COMPLEXK) :: alpha,beta

  alpha=CMPLX(2D0,0d0,complexk)
  beta =CMPLX(0d0,0d0,complexk)
  
  i=1
  j=1
      !write(lupri,*) 'C coefficients in density comp'
      !call write_zmatrix(C_tmp,nbast,nbast,lupri)

!     allocate(density_tmp(nbast,nbast))
     call mem_alloc(density_tmp,nbast,nbast)
     density_tmp=0d0
     !call zgemm('n','c',nbast,nbast,nbast,alpha,c_tmp,nbast,c_tmp,&
     !           nbast,beta,density_tmp,nbast)
     !write(lupri,*) 'c*c'

     DO i=1,nkmobas
      DO mu=1,nbast
       DO nu=1,nbast
         density_tmp(mu,nu)=density_tmp(mu,nu)+2D0*c_tmp(mu,i)*conjg(c_tmp(nu,i))
        ! write(lupri,*)c_tmp(mu,i)*conjg(c_tmp(nu,i))
       ENDDO
      ENDDO
     ENDDO


      DO i=1,nbast
        DO j=1,nbast
           ddensity(i,j)=density_tmp(i,j)
        ENDDO    
      ENDDO  
      !write(lupri,*) 'dk'
      !call write_zmatrix(ddensity,nbast,nbast,lupri)
      !deallocate(density_tmp)
      call mem_dealloc(density_tmp)




END SUBROUTINE pbc_get_ddensity



SUBROUTINE pbc_diisweights(errdim,error,weight,it,lupri)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: errdim,lupri,it
  REAL(realk),INTENT(INOUT) :: error(7,errdim)
  REAL(realk),INTENT(INOUT) :: weight(it)
  !LOCAL VARIABLES
  !REAL(realk) :: B_mat(it+1,it+1),weight_tmp(it+1,1)
  REAL(realk) :: B_mat(it,it),weight_tmp(it,1)
  REAL(realk) :: Sdgelss(it+1),rcond,normfac
  INTEGER :: i,j,info, N,m,rank,lwork
!  TYPE(matrix) :: Bmat_t
  !INTEGER :: solve(it+1)
  INTEGER :: solve(it)
  REAL(realk),pointer :: work(:)
  info=0

  B_mat(:,:)=0d0
  DO i=1,it
    write(lupri,*) 'error(',i,')'
    write(lupri,*) error(i,:)
    !B_mat(i,i)= 0.05D0
   DO j=1,it
    B_mat(i,j)=B_mat(i,j)+dot_product(error(i,:),error(j,:))! error(i)*error(j)
   ENDDO
  ENDDO
!  call mat_init(Bmat_t,it,it)
!  call mat_zero(Bmat_t)
!  call mat_set_from_full(B_mat,1D0,Bmat_t)

!  B_mat(it+1,it+1)=0D0
  solve=0
  !weight_tmp(it+1,1)=-1D0
  weight_tmp(:,1)=1!D0

  !N=it+1
  N=it
  m=1
!  write(*,*) 'Debug b_mat'
!  call write_matrix(B_mat,it+1,it+1)
  !write(lupri,*) 'Diis B matrix'
  !call write_matrix(B_mat,it+1,it+1,lupri)
  !call write_matrix(B_mat,it,it,lupri)
!  write(*,*) 'DEBUG 2: get weights '
!  call dposv('U',N,m,B_mat,N,weight_tmp,N,info)
!  info=1
!  IF(info .ne. 0) THEN
!    write(*,*) 'Calls dgesv instead'
    call dgesv(N,m,B_mat,N,solve,weight_tmp,N,info)
    IF(info .ne. 0) THEN
!      write(lupri,*) 'ERROR ERROR, diis matrix below'
!      write(*,'(X,A21,X,I3.2)') 'INFO NOT ZERO, INFO = ',info
!      write(*,*) 'if INFO = -i, the i-th argument had an illegal value'
!      write(*,*) 'if INFO = i, U(i,i) is exactly zero. The factorization has been'
!      write(*,*) 'completed, but the factor U is exactly singular,'
!      write(*,*) 'so the solution  could not  be computed'
!      write(*,*) 
!      write(*,*) 'Iteration number ', it
      write(*,*) 'Calls dgelss instead'
      write(lupri,*) 'Calls dgelss instead'
      DO i=1,it
       DO j=1,it
        B_mat(i,j)=dot_product(error(i,:),error(j,:))! error(i)*error(j)
        B_mat(j,i)=B_mat(i,j)
       ENDDO
      ENDDO
      DO j=1,it
      ! B_mat(it+1,j)=-1.D0
      ! B_mat(j,it+1)=-1D0
       !weight_tmp(j,1)=0d0
       weight_tmp(j,1)=1d0
      ENDDO
      !B_mat(it+1,it+1)=0D0
      !weight_tmp(it+1,1)=-1D0
      lwork=3*n+2*n+1000
      rcond= 1.0D-10
      call mem_alloc(work,lwork)
      call dgelss(n,n,m,B_mat,n,weight_tmp,n,Sdgelss,rcond,rank,work,lwork,info)
      !call dgels('N',n,n,1,B_mat,n,weight_tmp,n,work,lwork,info)
      !write(*,*) 'lwork',lwork,sdgelss(1),sdgelss(1)*rcond
      if(info .ne. 0) then
          call LSQUIT('pbc_diisweights: INFO not zero, while solving eigenvalue',&
                      lupri)
      endif
      call mem_dealloc(work)
    endif
!  ENDIF

!  write(*,*) 'DEBUG 3: get weights '
  DO i=1,it
   weight(i)=weight_tmp(i,1)
!   write(*,*) weight(i)
  ENDDO
  
  !Normalizing the weights
  normfac=sum(weight(1:it))
  !write(lupri,*) 'norm',normfac,weight
  weight(1:it)=weight(1:it)/normfac

!  call mat_free(Bmat_t)

END SUBROUTINE pbc_diisweights

!Sums over former fock matrices
!with the corresponing weights.
SUBROUTINE pbc_get_weighted_fock(i,ndiis,stdiis,nbast,weight,Aop,lupri)
  IMPLICIT NONE
  !DUMMY variables
  INTEGER,INTENT(IN)              :: i,ndiis,stdiis,nbast,lupri
  REAL(realk),INTENT(IN)          :: weight(i)
  TYPE(lvec_list_t),INTENT(INOUT) :: Aop
               !LOCAL VARIABLES
  INTEGER                         :: j,layer,k
  INTEGER                         :: l1,l2,l3
  CHARACTER(len=12)               :: diis,tmpdiis
  CHARACTER(len=7)                :: stiter
  TYPE(lvec_list_t)               :: tmp_mat


  tmp_mat%max_layer=Aop%max_layer
  tmp_mat%nneighbour=Aop%nneighbour
  tmp_mat%ldef%is_active(:)=Aop%ldef%is_active(:)
  call build_lvec_list(tmp_mat,nbast) 
  tmp_mat%fc1=Aop%fc1
  tmp_mat%fc2=Aop%fc2
  tmp_mat%fc3=Aop%fc3
  tmp_mat%oneop1 = Aop%oneop1
  tmp_mat%oneop2 = Aop%oneop2
  tmp_mat%oneop3 = Aop%oneop3
  tmp_mat%col1 =  Aop%col1
  tmp_mat%col2 =  Aop%col2
  tmp_mat%col3 =  Aop%col3
  tmp_mat%Kx1 =  Aop%Kx1
  tmp_mat%Kx2 =  Aop%Kx2
  tmp_mat%Kx3 =  Aop%Kx3
  !call mat_init(tmp_mat,nbast,nbast)
  !call mat_zero(tmp_mat)

  DO layer=1,size(Aop%lvec)
     call mat_init(tmp_mat%lvec(layer)%oper(2),nbast,nbast)
     call mat_zero(tmp_mat%lvec(layer)%oper(2))
     tmp_mat%lvec(layer)%g2_computed=Aop%lvec(layer)%g2_computed
     tmp_mat%lvec(layer)%f1_computed=Aop%lvec(layer)%f1_computed
     tmp_mat%lvec(layer)%ovl_computed=Aop%lvec(layer)%ovl_computed
  ENDDO

  if(ndiis .ge. 1) then
    write(stiter,'(I5)') ndiis
    stiter=adjustl(stiter)
    !write(*,*) 'stiter  ',stiter
    diis='diis_'//trim(stiter)//'_'
  endif

  Do layer=1,size(Aop%lvec)
     l1=int(Aop%lvec(layer)%lat_coord(1))
     l2=int(Aop%lvec(layer)%lat_coord(2))
     l3=int(Aop%lvec(layer)%lat_coord(3))
     if((abs(l1) .le. Aop%fc1 .and. abs(l2) .le. Aop%fc2) .and. abs(l3) .le. Aop%fc3)then
       call mat_zero(Aop%lvec(layer)%oper(2))
     endif
  ENDDO

  if(ndiis .le. stdiis) then
    do j=1,i
    write(stiter,'(I5)') j
    stiter=adjustl(stiter)
    tmpdiis='diis_'//trim(stiter)//'_'

    call pbc_read_matrix(tmp_mat,nbast,nbast,7,2,tmpdiis)

    Do layer=1,size(Aop%lvec)
     l1=int(Aop%lvec(layer)%lat_coord(1))
     l2=int(Aop%lvec(layer)%lat_coord(2))
     l3=int(Aop%lvec(layer)%lat_coord(3))
     if((abs(l1) .le. Aop%fc1 .and. abs(l2) .le. Aop%fc2) .and. abs(l3) .le. Aop%fc3)then
       call mat_daxpy(weight(j),tmp_mat%lvec(layer)%oper(2),&
       Aop%lvec(layer)%oper(2))
     endif
    enddo!layer
    enddo!j
  else
    k=0
    do j= ndiis-stdiis+1,ndiis
    k=k+1
    write(stiter,'(I5)') j
    stiter=adjustl(stiter)
    tmpdiis='diis_'//trim(stiter)//'_'
    !write(*,*) 'Debug 3 inside get_weights'
    write(lupri,*) 'Filename for fock matrix, ', tmpdiis

    call pbc_read_matrix(tmp_mat,nbast,nbast,7,2,tmpdiis)
    !write(*,*) 'Debug 4 inside get_weights'

    Do layer=1,size(Aop%lvec)
     l1=int(Aop%lvec(layer)%lat_coord(1))
     l2=int(Aop%lvec(layer)%lat_coord(2))
     l3=int(Aop%lvec(layer)%lat_coord(3))
     if((abs(l1) .le. Aop%fc1 .and. abs(l2) .le. Aop%fc2) .and. abs(l3) .le. Aop%fc3)then
       call mat_daxpy(weight(k),tmp_mat%lvec(layer)%oper(2),&
       Aop%lvec(layer)%oper(2))
     endif
    enddo!layer
    enddo!j
  endif

  Do layer=1,size(Aop%lvec)
   l1=int(Aop%lvec(layer)%lat_coord(1))
   l2=int(Aop%lvec(layer)%lat_coord(2))
   l3=int(Aop%lvec(layer)%lat_coord(3))
   if((abs(l1) .le. Aop%fc1 .and. abs(l2) .le. Aop%fc2) .and. abs(l3) .le. Aop%fc3)then
!!   ! write(*,*) 'Debug 5 inside get_weights',l1
     call pbc_get_file_and_write(Aop,nbast,nbast,layer,7,2,'            ')
!     call pbc_get_file_and_write(Aop,nbast,nbast,layer,7,2,diis)
   endif
   call mat_free(tmp_mat%lvec(layer)%oper(2))
  Enddo
   !  call pbc_fockmat_write(Aop,nbast,nbast,7,2)
   ! write(*,*) 'Debug 6 inside get_weights'

  deallocate(tmp_mat%lvec)
END SUBROUTINE pbc_get_weighted_fock


SUBROUTINE pbc_trans_k_energy(ll,cenergies,nvecsrs,nbast,nelectrons,bz)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvecsrs,nbast,nelectrons
  TYPE(lvec_list_t),intent(IN) :: ll
  TYPE(BZgrid_t),intent(in) :: bz
!  TYPE(pbc_elstr_t),INTENT(IN) :: kdep_tmp(bz%nk)
  REAL(realk),intent(INOUT) :: cenergies(nbast)
  !LOCAL VARIABLElS
  INTEGER :: lattindex(3)!,irealspc
  INTEGER :: i,kpt
  REAL(Realk) :: ehomo,elumo,etmph1,etmph2,etmpl1,etmpl2

  lattindex(1)=0
  cenergies =0.d0
  !write(*,*) ll%ldef%is_active(1)
  !write(*,*) ll%ldef%is_active(2)
  !write(*,*) ll%ldef%is_active(3)

     lattindex(1)=0
     lattindex(2)=0
     lattindex(3)=0
     
   
    DO kpt=2,bz%nk
       etmph1=bz%keigv((kpt-2)*nbast+nelectrons/2)
       etmph2=bz%keigv((kpt-1)*nbast+nelectrons/2)
       etmpl1=bz%keigv((kpt-2)*nbast+nelectrons/2+1)
       etmpl2=bz%keigv((kpt-1)*nbast+nelectrons/2+1)
       !etmph1=kdep_tmp(kpt-1)%keigv(nelectrons/2)
       !etmph2=kdep_tmp(kpt)%keigv(nelectrons/2)
       !etmpl1=kdep_tmp(kpt-1)%keigv(nelectrons/2+1)
       !etmpl2=kdep_tmp(kpt)%keigv(nelectrons/2+1)
       ehomo=max(etmph1,etmph2)
       elumo=min(etmpl1,etmpl2)
    ENDDO
       cenergies(1)=ehomo
       cenergies(2)=elumo

!     DO kpt=1,bz%nk
!      do i=1,nelectrons/2
!       cenergies(3)=cenergies(3)+2d0*kdep_tmp(kpt)%keigv(i)*&
!       &bz%kpnt(kpt)%weight/BZ%NK_nosym
!    enddo
!   enddo
!!       cenergies(3)=behomo
!       cenergies(4)=belumo


END SUBROUTINE  pbc_trans_k_energy

SUBROUTINE pbc_trans_to_realspc(ll,nfdensity,nvecsrs,nbast,bz,kdep_tmp,lupri)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvecsrs,nbast,lupri
  TYPE(lvec_list_t),intent(IN) :: ll
  TYPE(BZgrid_t),intent(in) :: bz
  TYPE(pbc_elstr_t),INTENT(IN) :: kdep_tmp(bz%nk)
  TYPE(matrix),intent(INOUT) :: nfdensity(nvecsrs)
  !LOCAL VARIABLElS
  real(realk) :: dkmat(nbast,nbast)
  INTEGER :: irealspc,idum
  INTEGER :: l1,l2,l3,fdim(3)!,i,j
  character(len=20) :: mattxt,numtostring1,numtostring2
  INTEGER, SAVE :: dmt=1
  dmt=dmt+1

  !write(*,*) ll%ldef%is_active(1)
  !write(*,*) ll%ldef%is_active(2)
  !write(*,*) ll%ldef%is_active(3)
  DO irealspc=1,nvecsrs

     call find_latt_vectors(irealspc,l1,l2,l3,fdim,ll)
     call mat_zero(nfdensity(irealspc))
     if(abs(l1) .gt. ll%ndmat) CYCLE
     if(abs(l2) .gt. ll%ndmat) CYCLE
     if(abs(l3) .gt. ll%ndmat) CYCLE
     call transformk_2_realmat(kdep_tmp,bz,&
          dkmat,nbast,ll%lvec(irealspc)%std_coord,lupri)
     
     call mat_set_from_full(dkmat,1D0,nfdensity(irealspc))
     
       !write(lupri,*) 'Density matrix after first solution, for l1=',l1
       !call mat_print(rsdensity(irealspc),nbast,nbast,nbast,nbast,lupri)
       !write(*,*) 'density realspace'
       !do j=1,nbast
       ! write(lupri,*) (rsdensity(irealspc)%elms(j+(i-1)*nbast),i=1,nbast)
       write(numtostring1,'(I5)')  l1
       write(numtostring2,'(I5)')  dmt
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       mattxt='mindmt'//trim(numtostring2)//trim(numtostring1)//'.dat'
       !enddo
       idum=-1
       call LSOPEN(idum,mattxt,'UNKNOWN','FORMATTED')
       write(idum,*) 4
       call write_tmatrix(dkmat,nbast,nbast,idum)
       call LSCLOSE(idum,'KEEP')

  ENDDO


END SUBROUTINE  pbc_trans_to_realspc


SUBROUTINE pbc_startddiis()
  IMPLICIT NONE
END SUBROUTINE pbc_startddiis

SUBROUTINE pbc_finish_zdiis()
  IMPLICIT NONE
END SUBROUTINE pbc_finish_zdiis

SUBROUTINE pbc_finsih_ddiis()
  IMPLICIT NONE
END SUBROUTINE pbc_finsih_ddiis

! The LAPACK routine ZGGEV gives back gen. eig. vectors normalized in an
! odd way. This subroutine simply fixes the normalization.
subroutine pbc_fixzggevnorm(siz,cocoeff,metric,alphavec,betavec,lupri)
  implicit none
  ! input and output arguments
  integer, intent(in) :: siz,lupri
  complex(COMPLEXK), intent(inout) :: cocoeff(siz,siz)
  complex(COMPLEXK), intent(in) :: metric(siz,siz)
  complex(COMPLEXK), intent(inout) :: alphavec(siz), betavec(siz)
  ! local variables
  integer :: ivec, i, j, mu, nu, b
  real(realk) :: normfac, coeff_sum, coeff_max
  complex(COMPLEXK) :: dmat_elem, n_occ

  do ivec = 1,siz
     normfac = 0.0D0
     coeff_sum = 0.0D0
     coeff_max = 0.0D0
     do i = 1,siz
        call lsquit('FIXME: coeff_max real while cocoeff(i,ivec) is complex',-1)
!        maybe but abs(cocoeff(i,ivec)) into a real(realk) before calling max
!        coeff_max = max(coeff_max, abs(cocoeff(i,ivec)))
        do j = 1,siz
           normfac = normfac + real(conjg(cocoeff(i,ivec)) &
                &  * metric(i,j) * cocoeff(j,ivec))
           coeff_sum = coeff_sum + real(conjg(cocoeff(i,ivec)) &
                &    * cocoeff(j,ivec))
        end do
     end do
     !write(lupri,*) 'normfac, coff_sum',normfac,coeff_sum

#if 1
     ! crude fix for linear dependence...
     if (normfac .lt. 1.0D-6) then
        write(LUPRI,*) 'Warning: Zero norm orbital (',ivec,normfac,').'
        !write(LUPRI,*) 'Action: Avoid occupation by setting orb. energy high.'
        !alphavec(ivec) = 1.0D4
        !betavec(ivec) = 1.0D0
        normfac = 1.0D50
        !call lsquit('Zero norm vector - eeeeeh',-1)
     end if
#endif

     cocoeff(:,ivec) = cocoeff(:,ivec) / sqrt(normfac)
  end do

!Careful
  n_occ = ( 0.0D0, 0.0D0 )
  loop_ao1: do mu = 1,siz
     loop_ao2: do nu = 1,siz
        dmat_elem = ( 0.0D0, 0.0D0 )
        loop_band: do b = 1,siz
           dmat_elem = dmat_elem + cocoeff(mu,b) * conjg(cocoeff(nu,b))
        end do loop_band
        n_occ = n_occ + dmat_elem * conjg(metric(mu,nu))
     end do loop_ao2
  end do loop_ao1
  write(LUPRI,*) ' k-space occupation, n_occ(k) = ',n_occ
!Careful
end subroutine pbc_fixzggevnorm



END MODULE pbc_scfdiis
#endif

