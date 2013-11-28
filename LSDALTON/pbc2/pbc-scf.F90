#ifdef MOD_UNRELEASED
MODULE pbc_scfdiis
	USE TYPEDEF
	USE precision
	USE matrix_module
	USE lattice_vectors
	USE lattice_type
	!  USE multipole_pbc
	!  USE harmonics_pbc
	USE pbc_matrix_operations
	USE pbc_interactions
	USE pbcffdata
	USE PBC_MSC
	USE PBC_kspce_rspc_operations
	USE pbc_ff_contrib

	PRIVATE
	PUBLIC :: pbc_startzdiis

CONTAINS

! Todo not in use. (delete??) 
SUBROUTINE pbc_zeigsolve(A,B,N,M,eigv,lupri)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: N,M,lupri
	COMPLEX(complexk), INTENT(INOUT) :: A(N,M),B(N,M)
	REAL(realk),INTENT(INOUT) :: eigv(N)
	! local
	COMPLEX(complexk), POINTER :: work(:)
	REAL(realk),POINTER :: rwork(:)
	INTEGER :: info,lwork,lrwork,iwork,liwork
	
	lwork=2*n-1
	lwork=2*n+n*n
	lrwork=1+5*n+2*n*n
	liwork=3+5*n

	call mem_alloc(work,lwork)
	call mem_alloc(rwork,lrwork)
	call zhegvd(1,'V','U',n,A,n,B,n,eigv,work,lwork,rwork,lrwork,iwork,liwork,info)

	if(info .ne. 0) THEN
		write(lupri,*) 'ERROR: zhegv problems, info=', info
		call write_zmatrix(B,n,n)
		write(*,*) 'ERROR: zhegv problems, info=', info
		call LSQUIT('pbc_zeigsolve: INFO not zero, while solving eigenvalue',lupri)
	endif

	call mem_dealloc(work)
	call mem_dealloc(rwork)

END SUBROUTINE pbc_zeigsolve


!Todo not in use. (delete?)
!> \author 
!> \date
!> \brief
!> \param
!> \param
SUBROUTINE pbc_zggeigsolve(kindex,a,b,smatk,n,m,eigv,lupri)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: kindex,n,m,lupri
	COMPLEX(complexk), INTENT(INOUT) :: a(n,m),b(n,m),smatk(n,m)
	REAL(realk),INTENT(INOUT) :: eigv(N)
	! local
	COMPLEX(complexk), POINTER :: work(:)
	REAL(realk),POINTER :: rwork(:)
	COMPLEX(complexk) :: alphavec(n),betavec(n)
	COMPLEX(complexk) :: vl(n,m),ccoeff(n,m)
	COMPLEX(complexk) :: ztemp
	INTEGER :: info,iband,i,j,lwork
	INTEGER,SAVE :: ncalls=0

	ncalls=ncalls+1
	lwork=2*n*n

	call mem_alloc(work,lwork)
	call mem_alloc(rwork,8*N)

	!write(*,*) 'debug after read tlatt',kindex
	call zggev('N','V',n,A,n,B,n,alphavec,betavec,vl,1,ccoeff,n,work,2*n,rwork,info)
	!write(*,*) 'debug after read tlatt',kindex

	if(info .ne. 0) THEN
		write(lupri,*) 'ERROR: zggev problems, info=', info
		write(*,*) 'ERROR: zggev problems, info=', info
		write(*,*) 'ERROR: zggev problems, for kindex = ', kindex
		write(lupri,*) 'ERROR: zggev problems, for kindex = ', kindex
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
		ccoeff(:,iband) = 0.0_realk
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

!> \author JR 
!> \date 2013
!> \brief For a matrix pos. semdefinite matrix S : calculates U = V\sigma^{-1/2} 
!> \brief where V\sigmaV^H = S. Removes singularities that are smaller than 
!> \brief singular_threshh.
!> \param s_input 				ndim x ndim complex matrix
!> \param umat 					Changed to V\sigma
!> \param is_singular 			Set to true if singularities found
!> \param ndim 					Matrix dim
!> \param nsingular 				ndim - rank(Sabk)
!> \param singular_threshh 	Singularity threshhold. 
!> \param lupri 					Logical print unit
SUBROUTINE pbc_spectral_decomp_ovl(s_input,umat,is_singular,ndim,nsingular,& 
		& singular_threshh,lupri)
	IMPLICIT NONE
	! input
	INTEGER,INTENT(IN) :: ndim,lupri
	INTEGER,INTENT(INOUT) :: nsingular
	COMPLEX(complexk),INTENT(IN) :: s_input(Ndim,Ndim)
	COMPLEX(complexk),INTENT(INOUT),pointer :: umat(:,:)
	REAL(realk), INTENT(IN) :: singular_threshh
	LOGICAL,INTENT(INOUT) :: is_singular
	! local
	REAL(realk),POINTER :: rwork(:), eigv(:)
	REAL(realk) :: rtemp
	COMPLEX(complexk) ,POINTER :: work(:), vmat(:,:),eigv_inv(:),ctmp(:)
	COMPLEX(complexk) :: alpha, beta
	INTEGER :: info,i,j,lwork,nonsingdim

	nsingular = 0
	lwork = 10*ndim-1
	is_singular = .false.

	call mem_alloc(work,max(1,lwork))
	call mem_alloc(rwork,max(1,3*ndim-2))
	call mem_alloc(eigv,ndim)
	call mem_alloc(vmat,ndim,ndim)
	call mem_alloc(eigv_inv,ndim)
	call mem_alloc(ctmp,ndim)

	alpha=cmplx(1._realk,0._realk,complexk)
	beta=cmplx(0._realk,0._realk,complexk)

	vmat(:,:)=s_input(:,:)
	eigv_inv(:)=cmplx(0._realk,0._realk,complexk)
	! calculate spectral decomposition s = v^h\sigma v  
	call zheev('v','u',ndim,vmat,ndim,eigv,work,lwork,rwork,info)

	do i=1,ndim
		if(eigv(i) .lt. 0._realk) then
			write(*,*) 'Error: Eigenvalue of overlap matrix is negative. &
				& Eigv(',i,') =', eigv(i)
			write(lupri,*) 'Error: Eigenvalue of overlap matrix is negative. &
				& Eigv(',i,') =', eigv(i)
		endif
	end do

	! reorder elements s.t. largest eigenv comes first
	do i=1,(ndim+1)/2 
		j=ndim-i+1
		! swap eigenvectors
		ctmp(:)=vmat(:,i)
		vmat(:,i)=vmat(:,j)
		vmat(:,j)=ctmp(:)
		! swap eigenvalues
		rtemp=eigv(i)
		eigv(i)=eigv(j)
		eigv(j)=rtemp
	end do

	! remove linear dependencies in v (vmat)
	do i=1,ndim
		if(eigv(i) .lt. singular_threshh) then
			nsingular=nsingular+1
			vmat(:,i)=cmplx(0._realk,0._realk,complexk) 
			eigv_inv(i)=cmplx(0._realk,0._realk,complexk)
			is_singular=.true.
		else
			eigv_inv(i)=cmplx(1._realk/sqrt(eigv(i)),0._realk,complexk)
		endif
	enddo

#if debugpbc
	write(*,*) 'number of singulars',nsingular
	write(lupri,*) 'number of singulars',nsingular
#endif
	! calculate u = v\sigma
	nonsingdim=ndim-nsingular
	call mem_alloc(umat,ndim,nonsingdim)
	do i = 1, nonsingdim
		umat(:,i) = vmat(:,i)*eigv_inv(i)
	end do

	call mem_dealloc(work)
	call mem_dealloc(rwork)
	call mem_dealloc(eigv)
	call mem_dealloc(vmat)
	call mem_dealloc(eigv_inv)

END SUBROUTINE pbc_spectral_decomp_ovl

!> \author JR
!> \date 2013
!> \brief Do the unitary transformation Ft = Uk^H F U 
!> \param ndim 		matrix dim
!> \param nosingdim 	matrix dim
!> \param kfock 		to be transformed (ndim x ndim)
!> \param Uk 			unitary mat (ndim x nonsingdim)
!> \param tfock 		output transformed mat (nonsingdim x nonsingdim)
SUBROUTINE pbc_unitary_transform(Ndim,nonsingdim,kfock,Uk,tfock)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,nonsingdim
  COMPLEX(complexk),INTENT(IN) :: kfock(ndim,ndim),Uk(ndim,nonsingdim)
  COMPLEX(complexk),INTENT(INOUT) :: tfock(Nonsingdim,Nonsingdim)
  ! local
  COMPLEX(complexk),pointer :: tmp(:,:)
  COMPLEX(complexk) :: alpha,beta

  alpha=CMPLX(1._realk,0._realk,complexk)
  beta=CMPLX(0._realk,0._realk,complexk)
  call mem_alloc(tmp,ndim,nonsingdim)
  call zgemm('N','N',Ndim,Nonsingdim,Ndim,alpha,kfock,Ndim, &
             & Uk,Ndim,beta,tmp,Ndim)
  call zgemm('C','N',Nonsingdim,Nonsingdim,Ndim,alpha,Uk,Ndim, &
             & tmp,Ndim,beta,tfock,Nonsingdim)
  call mem_dealloc(tmp)

END SUBROUTINE pbc_unitary_transform


!todo not in use (delete?)
!> \author JR
!> \date
!> \brief
!> \param
!> \param
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


!todo not in use (delete?)
!> \author 
!> \date
!> \brief
!> \param
!> \param
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
    D(i,j) =cmplx(0._realk,0._realk,complexk)
   ENDDO
  ENDDO

  alpha=cmplx(1._realk,0._realk,complexk)
  beta=cmplx(0._realk,0._realk,complexk)
  call zgemm('C','N',ndim,ndim,ndim,alpha,U,ndim,A,ndim,&
        beta,D,ndim)
  call zgemm('N','N',ndim,ndim,ndim,alpha,D,ndim,U,ndim,&
        beta,A,ndim)

END SUBROUTINE pbc_zdiagonalize


!todo not in use (delete?)
!> \author 
!> \date
!> \brief
!> \param
!> \param
SUBROUTINE ztransformbackC(C,ndim,U)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: ndim
	COMPLEX(COMPLEXK),INTENT(INOUT) :: C(ndim,ndim)
	COMPLEX(COMPLEXK),INTENT(IN) :: U(ndim,ndim)
	COMPLEX(COMPLEXK) :: ctmp(ndim,ndim)
	COMPLEX(COMPLEXK) :: alpha,beta

	alpha=cmplx(1._realk,0._realk,complexk)
	beta =cmplx(0._realk,0._realk,complexk)
	ctmp(:,:)=C(:,:)
	call zgemm('N','N',ndim,ndim,ndim,alpha,U,ndim,ctmp,ndim, &
		& beta,C,ndim)

END SUBROUTINE ztransformbackC

!> \author JR
!> \date 2013
!> \brief
!> \param
!> \param
SUBROUTINE transform_toMOfock(C_tmp,fock,fockMO,ndim,lupri)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: ndim,lupri
	COMPLEX(complexk),INTENT(IN) :: C_tmp(ndim,ndim),fock(ndim,ndim)  
	COMPLEX(complexk),INTENT(OUT) :: fockMO(ndim,ndim)
	! local
	COMPLEX(complexk), POINTER :: fockMOtmp(:,:)
	COMPLEX(complexk) :: alpha,beta

	call mem_alloc(fockMOtmp,ndim,ndim)

	alpha=CMPLX(1._realk,0._realk,complexk)
	beta =CMPLX(0._realk,0._realk,complexk)

	!fockMO(:,:)=CMPLX(0._realk,0._realk,complexk) !not necc
	!fockMOtmp(:,:)=CMPLX(0._realk,0._realk,complexk) !not necc

	call zgemm('C','N',ndim,ndim,ndim,alpha,C_tmp,ndim,fock,ndim, &
		& beta,fockMOtmp,ndim)
	call zgemm('N','N',ndim,ndim,ndim,alpha,fockMOtmp,ndim,C_tmp,ndim, &
		& beta,fockMO,ndim)

	call mem_dealloc(fockMOtmp)

END SUBROUTINE transform_toMOfock

!> \author 
!> \date
!> \brief
!> \param
!> \param
SUBROUTINE  pbc_geterrorvec(error,fockMO,ndim,errlm,nelectrons)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: ndim,errlm,nelectrons
	REAL(realk),INTENT(INOUT) :: error(errlm)
	COMPLEX(complexk),INTENT(IN) :: fockMO(ndim,ndim)
	! local
	INTEGER :: i,j,k

	k=0

	DO i=1,nelectrons/2
		DO j=nelectrons/2+1,ndim
			k=k+1
			error(k)=abs(real(fockMO(i,j),realk))
			error(k)=error(k)+abs(real(fockMO(j,i),realk))
			error(k)=error(k)*0.5
		ENDDO
	ENDDO

END SUBROUTINE  pbc_geterrorvec

! todo not in use (delete?) 
!> \author 
!> \date
!> \brief
!> \param
!> \param
SUBROUTINE pbc_kcomputeenergy()
  IMPLICIT NONE
END SUBROUTINE pbc_kcomputeenergy

! todo not in use (delete?) 
!> \author 
!> \date
!> \brief
!> \param
!> \param
SUBROUTINE pbc_dcomputeenergy()
  IMPLICIT NONE
END SUBROUTINE pbc_dcomputeenergy

!> \author JR 
!> \date 2013
!> \brief Solve the generalized eigenvalueproblem for FC=SCF
!> S is not input but Uk where S^-1=Uk Uk^H. The generalised e.v. problem is 
!> solved by calculating f = Uk^H F Uk and then solving fc = ce. C is now found by
!> calculating Uk c. See fex. Golub and Loan matrix comput on the symmetric
!> generalized eigenvalueproblem for a better explanation.
!> \param ndim 		matrix dim.
!> \param fock_old 	fock matrix (F).
!> \param c_tmp 		Output C.
!> \param Uk 			The unitary transformation.
!> \param eigv 		Output eigenvalues.
!> \param nsingular 	dim - rank of Uk.
!> \param lupri 		Logical print unit.
SUBROUTINE solve_kfcsc_mat(ndim,fock_old,C_tmp,Uk,eigv,nsingular,lupri)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: ndim,lupri,nsingular
	COMPLEX(complexk),INTENT(IN) :: fock_old(ndim,ndim)
	COMPLEX(complexk), INTENT(IN) :: Uk(ndim-nsingular,ndim-nsingular)
	COMPLEX(complexk),intent(INOUT) :: C_tmp(ndim,ndim)
	REAL(realk),intent(INOUT) :: eigv(ndim-nsingular)
	! local
	INTEGER :: lwork,info,lrwork,nonsingdim
	COMPLEX(complexk),POINTER :: tfock(:,:),work(:)
	REAL(realk),POINTER :: rwork(:)
	COMPLEX(COMPLEXK) :: alpha,beta

	nonsingdim=ndim-nsingular

	alpha=CMPLX(1._realk,0._realk,complexk)
	beta=CMPLX(0._realk,0._realk,complexk)

	lwork=2*Ndim-1
	lrwork=3*Ndim-2

	call mem_alloc(rwork,lrwork)
	call mem_alloc(work,max(1,lwork))
	call mem_alloc(tfock,nonsingdim,nonsingdim)

	! Transform f = U^T F U
	call pbc_unitary_transform(Ndim,nonsingdim,fock_old,Uk,tfock)
	! solve fc = ce 
	call zheev('V','U',nonsingdim,tfock,nonsingdim,eigv,work,lwork,rwork,info)
	if(info .ne. 0) then
		write(lupri,*) 'ERROR: zheev problems, info=', info
		call write_zmatrix(tfock,nonsingdim,nonsingdim)
		write(*,*) 'ERROR: zheev problems, info=', info
		call LSQUIT('pbc_solve_kfcsc_mat: INFO not zero,', &
			& ' while solving eigenvalue',lupri)
	endif
	!Transform back to C=Uc
	call zgemm('N','N',ndim,nonsingdim,nonsingdim,alpha,Uk,ndim,&
		tfock,nonsingdim,beta,C_tmp(:,1:nonsingdim),ndim)
	c_tmp(:,nonsingdim+1:ndim)=cmplx(0._realk,0._realk,complexk)

	call mem_dealloc(work)
	call mem_dealloc(rwork)
	call mem_dealloc(tfock)

END SUBROUTINE solve_kfcsc_mat

!Todo not in use. Only called from pbc_k_energy which is not in use. delete?
!> \author JR 
!> \date 2013
!> \brief Compute the energy for a k vector
!> \param
!> \param
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

!Todo not in use. (delete?)
!> \author JR 
!> \date 2013
!> \brief Compute the energy for one cell
!> \param
!> \param
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
 onepenergy=onepenergy+2.*bz%kpnt(kpt)%eigv(i)
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

!> \author JR 
!> \date 2013
!> \brief Does the SCF iterations and solces the F(k)C(k)=eps(k)C(k) for each k
!> \param
!> \param
SUBROUTINE pbc_startzdiis(molecule,setting,ndim,lattice,numrealvec,&
           maxmultmom,bz,dmat0,lupri,luerr)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,lupri,luerr,numrealvec,maxmultmom
  TYPE(lvec_list_t),INTENT(INOUT) :: lattice
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule
  TYPE(LSSETTING) :: setting
  TYPE(BZgrid_t),intent(inout) :: bz
  TYPE(matrix),INTENT(IN) :: dmat0
!  REAL(realk),INTENT(IN) :: E_nuc
!  TYPE(matrix),target,intent(inout) :: f_1(numrealvec),ovl(numrealvec)
!  TYPE(matrix),target,intent(inout) :: g_2(numrealvec)
  ! local
  COMPLEX(COMPLEXK),pointer :: fockMO(:,:),fock(:,:),smatk(:,:)
  COMPLEX(COMPLEXK),pointer :: C_k(:,:),D_k(:,:),C_0(:,:)
  REAL(realk),pointer :: cellenergies(:)
  REAL(realk),pointer :: error(:,:),nucmom(:)
  INTEGER :: i,errlm,tol!tol should be a real but now I just have it for test
  INTEGER :: k,kpt,n1,fdim(3),layer
  INTEGER :: realcut(3)
  real(realk) :: kvec(3),Ecell,E_nn,E_1,errortest
  REAL(realk),pointer :: tlat(:,:),weight(:)
  TYPE(matrix), pointer :: nfdensity(:)
  TYPE(moleculeinfo),pointer :: latt_cell(:)
  TYPE(moleculeinfo) :: refcell
  CHARACTER(LEN=10) :: stiter
  CHARACTER(LEN=12) :: diis,diismats
  LOGICAL :: diis_exit
  REAL(realk) :: E_J,E_K,E_XC,E_ff,E_cell
  REAL(realk) :: E_en,E_nuc,E_kin,E_nnff
  real(realk)         :: TS,TE,TST,TET,TOT,TWT !For finding time usage
  TYPE(matrix),pointer :: f_1(:),ovl(:)
  TYPE(matrix),pointer :: g_2(:)


  ! threshhold for removing singularities in overlap matrix s. If any of the
  ! eigenvalues of the overlapmtrx S are smaller than this number the
  ! corresponding eigenvectors are removed.
  REAL(realk) :: singular_threshh = 1e-6_realk 

  write(lupri,*) 'Entering routine startzdiis'

    write(stiter,'(I5)') 1
    stiter=adjustl(stiter)
    diis='diis_'//trim(stiter)//'_'

  call mem_alloc(tlat,(maxmultmom+1)**2,(maxmultmom+1)**2)

  !read the tlatticetensor to build up the multipole moments
!  write(lupri,*) 'Debug again'
!  call read_pbc_tlattice(tlat,maxmultmom,'Tlatticetensor.dat',lupri)
  !config%molecule%nelectrons

!  call mem_alloc(nfdensity,numrealvec)
  !allocate(nfdensity(numrealvec))

  call mem_alloc(nfdensity,numrealvec)
  call find_latt_index(n1,0,0,0,fdim,lattice,lattice%max_layer)
  call mat_init(nfdensity(n1),ndim,ndim)
  call mat_copy(1.0_realk,Dmat0,nfdensity(n1)) 
  call mem_alloc(f_1,numrealvec)
  call mem_alloc(Ovl,numrealvec)
  call mem_alloc(g_2,numrealvec)
  call mem_alloc(fockMO,ndim,ndim)
  call mem_alloc(fock,ndim,ndim)
  call mem_alloc(smatk,ndim,ndim)
  call mem_alloc(cellenergies,ndim)
  call mem_alloc(C_k,ndim,ndim)
  call mem_alloc(C_0,ndim,ndim)
  call mem_alloc(D_k,ndim,ndim)

  allocate(latt_cell(numrealvec)) !fixme alloc ?? use mem alloc el.

  call set_refcell(refcell,molecule)
  call set_lattice_cells(latt_cell,numrealvec,molecule,lattice,lupri)
  errlm=molecule%nelectrons/2*(ndim-molecule%nelectrons/2.)
  call mem_alloc(error,lattice%num_store,errlm)


  call LSTIMER('START ',TS,TE,LUPRI)
  call pbc_overlap_k(lupri,luerr,setting,molecule,ndim,&
      lattice,latt_cell,refcell,numrealvec,ovl)
  call LSTIMER('pbc_overlap_k',TS,TE,LUPRI)
 
  !CALCULATES kinetic energy of electrons
  call LSTIMER('START ',TS,TE,LUPRI)
  call pbc_kinetic_k(lupri,luerr,setting,molecule,ndim,&
   lattice,latt_cell,refcell,numrealvec,nfdensity,f_1,E_kin)
  call LSTIMER('pbc_kinetic_k',TS,TE,LUPRI)

  !CALCULATES electron nuclei attraction
  call LSTIMER('START ',TS,TE,LUPRI)
  call pbc_nucattrc_k(lupri,luerr,setting,molecule,ndim,&
     lattice,latt_cell,refcell,numrealvec,nfdensity,f_1,E_en)
  call LSTIMER('pbc_nucattrc_k',TS,TE,LUPRI)

  !CALCULATES nuclear repulsion
  call LSTIMER('START ',TS,TE,LUPRI)
  call pbc_nucpot(lupri,luerr,setting,molecule,lattice,&
                  latt_cell,refcell,numrealvec,E_nuc)
  call LSTIMER('pbc_nucpot',TS,TE,LUPRI)

 !computing the nuclear moments
  call mem_alloc(nucmom,(1+maxmultmom)**2)
  call pbc_comp_nucmom(refcell,nucmom,maxmultmom,lupri)

  CALL LSTIMER('START ',TS,TE,LUPRI)
  call pbc_multipole_expan_k(lupri,luerr,setting,ndim,lattice,&
    &latt_cell,refcell,numrealvec,maxmultmom)
  CALL LSTIMER('pbc_multipole',TS,TE,LUPRI)
  
  call pbc_controlmm(20,Tlat,lattice%Tlmax,maxmultmom,.false.,lattice%ldef%avec,&
     ndim,lupri,nfdensity,numrealvec,lattice,E_ff,E_nnff,refcell)

!  call pbc_controlmm(20,Tlat,Tlmax,maxmultmom,.false.,lattice%ldef%avec,&
!     nbast,lupri,nfdensity,num_latvectors,nfsze,lattice,g_2,E_ff,E_nnff,refcell)

   DO k=1,BZ%nk
!
      !call mem_alloc(bz%kpnt(k)%Uk,ndim,ndim)
      !call mem_alloc(bz%kpnt(k)%Uinv,ndim,ndim)
      call pbc_get_kpoint(k,kvec)

      if(lattice%store_mats)then
        !get the overlap matrices S^0l
        call pbc_read_matrix(lattice,ndim,ndim,1,1,'            ')
        !Transform overlap matrix 0l to k
        call pbc_rspc_to_kspc_mat(lattice,Bz,ndim,kvec,1)
      else
        ! transforms overlap to kspace
        call pbc_trans_mat_to_kspc(ovl,numrealvec,lattice,Bz,ndim,kvec,realcut)
      endif
!
!      call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)
!
      call pbc_zdevectorize_mat(smatk,ndim,ndim,bz%smat%zelms)
!
      !diagonalizes Sk, Uk transform operators
      call pbc_spectral_decomp_ovl(smatk,bz%kpnt(k)%Uk, & 
           & bz%kpnt(k)%is_singular,Ndim,bz%kpnt(k)%nsingular,singular_threshh,lupri)

      call mem_alloc(bz%kpnt(k)%eigv,ndim-bz%kpnt(k)%nsingular)

   ENDDO
!

  k=0
  i=0 !either 0 or 1, check it, for zero no errors
  tol=0

  Ecell=0.0_realk
  diis_exit = .false. !when diis_exit the iterations are finished
  ! self consistent iterations
  DO WHILE(tol .le. lattice%num_its)! 20)!should have an input parameter
    call LSTIMER('START ',TOT,TWT,LUPRI)
    k=k+1
    i=i+1

    !We keep only data of lattice%num_store past iterations
    if(i .ge. lattice%num_store) Then !should have an input parameter to decide how man we store
      k=lattice%num_store-1
      i=lattice%num_store
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

    call LSTIMER('START ',TST,TET,LUPRI)
    call LSTIMER('START ',TS,TE,LUPRI)
    call pbc_electron_rep_k(lupri,luerr,setting,molecule,ndim,&
      lattice,latt_cell,refcell,numrealvec,nfdensity,g_2,E_J)
    call LSTIMER('pbc Coul',TS,TE,LUPRI)

    !if hybrid or HF, include parameter if hybrid
    call LSTIMER('START ',TS,TE,LUPRI)
    call pbc_exact_xc_k(lupri,luerr,setting,molecule,ndim,&
     lattice,latt_cell,refcell,numrealvec,nfdensity,g_2,E_K)
    call LSTIMER('pbc xchange',TS,TE,LUPRI)
    call LSTIMER('rep xchange',TST,TET,LUPRI)


    !KOHN sham

    lattice%fc1=max(lattice%oneop1,lattice%col1)
    lattice%fc1=max(lattice%fc1,lattice%Kx1)
    lattice%fc2=max(lattice%oneop2,lattice%col2)
    lattice%fc2=max(lattice%fc2,lattice%Kx2)
    lattice%fc3=max(lattice%oneop3,lattice%col3)
    lattice%fc3=max(lattice%fc3,lattice%Kx3)
    !lattice%compare_elmnts=.false.

    !Far-field contribution to fock
    call LSTIMER('START ',TS,TE,LUPRI)
    call pbc_ff_fck(maxmultmom,tlat,lattice%tlmax,ndim,lattice,nfdensity,nucmom,&
                   g_2,E_ff,E_nn,lupri)
    call LSTIMER('pbc farfield',TS,TE,LUPRI)

  !sums the parts for the fock matrices from  f^0l
  call pbc_get_fock_mat(lattice,g_2,f_1,ndim,realcut,numrealvec,diismats,lupri)
  if(.not.lattice%store_mats) call pbc_fockmat_write(lattice,ndim,ndim,7,2,diismats,lupri)
  !Obsolete 
  realcut(1)=lattice%oneop1
  realcut(2)=lattice%oneop2
  realcut(3)=lattice%oneop3

  DO layer=1,numrealvec  
     if(nfdensity(layer)%init_magic_tag.EQ.mat_init_magic_value) THEN
       call mat_zero(nfdensity(layer))
     endif
  ENDDO


  !No we have the fock matrices in real space
  !We transform to kspace and solve
  call LSTIMER('START',TST,TET,LUPRI)
  do kpt=1,bz%nk

	  call pbc_get_kpoint(kpt,kvec)
	  call zero_pbc_elstr(Bz%fck)
	  call zero_pbc_elstr(Bz%Smat)

	  ! We compute only diis weights for the gamma point
	  if(bz%kpnt(kpt)%is_gamma)then
		  call LSTIMER('START',TS,TE,LUPRI)
		  call mem_alloc(weight,i)
		  weight(:)=0._realk
		  call pbc_get_diisweights(lattice,Bz,weight,i,tol,kvec,ndim,C_0, &
			  & fockMO,fock,numrealvec,errortest,error,diis_exit,errlm, &
			  & molecule%nelectrons,lupri)
		  call LSTIMER('diis weights',TS,TE,LUPRI)
	  endif ! is_gamma

	  if(tol.gt.0)then
		  if(.not.diis_exit)then
			  !call LSTIMER('START',TS,TE,LUPRI)
			  call pbc_get_weighted_fock(i,tol,lattice%num_store,ndim,weight, &
				  & lattice,lupri)
			  !call LSTIMER('Weighted fock',TS,TE,LUPRI)
		  endif
	  endif

	  if(lattice%store_mats)then
		  !get the overlap matrices S^0l
		  call pbc_read_matrix(lattice,ndim,ndim,1,1,'            ')
		  !We need k-space overlap
		  call pbc_rspc_to_kspc_mat(lattice,Bz,ndim,kvec,1)
	  else
		  !We need k-space overlap
		  call pbc_trans_mat_to_kspc(ovl,numrealvec,lattice,Bz,ndim,kvec,realcut)
	  endif

	  !Put overlap in a matrix form
	  !write(lupri,*) 'Testing for C coefficients in gamma point'
	  call pbc_zdevectorize_mat(smatk,ndim,ndim,bz%smat%zelms)
	  !We need to get k point fock matrix
	  call pbc_rspc_to_kspc_mat(lattice,Bz,ndim,kvec,2)
	  !Put it in a matrix form
	  call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)
	  !solves F(k)C(k)=S(k)C(k)e(k)
	  call solve_kfcsc_mat(ndim,fock,C_k,bz%kpnt(kpt)%Uk,bz%kpnt(kpt)%eigv, &
		  & bz%kpnt(kpt)%nsingular,lupri)
	  !C_0 is used for finding the weights
	  if(bz%kpnt(kpt)%is_gamma ) C_0(:,:)=C_k(:,:)
	  !call pbc_zdevectorize_mat(smatk,ndim,ndim,bz%smat%zelms)

	  !gets D(k) from C(k) (only working for unrestricted)
	  call pbc_get_kdensity(D_k,C_k,ndim,molecule%nelectrons/2 &
		  ,bz%kpnt(kpt)%nsingular,smatk,lupri)
	  !Converts D(k) to D^0l
	  call kspc_2_rspc_loop_k(nfdensity,Bz%Nk,D_k,lattice,kvec,bz%kpnt(kpt)%weight,BZ%NK_nosym,ndim,kpt)
  enddo !kpt

  call LSTIMER('k point energy',TST,TET,LUPRI)

  if(associated(weight)) call mem_dealloc(weight)

  write(*,*)
  call pbc_densitymat_write(nfdensity,lattice,ndim,ndim,8,'            ')
  call pbc_free_read_matrices(lattice)
  call print_bands(bz,ndim,'band-energy') !prints band energy to file band-energy
  ! Get HOMO LUMO energy, change name 
  call pbc_trans_k_energy(lattice,cellenergies,ndim,molecule%nelectrons,bz)

  E_cell=E_1+E_j+E_K+E_ff+E_nuc
  write(lupri,*) 'E(HOMO) =', cellenergies(1), tol
  write(lupri,*) 'E(LUMO) =', cellenergies(2), tol
  write(lupri,*) 'Cell Energy electrons=', E_cell
  write(lupri,*) 'K energy', E_k
  write(lupri,*) 'J energy', E_J
  write(lupri,*) 'h_1=',E_1
  write(lupri,*) 'Nuclear=',E_nuc
  write(lupri,*) 'Far field=', E_ff, E_nn
  write(*,*) 'H_1=',E_1
  write(*,*) 'E(HOMO) =', cellenergies(1)
  write(*,*) 'E(LUMO) =', cellenergies(2)
  write(*,*) 'Cell Energy =', E_cell
  write(*,*) 'K energy', E_k
  write(*,*) 'J energy', E_J
  write(*,*) 'h_1=',E_1
  write(*,*) 'Nuclear=',E_nuc
  write(*,*) 'Far field=', E_ff

  if(diis_exit) exit

  tol=tol+1
  call LSTIMER('Diis Iteration',TOT,TWT,LUPRI)

  ENDDO

  do k=1,bz%Nk
     call mem_dealloc(bz%kpnt(k)%Uk)
#ifdef DEBUGPBC
     call mem_dealloc(bz%kpnt(k)%Uinv)
#endif
  enddo
    

  do i=1,numrealvec
	call free_Moleculeinfo(latt_cell(i))
      if(nfdensity(i)%init_magic_tag.NE.mat_init_magic_value) CYCLE
	call mat_free(nfdensity(i))
  enddo
  call free_Moleculeinfo(refcell)
  call mem_dealloc(tlat)
  call mem_dealloc(nfdensity)
  call mem_dealloc(nucmom)
  call mem_dealloc(error)
  call mem_dealloc(fockMO)
  call mem_dealloc(fock)
  call mem_dealloc(smatk)
  call mem_dealloc(C_k)
  call mem_dealloc(C_0)
  call mem_dealloc(D_k)
  deallocate(latt_cell)


	do i=1,numrealvec
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
  write(lupri,*) 'h_1=',E_1
  write(lupri,*) 'Nuclear=',E_nuc
  write(lupri,*) 'Far field=', E_ff,E_nn
  write(lupri,*) 'K energy', E_k
  write(lupri,*) 'J energy', E_J
  write(*,*) 'E(HOMO) =', cellenergies(1)
  write(*,*) 'E(LUMO) =', cellenergies(2)
  write(*,*) 'Cell Energy =', E_cell
  write(*,*) 'K energy', E_k
  write(*,*) 'J energy', E_J
  write(*,*) 'h_1=',E_1
  write(*,*) 'Nuclear=',E_nuc
  write(*,*) 'Far field=', E_ff,E_nn
  call mem_dealloc(cellenergies)

END SUBROUTINE pbc_startzdiis



!> \author JR 
!> \date 2013
!> \brief 
!> \param
!> \param
SUBROUTINE pbc_get_fock_mat(lattice,g_2,f_1,ndim,realcut,numrealvec,diismats,lupri)
  IMPLICIT NONE
  TYPE(lvec_list_t),INTENT(INOUT) :: lattice
  INTEGER,INTENT(IN) :: ndim,lupri,numrealvec
  TYPE(matrix),target,intent(inout) :: f_1(numrealvec)
  TYPE(matrix),target,intent(inout) :: g_2(numrealvec)
  CHARACTER(LEN=12) :: diismats
  INTEGER,intent(OUT) :: realcut(3)
  !LOCAL VARIABLES
  INTEGER             :: i,j
  real(realk)         :: focknorm !For finding time usage
  real(realk)         :: TS,TE !For finding time usage


  !get the fock matrices f^0l
  if(lattice%store_mats) then
    call LSTIMER('START',TS,TE,LUPRI)
    !call pbc_read_fock_matrix(lattice,ndim,ndim,diismats)
    call pbc_read_fock_matrix(lattice,ndim,ndim,'            ')
    call LSTIMER('Reading fock',TS,TE,LUPRI)
    call pbc_fockmat_write(lattice,ndim,ndim,7,2,diismats,lupri)
    !get the overlap matrices S^0l
    call pbc_read_matrix(lattice,ndim,ndim,1,1,'            ')
  else
    call pbc_add_fock_matrix(f_1,g_2,lattice,ndim,ndim,numrealvec) 
    !call pbc_fockmat_write(lattice,ndim,ndim,7,2,diismats,lupri) !fixme
  endif

END SUBROUTINE pbc_get_fock_mat

!Todo not in use. (delete?)
!> \author JR 
!> \date 2013
!> \brief 
!> \param
!> \param
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
w=0._realk
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

!winv(:,:)=0._realk
do j=1,nrow
winv(j,j,i)=CMPLX(1./sqrt(w(j,i)),0._realk,complexk)
enddo
!write(*,*) w(:,i)
alpha=CMPLX(1._realk,0._realk,complexk)
beta=CMPLX(0._realk,0._realk,complexk)
call zgemm('N','N',nrow,ncol,ncol,alpha,winv,nrow,kdep_tmp(i)%koverlapmat,&
                                                     nrow,beta,sk_tmp,nrow)

call zgemm('C','N',nrow,ncol,ncol,alpha,kdep_tmp(i)%koverlapmat,&
                        nrow,sk_tmp,nrow,beta,kdep_tmp(i)%kddensitymat,nrow)

!call zgemm('T','N',nrow,ncol,ncol,alpha,V_tmp,nrow,V_tmp,nrow,beta,DMAT,nrow)
!write(*,*) kdep_tmp(i)%kddensitymat
ENDDO

!write(*,*) 'nfsize in initial',nfsze
!call pbc_trans_to_realspc(lattice,nfdensity,nfsze,ndim,bz,kdep_tmp)
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
!alpha=1._realk
!beta=0._realk
!call dgemm('N','N',nrow,ncol,ncol,alpha,winv,nrow,V,nrow,beta,V_tmp,nrow)
!call dgemm('T','N',nrow,ncol,ncol,alpha,V,nrow,V_tmp,nrow,beta,V_tmp,nrow)
!call dgemm('T','N',nrow,ncol,ncol,alpha,V_tmp,nrow,V_tmp,nrow,beta,DMAT,nrow)
!call write_matrix(DMAT,nrow,ncol)

END SUBROUTINE pbc_get_initial_density

!> \author JR 
!> \date 2013
!> \brief Calc the electronic one body energy 
!> \param
!> \param
SUBROUTINE  pbc_get_onehamenergy(numvecs,f_1,nfdensity,E_1)
INTEGER,INTENT(IN) :: numvecs
TYPE(matrix),target,intent(IN) :: f_1(numvecs),nfdensity(numvecs)
REAL(REALK),INTENT(INOUT) ::  E_1
!LOCAL
INTEGER :: celli

 E_1=0
 do celli=1,numvecs
    if(f_1(celli)%init_magic_tag.EQ.mat_init_magic_value .and. &
     & nfdensity(celli)%init_magic_tag.EQ.mat_init_magic_value ) THEN
      E_1=E_1+mat_dotproduct(f_1(celli),nfdensity(celli))
    endif
 enddo

END SUBROUTINE  pbc_get_onehamenergy

!> \author JR 
!> \date 2013
!> \brief 
!> \param
!> \param
subroutine pbc_get_kdensity(ddensity,c_tmp,nbast,nkmobas,nsingular,smatk,lupri)
  implicit none
  integer,intent(in) :: nbast,nkmobas,lupri,nsingular
  complex(complexk),intent(inout) :: ddensity(nbast,nbast)
  complex(complexk),intent(in) :: c_tmp(nbast,nbast),smatk(nbast,nbast)
  real(realk) :: nelectrons
  ! local
  complex(complexk), pointer :: tmp(:,:)
  integer :: i,nosingdim
  complex(complexk) :: alpha,beta

  alpha=cmplx(2._realk,0._realk,complexk)
  beta =cmplx(0._realk,0._realk,complexk)

  nosingdim=nbast-nsingular
  ddensity(:,:)=cmplx(0._realk,0._realk,complexk)

  ! D = C C^H
  call zgemm('n','c',nbast,nbast,nkmobas,alpha,c_tmp,nbast,&
	  c_tmp,nbast,beta,ddensity,nbast)

  !==================================
  !==================================
  ! DEBUG (smathk,tmp, only used here)
  ! test number of electrons
  ! DELETE THE REST
  !==================================
  !==================================

  call mem_alloc(tmp,nbast,nbast)
  alpha=cmplx(1._realk,0._realk,complexk)

  call zgemm('n','n',nbast,nbast,nbast,alpha,ddensity,nbast,&
	  smatk,nbast,beta,tmp,nbast)
  nelectrons=0._realk
  do i=1,nbast
	  nelectrons=nelectrons+real(tmp(i,i))
  end do
  write(*,*) 'nelectrons =', nelectrons,nkmobas,nsingular
  write(lupri,*) 'nelectrons =', nelectrons,nkmobas,nsingular

  call mem_dealloc(tmp)

end subroutine pbc_get_kdensity

!> \author JR 
!> \date 2013
!> \brief 
!> \param
!> \param
SUBROUTINE pbc_get_diisweights(lattice,Bz,weight,its,tol,kvec,ndim,C_0,fockMO,fock,&
                numrealvec,errortest,error,diis_exit,errdim,nelectrons,lupri)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,lupri,numrealvec,its,tol,errdim
  INTEGER,INTENT(IN) :: nelectrons
  TYPE(lvec_list_t),INTENT(INOUT) :: lattice
  Real(realk),intent(OUT) :: weight(its)
  REAL(realk),intent(inout) :: error(lattice%num_store,errdim)
  LOGICAL,intent(out) :: diis_exit
  REAL(realk),intent(out) :: errortest
  real(realk),intent(in) :: kvec(3)
  TYPE(BZgrid_t),intent(inout) :: bz
  COMPLEX(COMPLEXK) :: C_0(ndim,ndim)
  COMPLEX(COMPLEXK),intent(out) :: fockMO(ndim,ndim)
  COMPLEX(COMPLEXK),intent(inout) :: fock(ndim,ndim)
  !LOCAL
  !Real(realk) :: tol
  INTEGER :: j

      !We need k-space fock matrix in gamma point now
      call pbc_rspc_to_kspc_mat(lattice,Bz,ndim,kvec,2)

      call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)
        !write(*,*) pbc_it(i)%kdep_it(1)%kfockvec
        !write(*,*) 'fock(0)'
        !call write_zmatrix(fock,ndim,ndim)
        !write(lupri,*) 'fock(0)'
        !call write_zmatrix(fock,ndim,ndim,lupri)
!
        if(tol .ge. 1) then
          call transform_toMOfock(C_0,fock,fockMO(:,:),ndim,lupri)

        !call mem_alloc(weight,i)
        weight=0._realk
        !! k is still the gamma point
        if(tol .gt. lattice%num_store) then
          do j=1,lattice%num_store-1
           error(j,:)=error(j+1,:)
          enddo
          error(lattice%num_store,:)=0._realk
        endif

        !get the error vectors
        if(its .gt. tol) then ! this since I do not know the C0 matrix for it 0
          call pbc_geterrorvec(error(its-1,:),fockMO(:,:),ndim,errdim,nelectrons)
          errortest=dot_product(error(its-1,:),error(its-1,:))
        else
          call pbc_geterrorvec(error(its,:),fockMO(:,:),ndim,errdim,nelectrons)
          errortest=dot_product(error(its,:),error(its,:))
        endif
        errortest=sqrt(errortest)
        if(errortest .le. lattice%error) diis_exit=.true.

        !Get diis weights  !!!THIS HAS TO BE FIXED, its-1 only when its .gt.
        !tol
        if(its .gt. tol) then ! this since I do not know the C0 matrix for it 0
          call pbc_diisweights(errdim,error,weight,its-1,lattice%num_store,lupri)
        else
          call pbc_diisweights(errdim,error,weight,its,lattice%num_store,lupri)
        endif
        endif

        if(tol .eq. 0 .or. tol .eq. 1) weight(1)=1.0_realk

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
        !call pbc_get_weighted_fock(i,tol,7,ndim,weight,lattice)
        !write(*,*) 'Debug after get_weights'

        !We need to get gamma point fock matrix again
        !call pbc_rspc_to_kspc_mat(lattice,Bz,ndim,kvec,2)

        !call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)

        ! We have to reset the gamma point again since we use it again.
        call zero_pbc_elstr(Bz%fck)
        call zero_pbc_elstr(Bz%Smat)

END SUBROUTINE pbc_get_diisweights

!> \author JR 
!> \date 2013
!> \brief 
!> \param
!> \param
SUBROUTINE pbc_diisweights(errdim,error,weight,it,num_store,lupri)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: errdim,lupri,it,num_store
  REAL(realk),INTENT(INOUT) :: error(num_store,errdim)
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

  B_mat(:,:)=0._realk
  DO i=1,it
    write(lupri,*) 'error(',i,')'
    write(lupri,*) error(i,:)
    !B_mat(i,i)= 0.05_realk
   DO j=1,it
    B_mat(i,j)=B_mat(i,j)+dot_product(error(i,:),error(j,:))! error(i)*error(j)
   ENDDO
  ENDDO
!  call mat_init(Bmat_t,it,it)
!  call mat_zero(Bmat_t)
!  call mat_set_from_full(B_mat,1._realk,Bmat_t)

!  B_mat(it+1,it+1)=0._realk
  solve=0
  !weight_tmp(it+1,1)=-1._realk
  weight_tmp(:,1)=1._realk

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
      ! B_mat(it+1,j)=-1._realk
      ! B_mat(j,it+1)=-1._realk
       !weight_tmp(j,1)=0._realk
       weight_tmp(j,1)=1._realk
      ENDDO
      !B_mat(it+1,it+1)=0._realk
      !weight_tmp(it+1,1)=-1._realk
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

!> \author JR 
!> \date 2013
!> \brief Sums over former Fock matrices with the cerresponding weights. 
!> \param
!> \param
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
     if(Aop%lvec(layer)%f1_computed .or. Aop%lvec(layer)%g2_computed)then
       call mat_zero(Aop%lvec(layer)%oper(2))
     endif
  ENDDO

  if(ndiis .le. stdiis) then
    do j=1,ndiis
    write(stiter,'(I5)') j
    stiter=adjustl(stiter)
    tmpdiis='diis_'//trim(stiter)//'_'
!    write(*,*) tmpdiis

    call pbc_read_matrix(tmp_mat,nbast,nbast,7,2,tmpdiis)

    Do layer=1,size(Aop%lvec)
     l1=int(Aop%lvec(layer)%lat_coord(1))
     l2=int(Aop%lvec(layer)%lat_coord(2))
     l3=int(Aop%lvec(layer)%lat_coord(3))
     if(Aop%lvec(layer)%f1_computed .or. Aop%lvec(layer)%g2_computed)then
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
!    write(lupri,*) 'Filename for fock matrix, ', tmpdiis

    call pbc_read_matrix(tmp_mat,nbast,nbast,7,2,tmpdiis)
    !write(*,*) 'Debug 4 inside get_weights'

    Do layer=1,size(Aop%lvec)
     l1=int(Aop%lvec(layer)%lat_coord(1))
     l2=int(Aop%lvec(layer)%lat_coord(2))
     l3=int(Aop%lvec(layer)%lat_coord(3))
     if(Aop%lvec(layer)%f1_computed .or. Aop%lvec(layer)%g2_computed)then
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
   if(Aop%lvec(layer)%f1_computed .or. Aop%lvec(layer)%g2_computed)then
!!   ! write(*,*) 'Debug 5 inside get_weights',l1
     call pbc_get_file_and_write(Aop,nbast,nbast,layer,7,2,'            ')
!     call pbc_get_file_and_write(Aop,nbast,nbast,layer,7,2,diis)
   endif
   call mat_free(tmp_mat%lvec(layer)%oper(2))
  Enddo
   !  call pbc_fockmat_write(Aop,nbast,nbast,7,2)
   ! write(*,*) 'Debug 6 inside get_weights'

  !deallocate(tmp_mat%lvec)
  call mem_dealloc(tmp_mat%lvec)
END SUBROUTINE pbc_get_weighted_fock


!> \author JR 
!> \date 2013
!> \brief 
!> \param
!> \param
SUBROUTINE pbc_trans_k_energy(lattice,cenergies,nbast,nelectrons,bz)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nbast,nelectrons
  TYPE(lvec_list_t),intent(IN) :: lattice
  TYPE(BZgrid_t),intent(in) :: bz
!  TYPE(pbc_elstr_t),INTENT(IN) :: kdep_tmp(bz%nk)
  REAL(realk),intent(INOUT) :: cenergies(nbast)
  !LOCAL VARIABLElS
  INTEGER :: lattindex(3)!,irealspc
  INTEGER :: i,kpt
  REAL(Realk) :: ehomo,elumo,etmph1,etmph2,etmpl1,etmpl2

  lattindex(1)=0
  cenergies =0.0_realk
  ehomo=-huge(ehomo)
  elumo= huge(elumo)
  !write(*,*) lattice%ldef%is_active(1)
  !write(*,*) lattice%ldef%is_active(2)
  !write(*,*) lattice%ldef%is_active(3)

     lattindex(1)=0
     lattindex(2)=0
     lattindex(3)=0
     
     cenergies(1)=bz%kpnt(1)%eigv(nelectrons/2)
     cenergies(2)=bz%kpnt(1)%eigv(nelectrons/2+1)
   
    DO kpt=2,bz%nk
       etmph1=bz%kpnt(kpt)%eigv(nelectrons/2)
       etmph2=bz%kpnt(kpt-1)%eigv(nelectrons/2)
       etmpl1=bz%kpnt(kpt)%eigv(nelectrons/2+1)
       etmpl2=bz%kpnt(kpt-1)%eigv(nelectrons/2+1)
       !etmph1=kdep_tmp(kpt-1)%keigv(nelectrons/2)
       !etmph2=kdep_tmp(kpt)%keigv(nelectrons/2)
       !etmpl1=kdep_tmp(kpt-1)%keigv(nelectrons/2+1)
       !etmpl2=kdep_tmp(kpt)%keigv(nelectrons/2+1)
       ehomo=max(etmph1,etmph2)
       elumo=min(etmpl1,etmpl2)
    ENDDO
       cenergies(1)=max(ehomo,cenergies(1))
       cenergies(2)=min(elumo,cenergies(2))

!     DO kpt=1,bz%nk
!      do i=1,nelectrons/2
!       cenergies(3)=cenergies(3)+2._realk*kdep_tmp(kpt)%keigv(i)*&
!       &bz%kpnt(kpt)%weight/BZ%NK_nosym
!    enddo
!   enddo
!!       cenergies(3)=behomo
!       cenergies(4)=belumo


END SUBROUTINE  pbc_trans_k_energy


!Todo not in use (delete?)
!> \author JR 
!> \date 2013
!> \brief 
!> \param
!> \param
SUBROUTINE pbc_trans_to_realspc(lattice,nfdensity,nvecsrs,nbast,bz,kdep_tmp,lupri)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvecsrs,nbast,lupri
  TYPE(lvec_list_t),intent(IN) :: lattice
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

  !write(*,*) lattice%ldef%is_active(1)
  !write(*,*) lattice%ldef%is_active(2)
  !write(*,*) lattice%ldef%is_active(3)
  DO irealspc=1,nvecsrs

     call find_latt_vectors(irealspc,l1,l2,l3,fdim,lattice)
     call mat_zero(nfdensity(irealspc))
     if(abs(l1) .gt. lattice%ndmat) CYCLE
     if(abs(l2) .gt. lattice%ndmat) CYCLE
     if(abs(l3) .gt. lattice%ndmat) CYCLE
     call transformk_2_realmat(kdep_tmp,bz,&
          dkmat,nbast,lattice%lvec(irealspc)%std_coord,lupri)
     
     call mat_set_from_full(dkmat,1.0_realk,nfdensity(irealspc))
     
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


!> \author JR 
!> \date 2013
!> \brief The LAPACK routine ZGGEV gives back gen. eig. vectors normalized in an
!> \brief odd way. This subroutine simply fixes the normalization.
!> \param
!> \param
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
     normfac = 0.0_realk
     coeff_sum = 0.0_realk
     coeff_max = 0.0_realk
     do i = 1,siz
!        call lsquit('FIXME: coeff_max real while cocoeff(i,ivec) is complex',-1)
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
        !betavec(ivec) = 1.0_realk
        normfac = 1.0D50
        !call lsquit('Zero norm vector - eeeeeh',-1)
     end if
#endif

     cocoeff(:,ivec) = cocoeff(:,ivec) / sqrt(normfac)
  end do

!Careful
  n_occ = ( 0.0_realk, 0.0_realk )
  loop_ao1: do mu = 1,siz
     loop_ao2: do nu = 1,siz
        dmat_elem = ( 0.0_realk, 0.0_realk )
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

