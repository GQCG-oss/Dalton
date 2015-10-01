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
  USE pbc_timings
	
  PRIVATE  					
  PUBLIC :: pbc_startzdiis	

  CONTAINS

!> \author JR 
!> \date 2013
!> \brief For a matrix pos. semdefinite matrix S : calculates U = V\sigma^{-1/2} 
!> \brief where V\sigmaV^H = S. Removes singularities that are smaller than 
!> \brief singular_threshh.
!> \param Sabk 					ndim x ndim complex matrix
!> \param U 						Changed to V\sigma
!> \param is_singular 			Set to true if singularities found
!> \param ndim 					Matrix dim
!> \param nsingular 				ndim - rank(Sabk)
!> \param singular_threshh 	Singularity threshhold. 
!> \param lupri 					Logical print unit
SUBROUTINE pbc_spectral_decomp_ovl(Sabk,U,is_singular,Ndim,nsingular,& 
		& singular_threshh,lupri)
	IMPLICIT NONE
	! input
	INTEGER,INTENT(IN) :: Ndim,lupri
	INTEGER,INTENT(INOUT) :: nsingular
	COMPLEX(complexk),INTENT(IN) :: sabk(Ndim,Ndim)
	COMPLEX(complexk),INTENT(INOUT),pointer :: U(:,:)
	REAL(realk), INTENT(IN) :: singular_threshh
	LOGICAL,INTENT(INOUT) :: is_singular
	! local
	REAL(realk),POINTER :: rwork(:),w(:)
	REAL(realk) :: wtemp
	COMPLEX(complexk) ,POINTER :: Work(:),sigma(:,:)
	COMPLEX(complexk), POINTER :: Sk(:,:),diag(:)
	COMPLEX(complexk),POINTER  :: sabk_tmp(:,:),sabk2(:,:)
	COMPLEX(complexk) :: alpha,beta
	INTEGER :: info,i,j,lwork,nonsingdim
        INTEGER,save :: kp=0
        kp=kp+1

	nsingular = 0
	lwork=10*Ndim-1
	is_singular = .false.

	call mem_alloc(work,max(1,lwork))
	call mem_alloc(rwork,max(1,3*Ndim-2))
	call mem_alloc(w,ndim)
	call mem_alloc(Sk,ndim,ndim)
	call mem_alloc(diag,ndim)

	alpha=CMPLX(1.0_realk,0.0_realk,complexk)
	beta=CMPLX(0.0_realk,0.0_realk,complexk)

	Sk(:,:)=sabk(:,:)

       ! if(kp==1) then
       !   write(*,*) 'S gamma ='
       !   call write_zmatrix(sabk,ndim,ndim)
       ! endif

	diag(:)=CMPLX(0.0_realk,0.0_realk,complexk)
	! calculate spectral decomposition S = V^H\sigma V  
	call zheev('V','U',Ndim,Sk,Ndim,w,work,Lwork,Rwork,info)

	do i=1,ndim
		if(w(i) .lt. 0._realk) then
			write(*,*) 'Error: Eigenvalue of overlap matrix is negative. &
				& Eigv(',i,') =', w(i)
			write(lupri,*) 'Error: Eigenvalue of overlap matrix is negative. &
				& Eigv(',i,') =', w(i)
		endif
	end do

	call mem_dealloc(rwork)
	call mem_dealloc(work)

	! reorder elements s.t. largest eigenv comes first
	do i=1,(ndim+1)/2 
		j=ndim-i+1
		!swap eigenvectors
		diag(:)=Sk(:,i)
		Sk(:,i)=Sk(:,j)
		Sk(:,j)=diag(:)
		!swap eigenvalues
		wtemp=w(i)
		w(i)=w(j)
		w(j)=wtemp
	end do


	! remove linear dependencies in v (sk)
	do i=1,ndim
		if(w(i) .lt. singular_threshh) then
			nsingular=nsingular+1
			sk(:,i)=cmplx(0.0_realk,0.0_realk,complexk) 
			diag(i)=cmplx(0.0_realk,0.0_realk,complexk)
			is_singular=.true.
		else
			diag(i)=cmplx(1.0_realk/sqrt(w(i)),0.0_realk,complexk)
		endif
                !if(kp==1) write(*,*) w(i),'big sigma'
	enddo

	write(*,*) 'number of singulars',nsingular
	write(lupri,*) 'number of singulars',nsingular

	! calculate u = v\sigma
	nonsingdim=ndim-nsingular
	call mem_alloc(u,ndim,nonsingdim)
	do i = 1, nonsingdim
		u(:,i) = sk(:,i)*diag(i)
	end do

	if(nsingular .gt. 0) then
		write(*,*) 'Number of singulars',nsingular,&
			&'lowest eigenvalue of S(k)',w(ndim-nsingular)
		write(lupri,*) 'Number of singulars of S(k)',nsingular,&
			&'lowest eigenvalue',w(ndim-nsingular)
	endif
	call mem_dealloc(diag)

	! Test the difference between the original S matrix and the one with the
	! singularities removed. todo Is this test necc?
	if (is_singular) then

		call mem_alloc(sigma,ndim,ndim)
		call mem_alloc(sabk_tmp,ndim,ndim)
		call mem_alloc(sabk2,ndim,ndim)
		!sabk_tmp(:,:)=CMPLX(0._realk,0._realk,complexk)
		!sabk2(:,:)=CMPLX(0._realk,0._realk,complexk)
		sigma(:,:)=CMPLX(0._realk,0._realk,complexk)
		do i=1,nonsingdim
			sigma(i,i)=cmplx(w(i),0._realk,complexk)
		enddo

		call zgemm('N','N',ndim,ndim,ndim,alpha,sk,ndim,&
			& sigma,ndim,beta,sabk_tmp,ndim)
		call zgemm('N','C',ndim,ndim,ndim,alpha,sabk_tmp,ndim,&
			& sk,ndim,beta,sabk2,ndim)

		write(*,*) wtemp**2*singular_threshh,'should be max error in S matrix'
		do i=1,ndim
			do j=1,ndim
				beta=sabk(i,j) - sabk2(i,j)
				if(abs(beta) .gt. wtemp**2*singular_threshh) then
					write(*,*) 'S-Sprime',sabk(i,j)-sabk2(i,j) 
					write(*,*) 'S-Sprime,singular,i,j',is_singular,i,j
				endif
			enddo
		enddo
		call mem_dealloc(sigma)
		call mem_dealloc(sabk_tmp)
		call mem_dealloc(sabk2)
	endif

	call mem_dealloc(w)
	call mem_dealloc(sk)

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
  INTEGER,INTENT(in) :: ndim,nonsingdim
  COMPLEX(complexk),INTENT(IN) :: kfock(Ndim,Ndim),Uk(Ndim,nonsingdim)
  COMPLEX(complexk),INTENT(INOUT) :: tfock(Nonsingdim,Nonsingdim)
  ! local
  COMPLEX(complexk),POINTER :: tmp(:,:)
  COMPLEX(complexk) :: alpha,beta
  INTEGER :: i,j

  alpha=CMPLX(1._realk,0._realk,complexk)
  beta=CMPLX(0._realk,0._realk,complexk)
  call mem_alloc(tmp,ndim,nonsingdim)
  call zgemm('N','N',Ndim,Nonsingdim,Ndim,alpha,kfock,Ndim,&
             & Uk,Ndim,beta,tmp,Ndim)
  call zgemm('C','N',Nonsingdim,Nonsingdim,Ndim,alpha,Uk,Ndim,&
             & tmp,Ndim,beta,tfock,Nonsingdim)
  call mem_dealloc(tmp)

END SUBROUTINE pbc_unitary_transform

!> \author JR
!> \date 2013
!> \brief Transform Fock matrix to MO basis. F = C^HFC.
!> \param Ctmp 	HF exp. coeff.
!> \param fock 	Fock matrix in AO basis.
!> \param fockMO 	Fock matrix in MO basis.
!> \param ndim 	Mtrx. dim.
!> \param lupri 	Logical print unit.
SUBROUTINE transform_toMOfock(C_tmp,fock,fockMO,ndim,lupri)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: ndim,lupri
	COMPLEX(complexk),INTENT(IN) :: C_tmp(ndim,ndim),fock(ndim,ndim)  
	COMPLEX(complexk),INTENT(OUT) :: fockMO(ndim,ndim)
	! local 
	COMPLEX(complexk),pointer :: fockMOtmp(:,:)
	COMPLEX(complexk) :: alpha,beta

	alpha=CMPLX(1.0_realk,0.0_realk,complexk)
	beta =CMPLX(0.0_realk,0.0_realk,complexk)
	call mem_alloc(fockMotmp,ndim,ndim)
	fockMO=cmplx(0.0_realk,0.0_realk,complexk)
	fockMOtmp=cmplx(0.0_realk,0.0_realk,complexk)
	call zgemm('C','N',ndim,ndim,ndim,alpha,C_tmp,ndim,fock,ndim,&
		& beta,fockMOtmp,ndim)
	call zgemm('N','N',ndim,ndim,ndim,alpha,fockMOtmp,ndim,C_tmp,ndim,&
		& beta,fockMO,ndim)

	call mem_dealloc(fockmotmp)

END SUBROUTINE transform_toMOfock

!> \author JR 
!> \date 2013
!> \brief 
!> \param error
!> \param fockMO
!> \param ndim
!> \param errlm
!> \param nelectrons
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
			!write(*,*) error(k)
		ENDDO
	ENDDO

END SUBROUTINE  pbc_geterrorvec

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
        integer, save :: kp=0
        kp=kp+1

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
	!write(*,*) 'tfock'
	!call write_zmatrix(tfock,nonsingdim,nonsingdim) !todo necc?
	! solve fc = ce 
	call zheev('V','U',nonsingdim,tfock,nonsingdim,eigv,work,lwork,rwork,info)
	if(info .ne. 0) then
		write(lupri,*) 'ERROR: zheev problems, info=', info,kp
		call write_zmatrix(tfock,nonsingdim,nonsingdim)
		write(*,*) 'ERROR: zheev problems, info=', info,kp
		!call write_zmatrix(uk,nonsingdim,nonsingdim)
		write(lupri,*) 'ERROR: zheev problems, info=', info,kp
		call LSQUIT('pbc_solve_kfcsc_mat: INFO not zero,', &
			& ' while solving eigenvalue',lupri)
	endif
	!Transform back to C=Uc
	call zgemm('N','N',ndim,nonsingdim,nonsingdim,alpha,Uk,ndim,&
		& tfock,nonsingdim,beta,C_tmp(:,1:nonsingdim),ndim)
	c_tmp(:,nonsingdim+1:ndim)=cmplx(0._realk,0._realk,complexk)

	call mem_dealloc(work)
	call mem_dealloc(rwork)
	call mem_dealloc(tfock)

END SUBROUTINE solve_kfcsc_mat

!> \author JR 
!> \date 2013
!> \brief Does the SCF iterations and solves the F(k)C(k)=eps(k)C(k) for each k.
!> \param molecule 		Stores info about the reference cell.
!> \param setting 		Stores info about the integral routines. 
!> \param lattice 		Stores info about the pbc lattice.
!> \param numrealvec 	Number of unit vectors.
!> \param maxmultmom
!> \param bz 				Stores info about the reciprocal grid.
!> \param dmat0 			The initial density matrix.
!> \param lupri 			Logical print unit.
!> \param luerr 			
SUBROUTINE pbc_startzdiis(molecule,setting,ndim,lattice,numrealvec,&
		& maxmultmom,bz,dmat0,lupri,luerr)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: ndim,lupri,luerr,numrealvec,maxmultmom
	! input
	TYPE(lvec_list_t),INTENT(INOUT) :: lattice
	TYPE(moleculeinfo),INTENT(INOUT) :: molecule
	TYPE(LSSETTING) :: setting
	TYPE(BZgrid_t),INTENT(INOUT) :: bz
	TYPE(matrix),INTENT(IN) :: dmat0
	! local
	COMPLEX(COMPLEXK),POINTER :: fockMO(:,:),fock(:,:),smatk(:,:)
	COMPLEX(COMPLEXK),POINTER :: C_k(:,:),D_k(:,:),C_0(:,:)
	REAL(realk),POINTER :: cellenergies(:)
	REAL(realk),POINTER :: error(:,:),nucmom(:)
	INTEGER :: i,j,errlm,tol!tol should be a real but now I just have it for test
	INTEGER :: k,kpt,n1,fdim(3),layer
	INTEGER :: realcut(3),nonsingdim
	real(realk) :: kvec(3),Ecell,E_nn,E_1,errortest
	REAL(realk),POINTER :: tlat(:,:),weight(:)
	TYPE(matrix), POINTER :: nfdensity(:)
	TYPE(moleculeinfo) :: refcell
	CHARACTER(LEN=10) :: stiter
	CHARACTER(LEN=12) :: diis,diismats
	CHARACTER(LEN=20) :: mattxt
	LOGICAL :: diis_exit
	REAL(realk) :: E_J,E_K,E_XC,E_ff,E_cell
	REAL(realk) :: E_nuc,E_nnff
	REAL(realk)         :: TS,TE,TST,TET,TOT,TWT !For finding time usage
	TYPE(matrix),POINTER :: f_1(:),ovl(:)
	TYPE(matrix),POINTER :: g_2(:)

	! threshhold for removing singularities in overlap matrix s
	REAL(realk) :: singular_threshh = 1e-4_realk 
	
	! timings
	integer timer_exc, timer_coul, timer_fmm, timer_h1, timer_ovl 
	integer timer_fmm_mm, timer_kptenergy

	type(timing_info), pointer :: timings(:)
	timer_exc = 1
	timer_coul = 2 
	timer_fmm = 3 
	timer_h1 = 4
	timer_ovl = 5
	timer_fmm_mm = 6
	timer_kptenergy = 7

	! initialize timings module
	call pbc_timings_init(timings, 7)
	call pbc_timings_name_timer(timings, 'Exchange integrals  ', timer_exc)
	call pbc_timings_name_timer(timings, 'Coulomb integrals   ', timer_coul)
	call pbc_timings_name_timer(timings, 'Farfield integrals  ', timer_fmm)
	call pbc_timings_name_timer(timings, '1 particle integrals', timer_h1)
	call pbc_timings_name_timer(timings, 'Overlap Matrices    ', timer_ovl)
	call pbc_timings_name_timer(timings, 'FMM-multip. moments ', timer_fmm_mm)
	call pbc_timings_name_timer(timings, 'K-point energy      ',timer_kptenergy)

	write(lupri,*) 'Entering routine startzdiis'

	write(stiter,'(I5)') 1
	stiter=adjustl(stiter)
	diis='diis_'//trim(stiter)//'_'

	call mem_alloc(tlat,(maxmultmom+1)**2,(maxmultmom+1)**2)
	call mem_alloc(nfdensity,numrealvec)
	call find_latt_index(n1,0,0,0,lattice,lattice%max_layer)
	call mat_init(nfdensity(n1),ndim,ndim)
	lattice%lvec(n1)%dm_computed=.true.
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

	write(lupri,*) 'Density first'
	call mat_print(nfdensity(n1),1,ndim,1,ndim,lupri)

	call set_refcell(refcell,molecule)
	call set_lattice_cells(numrealvec,molecule,lattice,lupri)
	errlm=molecule%nelectrons/2*(ndim-molecule%nelectrons/2.)
	call mem_alloc(error,lattice%num_store,errlm)

	call pbc_timings_start(timings, timer_ovl)
	call pbc_overlap_k(lupri,luerr,setting,molecule%natoms,ndim,&
		& lattice,refcell,numrealvec,ovl)
	call pbc_timings_stop(timings, timer_ovl, 2, lupri)

	!calculate 1 particle operators
	call pbc_timings_start(timings, timer_h1)
	
	! - kinetic energy of electrons
	call pbc_kinetic_k(lupri,luerr,setting,molecule%natoms,ndim,&
		& lattice,refcell,numrealvec,nfdensity,f_1)

	! - electron nuclei attraction
	call pbc_nucattrc_k(lupri,luerr,setting,molecule%natoms,ndim,&
		& lattice,refcell,numrealvec,nfdensity,f_1)

	! - nuclear repulsion
	CALL pbc_nucpot(lupri,luerr,setting,molecule%natoms,lattice,&
		& refcell,numrealvec,E_nuc)

	call pbc_timings_stop(timings, timer_h1, 2, lupri)

	!computing the nuclear moments
	call mem_alloc(nucmom,(1+maxmultmom)**2)
	call pbc_comp_nucmom(refcell,nucmom,maxmultmom,lupri)

	call pbc_timings_start(timings, timer_fmm_mm)
	call pbc_multipole_expan_k(lupri,luerr,setting,ndim,lattice,&
		&refcell,numrealvec,maxmultmom)
	call pbc_timings_stop(timings, timer_fmm_mm, 2, lupri)

	call pbc_controlmm(20,Tlat,lattice%Tlmax,maxmultmom,.false.,lattice%ldef%avec,&
		& ndim,lupri,nfdensity,numrealvec,lattice,E_ff,E_nnff,refcell)

	do k=1,bz%nk
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
		call pbc_zdevectorize_mat(smatk,ndim,ndim,bz%smat%zelms)
		!diagonalizes Sk, Uk transform operators
		call pbc_spectral_decomp_ovl(smatk,bz%kpnt(k)%Uk, & 
			& bz%kpnt(k)%is_singular,Ndim,bz%kpnt(k)%nsingular,singular_threshh,lupri)
		call mem_alloc(bz%kpnt(k)%eigv,ndim-bz%kpnt(k)%nsingular)
	enddo

	k=0
	i=0 !either 0 or 1, check it, for zero no errors
	tol=0

	Ecell=0.0_realk
	diis_exit = .false. !when diis_exit the iterations are finished
	! self consistent iterations
	DO WHILE(tol .le. lattice%num_its)! 20)!should have an input parameter

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

		call pbc_get_onehamenergy(numrealvec,lattice,f_1,nfdensity,E_1)

		! coulomb integrals
		call pbc_timings_start(timings, timer_coul)
		call pbc_electron_rep_k(lupri,luerr,setting,molecule%natoms,ndim,&
			& lattice,refcell,numrealvec,nfdensity,g_2,E_J)
		call pbc_timings_stop(timings, timer_coul, 2, lupri)

		! exchange integrals
		call pbc_timings_start(timings, timer_exc)
		call pbc_exact_xc_k(lupri,luerr,setting,molecule%natoms,ndim,&
			& lattice,refcell,numrealvec,nfdensity,g_2,E_K)
		call pbc_timings_stop(timings, timer_exc, 2, lupri)

		lattice%fc1=max(lattice%oneop1,lattice%col1)
		lattice%fc1=max(lattice%fc1,lattice%Kx1)
		lattice%fc2=max(lattice%oneop2,lattice%col2)
		lattice%fc2=max(lattice%fc2,lattice%Kx2)
		lattice%fc3=max(lattice%oneop3,lattice%col3)
		lattice%fc3=max(lattice%fc3,lattice%Kx3)

		!Far-field contribution to fock
		call pbc_timings_start(timings, timer_fmm)
		call pbc_ff_fck(maxmultmom,tlat,lattice%tlmax,ndim,lattice,nfdensity,nucmom,&
			& g_2,E_ff,E_nn,lupri)
		call pbc_timings_stop(timings, timer_fmm, 2, lupri)

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
		call pbc_timings_start(timings, timer_kptenergy)
		do kpt=1,bz%nk

			call pbc_get_kpoint(kpt,kvec)
			call zero_pbc_elstr(Bz%fck)
			call zero_pbc_elstr(Bz%Smat)

			! We compute only weights for the gamma point
			if(bz%kpnt(kpt)%is_gamma )then
				call mem_alloc(weight,i)
				weight=0.0_realk
				call pbc_get_diisweights(lattice,Bz,weight,i,tol,kvec,ndim,C_0,fockMO,fock,errortest,error,&
					& diis_exit,errlm,molecule%nelectrons,lupri)
			endif ! is_gamma

			if(tol .gt. 0) then
				if(.not. diis_exit) then
					!CALL LSTIMER('START',TS,TE,LUPRI)
					call pbc_get_weighted_fock(i,tol,lattice%num_store,ndim,weight,lattice,lupri)
					!CALL LSTIMER('Weighted fock',TS,TE,LUPRI)
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
			call solve_kfcsc_mat(ndim,fock,C_k, &
				& bz%kpnt(kpt)%Uk,bz%kpnt(kpt)%eigv, &
				& bz%kpnt(kpt)%nsingular,lupri)

			!C_0 is used for finding the weights
			if(bz%kpnt(kpt)%is_gamma ) C_0(:,:)=C_k(:,:)
			!call pbc_zdevectorize_mat(smatk,ndim,ndim,bz%smat%zelms)

			!gets D(k) from C(k)
			call pbc_get_kdensity(D_k,C_k,ndim,molecule%nelectrons/2 &
				& ,bz%kpnt(kpt)%nsingular,smatk,lupri)

			!Converts D(k) to D^0l
			call kspc_2_rspc_loop_k(nfdensity,Bz%Nk,D_k,lattice,kvec,bz%kpnt(kpt)%weight,&
                          & BZ%NK_nosym,ndim,kpt,diis_exit,lupri)

		enddo !kpt
		!CALL LSTIMER('k point energy',TST,TET,LUPRI)
		call pbc_timings_stop(timings, timer_kptenergy, 2, lupri)

                if(diis_exit) then
                  write(lupri,*) 'Max Density elements for each layer'
                  call print_maxdens(nfdensity,lattice,lupri)
                endif

		if(associated(weight)) call mem_dealloc(weight)

		write(*,*)
		call pbc_densitymat_write(nfdensity,lattice,ndim,ndim,8,'            ')
		call pbc_free_read_matrices(lattice)
		call print_bands(bz,ndim,'band-energy') !prints band energy to file band-energy
		! Get HOMO LUMO energy, change name 
		call pbc_trans_k_energy(lattice,cellenergies,ndim,molecule%nelectrons,&
			& bz)

		E_cell=E_1+E_j+E_K+E_ff+E_nuc!+E_nn
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

	ENDDO

	do k=1,bz%Nk
		call mem_dealloc(bz%kpnt(k)%Uk)
	enddo

	do i=1,numrealvec
		call free_Moleculeinfo(lattice%lvec(i)%molecule)
		if(nfdensity(i)%init_magic_tag.NE.mat_init_magic_value) CYCLE
		call mat_free(nfdensity(i))
	enddo
		
	! print timings
	call pbc_timings_print(timings, lupri)
	call pbc_timings_print(timings, 6)
	! deallocate allocated memory in timings 
	call pbc_timings_destruct(timings)

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

	write(lupri,*) 'Final E(HOMO) =', cellenergies(1)
	write(lupri,*) 'Final E(LUMO) =', cellenergies(2)
	write(lupri,*) 'Final Cell Energy =', E_cell
	write(lupri,*) 'Final h_1=',E_1
	write(lupri,*) 'Final Nuclear=',E_nuc
	write(lupri,*) 'Final Far field=', E_ff,E_nn
	write(lupri,*) 'Final K energy', E_k
	write(lupri,*) 'Final J energy', E_J
	write(*,*) 'Final E(HOMO) =', cellenergies(1)
	write(*,*) 'Final E(LUMO) =', cellenergies(2)
	write(*,*) 'Final Cell Energy =', E_cell
	write(*,*) 'Final K energy', E_k
	write(*,*) 'Final J energy', E_J
	write(*,*) 'Final h_1=',E_1
	write(*,*) 'Final Nuclear=',E_nuc
	write(*,*) 'Final Far field=', E_ff,E_nn
	call mem_dealloc(cellenergies)

END SUBROUTINE pbc_startzdiis

!> \author JR 
!> \date 2013
!> \brief 
!> \param lattice
!> \param g_2
!> \param f_1
!> \param ndim
!> \param realcut
!> \param numrealvec
!> \param diismats
!> \param lupri
SUBROUTINE pbc_get_fock_mat(lattice,g_2,f_1,ndim,realcut,numrealvec,diismats,lupri)
	IMPLICIT NONE
	TYPE(lvec_list_t),INTENT(INOUT) :: lattice
	INTEGER,INTENT(IN) :: ndim,lupri,numrealvec
	TYPE(matrix),target,intent(inout) :: f_1(numrealvec)
	TYPE(matrix),target,intent(inout) :: g_2(numrealvec)
	CHARACTER(LEN=12) :: diismats
	INTEGER,intent(INOUT) :: realcut(3)
	! local
	INTEGER :: i,j
	REAL(realk) :: focknorm !for finding time usage
	REAL(realk) :: ts,te !for finding time usage

	!get the fock matrices f^0l
	if(lattice%store_mats) then
		call lstimer('START',ts,te,lupri)
		!call pbc_read_fock_matrix(lattice,ndim,ndim,diismats)
		call pbc_read_fock_matrix(lattice,ndim,ndim,'            ')
		call lstimer('Reading fock',ts,te,lupri)
		call pbc_fockmat_write(lattice,ndim,ndim,7,2,diismats,lupri)
		!get the overlap matrices S^0l
		call pbc_read_matrix(lattice,ndim,ndim,1,1,'            ')
	else
		call pbc_add_fock_matrix(f_1,g_2,lattice,ndim,ndim,numrealvec)
		!call pbc_fockmat_write(lattice,ndim,ndim,7,2,diismats,lupri)
	endif
END SUBROUTINE pbc_get_fock_mat



!> \author JR 
!> \date 2013
!> \brief Calc the electronic one body energy. 
!> \param numvecs
!> \param f_1
!> \param nfdensity
!> \param E1
SUBROUTINE  pbc_get_onehamenergy(numvecs,lattice,f_1,nfdensity,E_1)
	INTEGER,INTENT(IN) :: numvecs
	TYPE(matrix),TARGET,INTENT(IN) :: f_1(numvecs),nfdensity(numvecs)
	TYPE(lvec_list_t),INTENT(INOUT) :: lattice
	REAL(realk),INTENT(INOUT) ::  E_1
	!LOCAL
	INTEGER :: celli,l1,l2,l3

	E_1=0
	do celli=1,numvecs
		if( lattice%lvec(celli)%f1_computed &
			& .and. &
			& nfdensity(celli)%init_magic_tag .eq. mat_init_magic_value ) &
			& then
                      if(.not.lattice%lvec(celli)%is_redundant)then
                        l1=int(lattice%lvec(celli)%lat_coord(1))
                        l2=int(lattice%lvec(celli)%lat_coord(2))
                        l3=int(lattice%lvec(celli)%lat_coord(3))
		        if(l1**2+l2**2+l3**2 .gt. 0) then
                          E_1=E_1+2.0_realk*mat_dotproduct(f_1(celli),nfdensity(celli))
                        else
                          E_1=E_1+mat_dotproduct(f_1(celli),nfdensity(celli))
                        endif
                      endif
                    endif
	enddo

END SUBROUTINE  pbc_get_onehamenergy

!> \author JR 
!> \date 2013
!> \brief Find the density matrix in resiprocal space. Only for unres calc where
!> \brief alpha and beta densities are similar. 
!> \param ddensity
!> \param c_tmp
!> \param nbast
!> \param nkmobas
!> \param nsingular
!> \param smathk
!> \param lupri
SUBROUTINE pbc_get_kdensity(ddensity,C_tmp,nbast,nkmobas,nsingular,smatk,lupri)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: nbast,nkmobas,lupri,nsingular
	COMPLEX(complexk),INTENT(INOUT) :: ddensity(nbast,nbast)
	COMPLEX(complexk),INTENT(IN) :: C_tmp(nbast,nbast),smatk(nbast,nbast)
	! local
	REAL(realk) :: nelectrons,dummy1,dummy2
	COMPLEX(complexk), POINTER :: density_tmp(:,:)
	COMPLEX(complexk), POINTER :: tmp(:,:)
	INTEGER :: i,j,nosingdim
	COMPLEX(complexk) :: alpha,beta

	alpha=CMPLX(2.0_realk,0.0_realk,complexk)
	beta =CMPLX(0.0_realk,0.0_realk,complexk)

	nosingdim=nbast-nsingular
	ddensity(:,:)=CMPLX(0.0_realk,0.0_realk,complexk)

	call zgemm('N','C',nbast,nbast,nkmobas,alpha,c_tmp,nbast,&
		& c_tmp,nbast,beta,ddensity,nbast)

  !=========================================
  ! DEBUG below (smathk,tmp, only used here)
  ! test number of electrons etc
  !=========================================
	
	call mem_alloc(tmp,nbast,nbast)
  
	alpha=CMPLX(1.0_realk,0.0_realk,complexk)
	call zgemm('N','N',nbast,nbast,nbast,alpha,ddensity,nbast,&
		& smatk,nbast,beta,tmp,nbast)

	nelectrons=0._realk
	do i=1,nbast
		nelectrons=nelectrons+real(tmp(i,i))
	enddo

	!write(*,*) 'Nelectrons =', nelectrons,nkmobas,nsingular
	write(lupri,*) 'Nelectrons =', nelectrons,nkmobas,nsingular

	dummy2=-huge(dummy2)
	do i=1,nbast
		do j=1,nbast
			dummy1=ddensity(i,j)
			dummy2=max(dummy1,dummy2)
		enddo
	enddo

	!write(*,*) 'max value of D(k)',dummy2,'singularities', nsingular
	write(lupri,*) 'max value of D(k)',dummy2,'singularities', nsingular

	call mem_dealloc(tmp)

END SUBROUTINE pbc_get_kdensity

!> \author JR 
!> \date 2013
!> \brief ??? 
!> \param lattice 		Information about the pbc lattice.
!> \param Bz 				Information about the reciprocal lattice.
!> \param weight 			Diis weights.
!> \param its 				???
!> \param tol 				???
!> \param kvec 			???
!> \param ndim 			Fock matrix dim
!> \param C_0 				Fock exp. coeff. (??) basis.
!> \param fockMO 			Fock mat in MO basos
!> \param fock 			Fock mat in AO basis
!> \param errortest 		???
!> \param error 			???
!> \param diis_Exit 		???
!> \param errdim 			???
!> \param nelectrons 	Number of electrons total.
!> \param lupri 			Logical print unit.
SUBROUTINE pbc_get_diisweights(lattice,Bz,weight,its,tol,kvec,ndim,C_0,fockMO,fock,&
                & errortest,error,diis_exit,errdim,nelectrons,lupri)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim,lupri,its,tol,errdim
  INTEGER,INTENT(IN) :: nelectrons
  TYPE(lvec_list_t),INTENT(INOUT) :: lattice
  REAL(realk),INTENT(OUT) :: weight(its)
  REAL(realk),INTENT(INOUT) :: error(lattice%num_store,errdim)
  LOGICAL,INTENT(OUT) :: diis_exit
  REAL(realk),INTENT(OUT) :: errortest
  REAL(realk),INTENT(IN) :: kvec(3)
  TYPE(BZgrid_t),INTENT(INOUT) :: bz
  COMPLEX(complexk) :: C_0(ndim,ndim)
  COMPLEX(complexk),INTENT(OUT) :: fockMO(ndim,ndim)
  COMPLEX(complexk),INTENT(INOUT) :: fock(ndim,ndim)
  !LOCAL
  INTEGER :: j

      !We need k-space fock matrix in gamma point now
      call pbc_rspc_to_kspc_mat(lattice,Bz,ndim,kvec,2)
      call pbc_zdevectorize_mat(fock,ndim,ndim,bz%fck%zelms)

        if(tol .ge. 1) then
          call transform_toMOfock(C_0,fock,fockMO(:,:),ndim,lupri)

          weight(:)=0.0_realk
          ! k is still the gamma point
          if(tol .gt. lattice%num_store) then
            do j=1,lattice%num_store-1
            error(j,:)=error(j+1,:)
            enddo
            error(lattice%num_store,:)=0.0_realk
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
          write(lupri,*) 'check error', errortest
          write(*,*) 'check error', errortest
          if(errortest .le. lattice%error) diis_exit=.true.

          !Get diis weights TODO !!!THIS HAS TO BE FIXED, its-1 only when its .gt. tol
          if(its .gt. tol) then ! this since I do not know the C0 matrix for it 0
            call pbc_diisweights(errdim,error,weight,its-1,lattice%num_store,lupri)
          else
            call pbc_diisweights(errdim,error,weight,its,lattice%num_store,lupri)
          endif
        endif

        if(tol .le. 1) weight(1)=1.0_realk

        write(*,*) 
        write(*,*) 'Iteration nr. ', tol
        write(lupri,*) 
        write(lupri,*) 'Iteration nr. ', tol
        write(lupri,*) 'Weights'
        write(lupri,*) weight

        call zero_pbc_elstr(Bz%fck)
        call zero_pbc_elstr(Bz%Smat)

END SUBROUTINE pbc_get_diisweights

!> \author JR 
!> \date 2013
!> \brief 				???
!> \param errdim 		???
!> \param weight 		???
!> \param it 			???
!> \param num_store	???
!> \param lupri 		Logical print unit.
SUBROUTINE pbc_diisweights(errdim,error,weight,it,num_store,lupri)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: errdim,lupri,it,num_store
	REAL(realk),INTENT(INOUT) :: error(num_store,errdim)
	REAL(realk),INTENT(INOUT) :: weight(it)
	! local
	REAL(realk) :: B_mat(it,it),weight_tmp(it,1)
	REAL(realk) :: Sdgelss(it+1),rcond,normfac
	INTEGER :: i,j,info,N,m,rank,lwork
	INTEGER :: solve(it)
	REAL(realk),pointer :: work(:)

	B_mat(:,:)=0.0_realk
	DO i=1,it
		DO j=1,it
			B_mat(i,j)=dot_product(error(i,:),error(j,:))! error(i)*error(j)
		ENDDO
	ENDDO

	solve=0
	weight_tmp(:,1)=1.0_realk

	N=it
	m=1
	call dgesv(N,m,B_mat,N,solve,weight_tmp,N,info)
	IF(info .ne. 0) THEN
		write(*,*) 'Calls dgelss instead'
		write(lupri,*) 'Calls dgelss instead'
		DO i=1,it
			DO j=1,it
				B_mat(i,j)=dot_product(error(i,:),error(j,:))! error(i)*error(j)
				B_mat(j,i)=B_mat(i,j)
			ENDDO
		ENDDO
		DO j=1,it
			weight_tmp(j,1)=1.0_realk
		ENDDO
		lwork=3*n+2*n+1000
		rcond= 1.0e-10_realk
		call mem_alloc(work,lwork)
		call dgelss(n,n,m,B_mat,n,weight_tmp,n,Sdgelss,rcond,rank,work,lwork,info)
		if(info .ne. 0) then
			call LSQUIT('pbc_diisweights: INFO not zero, while solving eigenvalue',&
				& lupri)
		endif
		call mem_dealloc(work)
	endif

	DO i=1,it
		weight(i)=weight_tmp(i,1)
	ENDDO

	normfac=sum(weight(1:it))
	weight(1:it)=weight(1:it)/normfac

END SUBROUTINE pbc_diisweights

!> \author JR 
!> \date 2013
!> \brief Sums over former Fock matrices with the corresponding weights. 
!> \param i 			???
!> \param ndiis 		???
!> \param stdiis 		???
!> \param nbast 		Number of basis functions.
!> \param weight 		???
!> \param aop 			???
!> \param lupri 		LOgical print unit
SUBROUTINE pbc_get_weighted_fock(i,ndiis,stdiis,nbast,weight,aop,lupri)
	IMPLICIT NONE
	! input
	INTEGER,INTENT(IN) :: i,ndiis,stdiis,nbast,lupri
	REAL(realk),INTENT(IN) :: weight(i)
	TYPE(lvec_list_t),INTENT(INOUT) :: aop
	! local
	INTEGER :: j,layer,k
	INTEGER :: l1,l2,l3
	CHARACTER(LEN=12) :: diis,tmpdiis
	CHARACTER(LEN=7) :: stiter
	TYPE(lvec_list_t) :: tmp_mat

	tmp_mat%max_layer=aop%max_layer
	tmp_mat%nneighbour=aop%nneighbour
	tmp_mat%ldef%is_active(:)=aop%ldef%is_active(:)
	call build_lvec_list(tmp_mat,nbast) 
	tmp_mat%fc1=aop%fc1
	tmp_mat%fc2=aop%fc2
	tmp_mat%fc3=aop%fc3
	tmp_mat%oneop1 = aop%oneop1
	tmp_mat%oneop2 = aop%oneop2
	tmp_mat%oneop3 = aop%oneop3
	tmp_mat%col1 = aop%col1
	tmp_mat%col2 = aop%col2
	tmp_mat%col3 = aop%col3
	tmp_mat%kx1 = aop%kx1
	tmp_mat%kx2 = aop%kx2
	tmp_mat%kx3 = aop%kx3

	do layer=1,size(aop%lvec)
		call mat_init(tmp_mat%lvec(layer)%oper(2),nbast,nbast)
		call mat_zero(tmp_mat%lvec(layer)%oper(2))
		tmp_mat%lvec(layer)%g2_computed=aop%lvec(layer)%g2_computed
		tmp_mat%lvec(layer)%f1_computed=aop%lvec(layer)%f1_computed
		tmp_mat%lvec(layer)%ovl_computed=aop%lvec(layer)%ovl_computed
	enddo

	if(ndiis .ge. 1) then
		write(stiter,'(i5)') ndiis
		stiter=adjustl(stiter)
		diis='diis_'//trim(stiter)//'_'
	endif

	do layer=1,size(aop%lvec)
		l1=int(aop%lvec(layer)%lat_coord(1))
		l2=int(aop%lvec(layer)%lat_coord(2))
		l3=int(aop%lvec(layer)%lat_coord(3))
		if(aop%lvec(layer)%f1_computed .or. aop%lvec(layer)%g2_computed)then
			call mat_zero(aop%lvec(layer)%oper(2))
		endif
	enddo

	if(ndiis .le. stdiis) then
		do j=1,ndiis
			write(stiter,'(i5)') j
			stiter=adjustl(stiter)
			tmpdiis='diis_'//trim(stiter)//'_'
			!    write(*,*) tmpdiis

			call pbc_read_matrix(tmp_mat,nbast,nbast,7,2,tmpdiis)

			do layer=1,size(aop%lvec)
				l1=int(aop%lvec(layer)%lat_coord(1))
				l2=int(aop%lvec(layer)%lat_coord(2))
				l3=int(aop%lvec(layer)%lat_coord(3))
				if(aop%lvec(layer)%f1_computed .or. aop%lvec(layer)%g2_computed)then
					call mat_daxpy(weight(j),tmp_mat%lvec(layer)%oper(2),&
						& aop%lvec(layer)%oper(2))
				endif
			enddo
		enddo
	else
		k=0
		do j= ndiis-stdiis+1,ndiis
			k=k+1
			write(stiter,'(i5)') j
			stiter=adjustl(stiter)
			tmpdiis='diis_'//trim(stiter)//'_'

			call pbc_read_matrix(tmp_mat,nbast,nbast,7,2,tmpdiis)

			do layer=1,size(aop%lvec)
				l1=int(aop%lvec(layer)%lat_coord(1))
				l2=int(aop%lvec(layer)%lat_coord(2))
				l3=int(aop%lvec(layer)%lat_coord(3))
				if(aop%lvec(layer)%f1_computed .or. aop%lvec(layer)%g2_computed)then
					call mat_daxpy(weight(k),tmp_mat%lvec(layer)%oper(2),&
						& aop%lvec(layer)%oper(2))
				endif
			enddo
		enddo
	endif

	do layer=1,size(aop%lvec)
		l1=int(aop%lvec(layer)%lat_coord(1))
		l2=int(aop%lvec(layer)%lat_coord(2))
		l3=int(aop%lvec(layer)%lat_coord(3))
		if(aop%lvec(layer)%f1_computed .or. aop%lvec(layer)%g2_computed)then
			call pbc_get_file_and_write(aop,nbast,nbast,layer,7,2,'            ')
		endif
		call mat_free(tmp_mat%lvec(layer)%oper(2))
	enddo

	call mem_dealloc(tmp_mat%lvec)

END SUBROUTINE pbc_get_weighted_fock

!> \author JR 
!> \date 2013
!> \brief ???
!> \param lattice 		Information about the pbc lattice.
!> \param cenergies 		???
!> \param nbast 			Number of basis functions.
!> \param nelectrons 	Number of electrons.
!> \param bz 				Information about the reciprocal grid.
SUBROUTINE pbc_trans_k_energy(lattice,cenergies,nbast,nelectrons,bz)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nbast,nelectrons
	TYPE(lvec_list_t), INTENT(IN) :: lattice
	TYPE(BZgrid_t), INTENT(IN) :: bz
	REAL(realk), INTENT(INOUT) :: cenergies(nbast)
	! local
	INTEGER :: lattindex(3), kpt
	real(realk) :: ehomo,elumo,etmph1,etmph2,etmpl1,etmpl2

	lattindex(1)=0
	cenergies =0.0_realk
	ehomo=-huge(ehomo)
	elumo= huge(elumo)

	lattindex(1)=0
	lattindex(2)=0
	lattindex(3)=0

	cenergies(1)=bz%kpnt(1)%eigv(nelectrons/2)
	cenergies(2)=bz%kpnt(1)%eigv(nelectrons/2+1)

	do kpt=2,bz%nk
		etmph1=bz%kpnt(kpt)%eigv(nelectrons/2)
		etmph2=bz%kpnt(kpt-1)%eigv(nelectrons/2)
		etmpl1=bz%kpnt(kpt)%eigv(nelectrons/2+1)
		etmpl2=bz%kpnt(kpt-1)%eigv(nelectrons/2+1)
		ehomo=max(etmph1,etmph2)
		elumo=min(etmpl1,etmpl2)
	enddo
	cenergies(1)=max(ehomo,cenergies(1))
	cenergies(2)=min(elumo,cenergies(2))

END SUBROUTINE pbc_trans_k_energy

END MODULE pbc_scfdiis

