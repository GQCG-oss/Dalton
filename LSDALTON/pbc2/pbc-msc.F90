

MODULE PBC_MSC
use files
	use precision
	use typedef
	use matrix_module
	use lattice_type
	use lattice_vectors
	use pbc_matrix_operations
	IMPLICIT NONE

	type latticegeom
		!> indicates whether there are meaningful data in the lattice vector variables
		logical :: latvec_ok           
		!> indicates whether a lattice dimension is active (i.e. the system is to 
		!> be replicated along this dimension)
		logical :: dim_is_active(3)    
		!> latvec(k,:) is the k:th lattice vector
		real(realk) :: latvec(3,3)     
		!> reciprocal lattice vectors
		real(realk) :: reclatvec(3,3)  
		!> the matrix inverse of latvec(:,:)
		real(realk) :: invlatvec(3,3)  
		!> volume of the unit cell
		real(realk) :: cell_volume     
		!> number of k points to sample in BZ
		integer :: num_kpoints   
		!real(realk) :: realspc_thres   !> defines when to neglect Fourier comp.
		!integer :: nfdebug     !> manual NF size for debugging
		integer :: num_k1, num_k2, num_k3
		!> counts the number of updates of Qfict
		integer :: qfict_update_cnt 
		!> fictitious charges
		real(realk) :: qfict(4)           
		!> position (standard coord.)
		real(realk) :: qfict_pos_std(3,4) 
	end type latticegeom

	type(latticegeom), save :: lat_data

	type occ_scheme_typ
		logical insulator_occ
		logical force_symmetric
	end type occ_scheme_typ

	type(occ_scheme_typ), save :: occ_scheme

	type scf_scheme_typ
		logical :: use_tc_dmat
		logical :: use_Cmax2_tol
		integer :: nproject
		real(realk) :: lindep_tol
		real(realk) :: Cmax2_tol
	end type scf_scheme_typ

	type(scf_scheme_typ), save :: scf_scheme

	type lindep_data_t
		logical :: do_projection
		logical :: use_nproject
		logical :: be_silent
		integer :: num_occ
		integer :: nproject
		real(realk) :: Seig_tol
		real(realk) :: Cmax2_tol
	end type lindep_data_t

	! ========================================================
	! data structures to handle the k-point grid
	! ========================================================

	integer, parameter :: Max_bassiz = 1000
	integer, parameter :: Max_kpoints = 1500

	type pbc_elstr_t
		COMPLEX(complexk),pointer :: zelms(:)
		COMPLEX(complexk),pointer :: kfockvec(:)
		COMPLEX(complexk),pointer :: kfockmat(:,:)
		COMPLEX(complexk),pointer :: keigv(:)
		COMPLEX(complexk),pointer :: koverlapvec(:)
		COMPLEX(complexk),pointer :: koverlapmat(:,:)
		COMPLEX(complexk),pointer :: kcdensityvec(:)
		COMPLEX(complexk),pointer :: kcdensitymat(:,:)
		COMPLEX(complexk),pointer :: kddensityvec(:)
		COMPLEX(complexk),pointer :: kddensitymat(:,:)
	end type pbc_elstr_t

	!> This is a data structure to represent a single k-point
	type BZpoint_t
		logical :: self_dual, is_gamma,is_singular
		integer :: ix_orig
		integer :: n(3),nsingular
		integer :: ninv(3)
		real(realk) :: weight
		real(realk) :: lambda(3)
		real(realk), pointer :: eigv(:)
		COMPLEX(complexk),pointer :: Uk(:,:),Uinv(:,:)
	end type BZpoint_t

	!> This is a data structure to represent a sampling grid in
	!> the first Brillouin zone
	type BZgrid_t
		logical :: use_invsym
		integer :: Nk_dim1, Nk_dim2, Nk_dim3, Nk
		integer :: Nk_nosym
		real(realk) :: reclvec(3,3)
		type(BZpoint_t) :: kpnt(Max_kpoints)
		type(pbc_elstr_t) :: fck
		type(pbc_elstr_t) :: Smat
	end type BZgrid_t

	!> This data structure represents coarse-grained grids that
	!> when joined form the full fine-grained grid
	type splitBZgrid_t
		integer :: Nsplit(3)
		type(BZgrid_t), pointer :: subBZ(:,:,:)
		type(BZgrid_t), pointer :: fullBZ
	end type splitBZgrid_t


	type pbc_scfiterations_t
		type(pbc_elstr_t), pointer :: kdep_it(:)
	end type

CONTAINS

!> \author Johannes Rekkedal
!> \date 2013
!> \brief ???
!> \param rlvec ???
!> \param nk1 ???
!> \param nk2 ???
!> \param nk3 ???
!> \param symtxt ???
!> \param BZ ???
!> \param ndim1 ???
!> \param ndim2 ???
!> \param lupri ???
SUBROUTINE pbc_init_BZgrid(rlvec,nk1,nk2,nk3,symtxt,BZ,ndim1,ndim2,lupri)
	implicit none
	! input and output arguments
	integer,intent(in) :: lupri,ndim1,ndim2
	character*(*), intent(in) :: symtxt
	integer, intent(in) :: nk1, nk2, nk3
	real(realk), intent(in) :: rlvec(3,3) !reciprocal lattice vectors
	type(BZgrid_t), intent(inout) :: BZ
	! local variables
	logical :: sd_flg1, sd_flg12, is_gamma, self_dual
	integer :: m1, m2, m3, m1_inv, m2_inv, m3_inv
	integer :: kindex, k_cnt

	if (symtxt .eq. 'invsym') then
		BZ%use_invsym = .true.
	else if (symtxt .eq. 'nosym') then
		BZ%use_invsym = .false.
	else
		call lsquit('Invalid symtxt in init_BZgrid.',-1)
	end if

	! copy reciprocal lattice vectors, which are assumed to be
	! stored as column vectors in a 3 x 3 matrix
	BZ%reclvec(1:3,1:3) = rlvec(1:3,1:3)

	! set dimensions of sample grid
	BZ%Nk_dim1 = nk1
	BZ%Nk_dim2 = nk2
	BZ%Nk_dim3 = nk3
	BZ%Nk_nosym = nk1 * nk2 * nk3

	! set k-points
	kindex = 0
	k_cnt = 0
	do m1 = 0, nk1 - 1
		m1_inv = mod(nk1 - m1, nk1)
		sd_flg1 = (m1 .eq. m1_inv)
		if (BZ%use_invsym .and. m1 .gt. m1_inv) cycle

		do m2 = 0, nk2 - 1
			m2_inv = mod(nk2 - m2, nk2)
			sd_flg12 = sd_flg1 .and. (m2 .eq. m2_inv)
			if (BZ%use_invsym .and. sd_flg1 .and. m2 .gt. m2_inv) cycle

			do m3 = 0, nk3 - 1
				m3_inv = mod(nk3 - m3, nk3)
				if (BZ%use_invsym .and. sd_flg12 .and. m3 .gt. m3_inv) cycle

				! determine if the k-point is (i) the Gamma point,
				! (ii) self-dual so that +k is equivalent to -k
				is_gamma = (m1 .eq. 0 .and. m2 .eq. 0 .and. m3 .eq. 0)
				self_dual = (m1 .eq. m1_inv) .and. (m2 .eq. m2_inv) &
					& .and. (m3 .eq. m3_inv)

				! store k-point data
				kindex = kindex + 1
				if (kindex .gt. Max_kpoints) then
					call lsquit('Max k-points exceeded in pbc_init_BZgrid.',-1)
				end if

				BZ%kpnt(kindex)%ix_orig = -1
				BZ%kpnt(kindex)%n(1:3) = (/ m1, m2, m3 /)
				BZ%kpnt(kindex)%ninv(1:3) = (/ m1_inv, m2_inv, m3_inv /)
				BZ%kpnt(kindex)%lambda(1) = real(m1, realk) / nk1
				BZ%kpnt(kindex)%lambda(2) = real(m2, realk) / nk2
				BZ%kpnt(kindex)%lambda(3) = real(m3, realk) / nk3
				BZ%kpnt(kindex)%is_gamma = is_gamma
				BZ%kpnt(kindex)%self_dual = self_dual
				BZ%kpnt(kindex)%nsingular = 0
				if (BZ%use_invsym .and. .not. self_dual) then
					k_cnt = k_cnt + 2
					BZ%kpnt(kindex)%weight = 2.0_realk
				else
					k_cnt = k_cnt + 1
					BZ%kpnt(kindex)%weight = 1.0_realk
				end if
			end do
		end do
	end do

	BZ%Nk = kindex

	! error checking
	if (k_cnt .ne. BZ%Nk_nosym) then
		write(LUPRI,*) 'k_cnt, BZ%Nk_nosym = ',k_cnt,BZ%Nk_nosym
		call lsquit('Unexpected k-point count in init_BZgrid.',-1)
	end if

	call init_pbc_elstr(bz%fck,ndim1,ndim2)
	call init_pbc_elstr(bz%Smat,ndim1,ndim2)
	!call mem_alloc(bz%keigv,bz%nk*ndim1)

END SUBROUTINE pbc_init_BZgrid

!> \author Johannes Rekkedal
!> \date 2013
!> \brief ???
!> \param bz ???
SUBROUTINE pbc_end_BZgrid(BZ)
	implicit none
	! input and output arguments
	type(BZgrid_t), intent(inout) :: BZ
	! local variables
	integer :: i

	call free_pbc_elstr(bz%fck)
	call free_pbc_elstr(bz%Smat)
	do i=1,Bz%nk
		call mem_dealloc(bz%kpnt(i)%eigv)
	enddo
END SUBROUTINE pbc_end_BZgrid

!> \author Johannes Rekkedal
!> \date 2013 ???
!> \brief ???
!> \param indx ???
!> \param kvec ???
subroutine pbc_get_kpoint(indx,kvec)
	implicit none

	integer, intent(in) :: indx
	real(realk), intent(out) :: kvec(3)

	integer :: ix,n1,n2,n3
	!real(realk) :: lambda(3)

	!logical :: pbc_is_dim_active

	if ((indx .lt. 1) .or. (indx .gt. lat_data%num_kpoints)) then
		call lsquit('pbc_get_kpoint: Invalid k-point index.',-1)
		!else if (.not. lat_data%latvec_ok) then
		!   call lsquit('pbc_get_kpoint: No lattice vectors specified.',-1)
	else if (lat_data%num_kpoints .le. 0) then
		call lsquit('pbc_get_kpoint: No k-space grid specified.',-1)
	end if

	!  if (.not. pbc_is_dim_active(1) .and. lat_data%num_k1 .ne. 1) then
	!     write(LUPRI,*)'pbc_get_kpoint: Dim. 1 is inactive. ', &
	!          & 'Exactly 1 k-point along this dimension is required.'
	!     call quit('Invalid k point grid.')
	!  else if (.not. pbc_is_dim_active(2) .and. lat_data%num_k2 .ne. 1) then
	!     write(LUPRI,*)'pbc_get_kpoint: Dim. 2 is inactive. ', &
	!          & 'Exactly 1 k-point along this dimension is required.'
	!     call quit('Invalid k point grid.')
	!  else if (.not. pbc_is_dim_active(3) .and. lat_data%num_k3 .ne. 1) then
	!     write(LUPRI,*)'pbc_get_kpoint: Dim. 3 is inactive. ', &
	!          & 'Exactly 1 k-point along this dimension is required.'
	!     call quit('Invalid k point grid.')
	!  end if
	! convert 1-tuple index to 3-tuple index
	ix = indx - 1
	n3 = modulo(ix,lat_data%num_k3)
	ix = (ix - n3) / lat_data%num_k3
	n2 = modulo(ix,lat_data%num_k2)
	ix = (ix - n2) / lat_data%num_k2
	n1 = modulo(ix,lat_data%num_k1)
	if (ix .ne. n1) then
		write(*,*) 'pbc_get_kpoint: Error when computing k point index.'
		write(*,*) 'dims: ',lat_data%num_k1,lat_data%num_k2,lat_data%num_k3
		write(*,*) 'index, n1, n2, n3, ix',indx,n1,n2,n3,ix
		call lsquit('pbc_get_kpoint: Indexing error.',-1)
	end if

	! Translate k-points into the first Brillouin zone
	if (2*n1 .gt. lat_data%num_k1) &
		& n1 = n1 - lat_data%num_k1
	if (2*n2 .gt. lat_data%num_k2) &
		& n2 = n2 - lat_data%num_k2
	if (2*n3 .gt. lat_data%num_k3) &
		& n3 = n3 - lat_data%num_k3

	kvec(:) = n1 * lat_data%reclatvec(:,1) / lat_data%num_k1
	kvec(:) = kvec(:) + n2 * lat_data%reclatvec(:,2) / lat_data%num_k2
	kvec(:) = kvec(:) + n3 * lat_data%reclatvec(:,3) / lat_data%num_k3
end subroutine pbc_get_kpoint

!> \author Johannes Rekkedal
!> \date 2013
!> \brief ???
!> \param kdep ???
!> \param ndim1 ???
!> \param ndim2 ???
SUBROUTINE init_pbc_elstr(kdep,ndim1,ndim2)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: ndim1,ndim2
	TYPE(pbc_elstr_t),intent(INOUT) :: kdep

	allocate(kdep%kfockvec(ndim1*ndim2)) !todo use mem alloc
	allocate(kdep%kfockmat(ndim1,ndim2))
	allocate(kdep%koverlapvec(ndim1*ndim2))
	allocate(kdep%koverlapmat(ndim1,ndim2))
	allocate(kdep%kcdensityvec(ndim1*ndim2))
	allocate(kdep%kcdensitymat(ndim1,ndim2))
	allocate(kdep%kddensityvec(ndim1*ndim2))
	allocate(kdep%kddensitymat(ndim1,ndim2))
	allocate(kdep%zelms(ndim1*ndim2))
	allocate(kdep%keigv(ndim1))

	call zero_pbc_elstr(kdep)

END SUBROUTINE init_pbc_elstr

!> \author Johannes Rekkedal
!> \date 2013
!> \brief initialize the type(pbc_elstr_t)
!> \param kdep ???
SUBROUTINE zero_pbc_elstr(kdep)
	IMPLICIT NONE
	TYPE(pbc_elstr_t),intent(INOUT) :: kdep

	kdep%kfockvec=CMPLX(0.0_realk,0.0_realk,complexk)
	kdep%kfockmat=CMPLX(0.0_realk,0.0_realk,complexk)
	kdep%zelms=CMPLX(0.0_realk,0.0_realk,complexk)
	kdep%koverlapvec=CMPLX(0.0_realk,0.0_realk,complexk)
	kdep%keigv=CMPLX(0.0_realk,0.0_realk,complexk)
	kdep%koverlapmat=CMPLX(0.0_realk,0.0_realk,complexk)
	kdep%kddensityvec=CMPLX(0.0_realk,0.0_realk,complexk)
	kdep%kddensitymat=CMPLX(0.0_realk,0.0_realk,complexk)
	kdep%kcdensityvec=CMPLX(0.0_realk,0.0_realk,complexk)
	kdep%kcdensitymat=CMPLX(0.0_realk,0.0_realk,complexk)
END SUBROUTINE zero_pbc_elstr

!> \author Johannes Rekkedal
!> \date 2013
!> \brief free the type(pbc_elstr_t)
!> \param kdep ???
SUBROUTINE free_pbc_elstr(kdep)
	IMPLICIT NONE
	TYPE(pbc_elstr_t),intent(INOUT) :: kdep

	deallocate(kdep%kfockvec) !todo use mem_dealloc
	deallocate(kdep%kfockmat)
	deallocate(kdep%koverlapvec)
	deallocate(kdep%koverlapmat)
	deallocate(kdep%kcdensityvec)
	deallocate(kdep%kcdensitymat)
	deallocate(kdep%kddensityvec)
	deallocate(kdep%kddensitymat)
	deallocate(kdep%zelms)
	deallocate(kdep%keigv)
END SUBROUTINE free_pbc_elstr

!> \author Johannes Rekkedal
!> \date 2013
!> \brief  ???
!> \param mat ???
!> \param n ???
!> \param m ???
!> \param vec ???
SUBROUTINE pbc_ddevectorize_mat(MAT,n,m,VEC)
	IMPLICIT NONE
	! input variables
	INTEGER, INTENT(IN) :: n,m
	REAL(realk),INTENT(INOUT) :: MAT(n,m)
	REAL(realk),INTENT(IN) :: vec(n*m)
	! local variables
	INTEGER :: i,j,k

	j=1
	i=1
	DO k=1,n*m
		if(j .gt. m) THEN
			j=1
			i=i+1
		ENDIF
		MAT(i,j)=0.0_realk
		MAT(i,j)=VEC(k)
		j=j+1
	ENDDO

END SUBROUTINE pbc_ddevectorize_mat

!> \author Johannes Rekkedal
!> \date 2013
!> \brief  ???
!> \param mat ???
!> \param n ???
!> \param m ???
!> \param vec ???
SUBROUTINE pbc_zdevectorize_mat(MAT,n,m,VEC)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n,m
	COMPLEX(complexk),INTENT(INOUT) :: MAT(n,m)
	COMPLEX(complexk),INTENT(IN) :: vec(n*m)
	!LOCAL VARIABLES
	INTEGER :: i,j,k

	i=1
	j=1
	DO k=1,n*m
		if(j .gt. m) THEN
			j=1
			i=i+1
		ENDIF
		MAT(i,j)=CMPLX(0.0_realk,0.0_realk,complexk)
		MAT(i,j)=VEC(k)
		j=j+1
	ENDDO

END SUBROUTINE pbc_zdevectorize_mat

!> \author Johannes Rekkedal
!> \date 2013
!> \brief  ???
!> \param filename ???
!> \param matop ???
!> \param siz ???
!> \param reset_flag ???
!> \param transpose_flag ???
subroutine pbc_readopmat2(filename,matop,siz,reset_flag,transpose_flag)
	implicit none
	! Input and output variables
	logical, intent(in) :: reset_flag, transpose_flag
	integer, intent(in) ::  siz
	real(realk), intent(inout) :: matop(siz,siz)
	character*15,intent(in) :: filename

	integer, parameter :: oldintk = SELECTED_INT_KIND(9)

	! local variables
	integer, parameter :: buflen = 600
	integer :: length, I
	integer :: LUPBC_loc
	!logical :: existflag

	integer :: ibuf(buflen)
	real(realk) :: buf(buflen)

	integer :: row,col        !,mat_index

	! safety check
	if (siz .lt. 1 .or. siz .gt. 10000) then
		call lsquit('Unreasonable siz in pbc_readopmat2.',-1)
	end if

	! get file name
	!call pbcmatfilnam(filename,latpos1,latpos2,latpos3,optype)

	!  write(LUPRI,'(A,A,A,I3,I3,I3,A,A,A)') 'pbc_readopmat: Reading ',optype,'(l = [',latpos1,latpos2,latpos3,' ] ) from file ',filename

	!call lsINQ(filename,'EXIST',existflag)
	!if (existflag) then
	LUPBC_loc = -1
	call lsopen(LUPBC_loc,filename,'OLD','UNFORMATTED')
	!else
	!   write(*,'(/3A)') 'pbc_readopmat: Error ',filename,' does not exist'
	!   call lsquit('pbc_readopmat: Error, non-existent file.',-1)
	!end if

	! Read overlap
	rewind LUPBC_loc

	if (reset_flag) then
		matop(:,:) = 0.0_realk
	end if
	length = 1
	do while (length .ge. 0)
		read(LUPBC_loc) (buf(I),I=1,buflen),(ibuf(I),I=1,buflen),length
		do I = 1, length
#if 1
			if (transpose_flag) then
				col = 1 + mod(ibuf(I)-1,siz)
				row = 1 + (ibuf(I) - col) / siz
				!mat_index = (col-1)*siz + row
			else
				!mat_index = ibuf(I)
				col = 1 + mod(ibuf(I)-1,siz)
				row = 1 + (ibuf(I) - col) / siz
			end if

			if (reset_flag) then
				matop(row,col) = buf(I)
			else
				matop(row,col) = matop(row,col) + buf(I)
			end if
#else
			! older code
			if (reset_flag) then
				matop(ibuf(I),1) = buf(I)
			else
				matop(ibuf(I),1) = matop(ibuf(I),1) + buf(I)
			end if
#endif
		end do
	end do
	call lsclose(LUPBC_loc,'KEEP')

end subroutine pbc_readopmat2

!> \author Johannes Rekkedal
!> \date 2013
!> \brief  ???
!> \param fockmat ???
!> \param siz ???
!> \param lvec ???
!> \param timestep ???
!> \param reset_flag ???
!> \param sfac ???
!> \param lupri ???
subroutine pbc_read_diis_fock(fockmat,siz,lvec,timestep,reset_flag,sfac,lupri)
	implicit none
	integer, parameter :: intk = selected_int_kind(9)
	logical, intent(in) :: reset_flag
	integer,intent(in)  :: lupri
	integer, intent(in) :: timestep, siz, lvec(3)
	real(realk), intent(inout) :: fockmat(siz,siz), sfac

	logical :: existflag
	character(LEN=50) :: filename
	integer :: fnamlen, timestep_file, siz_file
	integer :: dummy_int, mu, nu
	integer :: ludiis
	real(realk) :: fock_elem

	! error check
	if (siz .lt. 1 .or. siz .gt. 5000) then
		call lsquit('Bad Fock matrix size in DIIS.',lupri)
	end if

	existflag=.true.
	! assign filename and open file
	call pbc_diis_filename(filename,fnamlen,timestep,lvec,'FOCK',lupri)
	!call GPINQ(filename(1:fnamlen),'EXIST',existflag)
	if (existflag) then
		LUDIIS = -1
		call lsopen(LUDIIS,filename(1:fnamlen),'OLD','UNFORMATTED')
	else
		write(LUPRI,'(/3A)') 'DIIS read error: ',filename,' does not exist'
		call lsquit('DIIS read error.',lupri)
	end if

	! rewind and read header info
	rewind LUDIIS
	read(LUDIIS) timestep_file, siz_file

	if (timestep_file .ne. timestep) then
		call lsquit('Inconsistent time index on DIIS Fock matrix file.',lupri)
	else if (siz_file .ne. siz) then
		call lsquit('Inconsistent size info on DIIS Fock matrix file.',lupri)
	end if
	! read Fock matrix
	loop_ao1: do mu = 1, siz
		loop_ao2: do nu = 1, siz
			read(LUDIIS) fock_elem
			if (reset_flag) then
				fockmat(mu,nu) = fock_elem * sfac
			else
				fockmat(mu,nu) = fockmat(mu,nu) + fock_elem * sfac
			end if
		end do loop_ao2
	end do loop_ao1

	! close file
	call lsclose(LUDIIS,'KEEP')

end subroutine pbc_read_diis_fock

!> \author Johannes Rekkedal
!> \date 2013
!> \brief  ???
!> \param fname ???
!> \param fnamlen ??? 
!> \param timestep ???
!> \param lvec ???
!> \param ctype ???
!> \param lupri ???
subroutine pbc_diis_filename(fname,fnamlen,timestep,lvec,ctype,lupri)
	implicit none
	character(LEN=21), intent(inout) :: fname
	integer,intent(in) :: lupri
	integer, intent(out) :: fnamlen
	integer, intent(in) :: timestep, lvec(1:3)
	character*(*), intent(in) :: ctype

	integer :: k

	if (ctype .eq. 'FOCK') then
		! create filename of the form 'diisfckIT_LX_LY_LZ'
		! example: IT = 5, LX = 2, LY = 0, LZ = -1
		!          gives 'diisfck05_p02_p00_n01'
		fname(1:7) = 'diisfck'
		call pbc_num2str(fname(8:9),timestep,2,lupri)
		fnamlen = 10

		do k = 1, 3
			if (lvec(k) .lt. 0) then
				fname(fnamlen:(fnamlen+1)) = '_n'
			else
				fname(fnamlen:(fnamlen+1)) = '_p'
			end if
			call pbc_num2str(fname((fnamlen+2):(fnamlen+3)),abs(lvec(k)),2,lupri)
			fnamlen = fnamlen + 4
		end do
		fnamlen = 21
	else if (ctype .eq. 'ERRVEC') then
		! create filename of the form 'diis_errvecIT'
		fname(1:11) = 'diis_errvec'
		call pbc_num2str(fname(12:13),timestep,2,lupri)
		fnamlen = 13
	else
		call lsquit('Invalid type supplied to pbc_diis_filename.',lupri)
	end if
end subroutine pbc_diis_filename

!> \author Johannes Rekkedal
!> \date 2013
!> \brief This subroutine converts a number to a zero-padded string
!> \brief of digits. E.g., "call pbc_num2str(txt,117,2)" writes '17' to txt.
!> \param str ???
!> \param num ???
!> \param numdigits ???
!> \param lupri ???
subroutine pbc_num2str(str,num,numdigits,lupri)
	implicit none

	integer,intent(in)  :: lupri
	integer, intent(in) :: num
	integer :: numdigits
	character(LEN=numdigits), intent(inout) :: str

	character(LEN=10), parameter :: digit_symbol = '0123456789'
	integer :: i, digit, new_num

	! Some safety checks
	if (numdigits .le. 1 .or. numdigits .gt. 4) then
		call lsquit('Unreasonably many digits in num2str conversion.',lupri)
	else if (num .lt. 0) then
		call lsquit('Negative numbers not supported in num2str conversion.',lupri)
	end if

	! loop from least to most significant digit
	new_num = num
	do i = numdigits, 1, -1
		digit = 1+mod(new_num,10)
		str(i:i) = digit_symbol(digit:digit)
		new_num = new_num / 10
	end do
end subroutine pbc_num2str

!> \author Johannes Rekkedal
!> \date 2013
!> \brief  ???
!> \param bz ???
!> \param nrealvec ???
!> \param ll ???
!> \param kdep ???
!> \param ndim ???
SUBROUTINE convert_2_kmat(bz,nrealvec,ll,kdep,ndim)
	INTEGER,intent(in):: ndim,nrealvec
	TYPE(BZgrid_t),intent(In) :: BZ
	TYPE(lvec_list_t),intent(in) :: ll
	TYPE(pbc_elstr_t),intent(inout) :: kdep(1)!(bz%nk)
	!Local variables
	!real(realk) :: pi=3.14159265
	real(realk) :: sumrealfock(ndim*ndim)
	INTEGER :: nlat,kpt,i,j,l1,l2,l3
	real(realk) :: phase1,phase2,phase3
	COMPLEX(complexk) :: phasetot
	sumrealfock=0.0
	DO kpt=1,1!bz%nk
		DO nlat=1,nrealvec
			call find_latt_vectors(nlat,l1,l2,l3,ll)
			if(abs(l1) .gt. ll%nneighbour) CYCLE
			if(abs(l2) .gt. ll%nneighbour) CYCLE
			if(abs(l3) .gt. ll%nneighbour) CYCLE
			phase1=bz%kpnt(kpt)%lambda(1)*ll%lvec(nlat)%lat_coord(1)
			phase2=bz%kpnt(kpt)%lambda(2)*ll%lvec(nlat)%lat_coord(2)
			phase3=bz%kpnt(kpt)%lambda(3)*ll%lvec(nlat)%lat_coord(3)
			phasetot=CMPLX(0.,(phase1+phase2+phase3)*2.*pi,complexk)
			!  do i=1,ndim*ndim
			!   kdep(kpt)%kfockvec(i)=kdep(kpt)%kfockvec(i)+&
			!   ll%lvec(nlat)%fck_vec(i)*exp(phasetot)
			!  enddo
		enddo
	enddo


	DO nlat=1,nrealvec
		call find_latt_vectors(nlat,l1,l2,l3,ll)
		if(abs(l1) .gt. ll%nneighbour) CYCLE
		if(abs(l2) .gt. ll%nneighbour) CYCLE
		if(abs(l3) .gt. ll%nneighbour) CYCLE
		!do i=1,ndim*ndim
		! sumrealfock(i)=sumrealfock(i)+ll%lvec(nlat)%fck_vec(i)
		!enddo
	enddo

	!write(*,*) 'sum over fock matrices'
	!DO j=1,ndim
	!       write(*,*) (sumrealfock(i+(j-1)*ndim),i=1,ndim)
	!ENDDO
	stop
END SUBROUTINE convert_2_kmat

!> \author Johannes Rekkedal
!> \date 2013
!> \brief ???
!> \param bz ???
!> \param nbast ???
!> \param mattxt ???
SUBROUTINE print_bands(bz,nbast,mattxt)
	IMPLICIT NONE
	! input and output
	INTEGER, INTENT(IN) :: nbast
	CHARACTER(len=11),INTENT(IN) :: mattxt
	TYPE(BZgrid_t),intent(iN) :: bz
	! local variables
	REAL(realk) :: kpt,kvec(3)
	INTEGER :: i,k,iunit,nonsingdim
	CHARACTER(LEN=3) :: nline


	iunit=-1
	CALL lsOPEN(IUNIT,mattxt,'unknown','FORMATTED')
	nline='no'
	DO k=1,(bz%nk+1)/2
                nonsingdim=nbast-bz%kpnt(k)%nsingular
		nline='no'
		call pbc_get_kpoint(k,kvec)
                kpt=kvec(1)**2+kvec(2)**2+kvec(3)**2
                kpt=sqrt(kpt)
		write(iunit,101,advance=nline) kpt !must convert to kpoint angle
		101 FORMAT(E12.4)
		DO i=1,nbast
			if(i .eq. nonsingdim) nline = 'yes'
			!write(iunit,100,advance=nline) real(kdep(k)%keigv(i))
			write(iunit,100,advance=nline) bz%kpnt(k)%eigv(i)
			100 FORMAT(E18.8)
		ENDDO
	ENDDO
	call lsclose(iunit,'KEEP')

END SUBROUTINE print_bands

END MODULE PBC_MSC

