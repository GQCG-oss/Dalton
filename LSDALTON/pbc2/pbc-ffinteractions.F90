#ifdef MOD_UNRELEASED

MODULE pbc_ff_contrib
use files
USE precision
USE TYPEDEF
USE lattice_type
USE lattice_vectors
USE matrix_module
USE matrix_operations
USE pbc_MSC
USE pbc_matrix_operations
USE ls_util
USE integralinterfaceMod
USE lstiming
CONTAINS

SUBROUTINE pbc_controlmm(num_its,Tlat,Tlmax,lmax,square_intermed,latvec, &
		& nbast,lupri,nfdensity,numvecs,ll,E_ff,E_nn,refcell)
IMPLICIT NONE

  TYPE(moleculeinfo) :: refcell
  INTEGER,INTENT(IN) :: lupri
  INTEGER, INTENT(IN) :: num_its,Tlmax
  INTEGER, INTENT(IN) :: lmax!,nfsze
  INTEGER,INTENT(IN) :: nbast,numvecs
  TYPE(lvec_list_t), INTENT(INOUT) :: ll
  REAL(realk),INTENT(INOUT) :: E_ff,E_nn
  REAL(realk), INTENT(INOUT) :: Tlat((1+lmax)**2,(1+lmax)**2)
  TYPE(matrix) :: nfdensity(numvecs)
  REAL(realk),INTENT(IN) ::latvec(3,3)
  LOGICAL, INTENT(IN) :: square_intermed
  REAL(realk),POINTER :: nucmom(:)
  REAL(realk) :: fock_tmp(nbast**2)
  CHARACTER(LEN=18) :: mattxt
  INTEGER :: i,j
  INTEGER :: num_latvec
  INTEGER :: iunit
!cell_mm is multipole moment as input
!it giveback the local moment of the central cell in a structure
!I am not sure whether I need it yet.
!call pbc_calc_cell_locmom(cell_mm,lmax) 


!Does the iterations in Eq. 24 J.Chem. Phys., Vol. 121 No. 7 (2004)
!This function calls pbc_WNF_matrix, pbc_TsupNF_matrix
!and pbc_iterate_T_lattice
call pbc_form_Tlattice(num_its,Tlat,tlmax,square_intermed,latvec,ll%nf,lupri)

!Does the iterations in Eq. 24 J.Chem. Phys., Vol. 121 No. 7 (2004)
!call pbc_get_nfsize(n1,n2,n3,ll%nneighbour,lupri)
!nfsze=(2*n1+1)*(2*n2+1)*(2*n3+1)
!write(*,*) 'nfsze',nfsze
!num_latvec = size(ll%lvec)
!call mem_alloc(nucmom,(lmax+1)**2)
!call pbc_comp_nucmom(refcell,nucmom,lmax,nfsze,lupri)

!fock_tmp=0.0_realk
!call pbc_fform_fck(Tlmax,tlat,lmax,nbast,nfsze,ll,nfdensity,nucmom,g_2,E_ff,&
!                   E_nn,lupri)

!call mem_dealloc(nucmom)

!WRITE TLAT OUT TO DISK AND THEN READ IT FROM FILE IN THE scf iterations
!  write(*,*) nbast,ndim
  iunit=1
  write(mattxt,'(A18)') 'Tlatticetensor.dat'
  CALL lsOPEN(IUNIT,mattxt,'UNKNOWN','UNFORMATTED')
  write(iunit) (1+lmax)**2!,nil
  DO j=1,(1+lmax)**2
    write(iunit) (Tlat(j,i),i=1,(1+lmax)**2)
  ENDDO
  CALL lsCLOSE(IUNIT,'KEEP')






END SUBROUTINE pbc_controlmm

SUBROUTINE &
&pbc_form_Tlattice(num_its,Tlat,lmax,square_intermed,latvec,nfs,lupri)
  implicit none

  integer,intent(in) :: lupri
  integer, intent(in) :: num_its,nfs
  integer, intent(in) ::  lmax
  real(realk), intent(inout) :: Tlat((1+lmax)**2,(1+lmax)**2)
  real(realk),intent(in) :: latvec(3,3)
  logical, intent(in) :: square_intermed

  integer :: siz, i, row, col, l, qmin, qmax
  integer :: nfsiz, nfsiz1, nfsiz2, nfsiz3
  real(realk) :: layers_in_T
  real(realk), allocatable :: W_NF(:,:), T_supNF(:,:), Tlat_aux(:,:)

  integer, parameter ::  debug_siz = 36
  integer, parameter :: debug_lmax = 5
  integer :: debug_nf
  real(realk) :: Tdebug(debug_siz,debug_siz), Tdebug2(debug_siz,debug_siz)
#ifdef DEBUGPBC
  real(realk) :: normsq
  INTEGER :: j
#endif
!  logical, external :: pbc_3dims_active

  siz = (1+lmax)**2
 
#ifdef DEBUGPBC
  write(*,*) 'In subroutine pbc_form_tlattice'
  write(*,*) 'lmax and size',lmax,siz
#endif

  if (num_its .lt. 0 .or. num_its .gt. 25) then
     call lsquit('pbc_form_Wlattice: Unreasonable number of iterations requested.',lupri)
  end if

  allocate(W_NF(siz,siz))
  allocate(T_supNF(siz,siz))
  allocate(Tlat_aux(siz,siz))
  if (.not. allocated(W_NF)) then
     call lsquit('pbc_form_Tlattice: Failed to allocate memory for a W tensor.'&
     &,lupri)
  else if (.not. (allocated(T_supNF) .and. allocated(Tlat_aux))) then
     call lsquit('pbc_form_Tlattice: Failed to allocate memory for a T tensor.'&
     &,lupri)
  end if
  ! For debugging purposes we can call pbc_TW_test here
  ! call pbc_TW_test(lmax,W_NF,Tlat_aux,Tlat,square_intermed)

  ! get the near-field size
  call pbc_get_nfsize(nfsiz1,nfsiz2,nfsiz3,nfs,lupri)
  if (nfsiz1*nfsiz2 .ne. 0 .and. nfsiz1 .ne. nfsiz2) then
     call lsquit('pbc_form_Tlattice: All NF sizes must be equal.',lupri)
  else if (nfsiz1*nfsiz3 .ne. 0 .and. nfsiz1 .ne. nfsiz3) then
     call lsquit('pbc_form_Tlattice: All NF sizes must be equal.',lupri)
  else if (nfsiz2*nfsiz3 .ne. 0 .and. nfsiz2 .ne. nfsiz3) then
     call lsquit('pbc_form_Tlattice: All NF sizes must be equal.',lupri)
  end if

  nfsiz = max(max(nfsiz1,nfsiz2),nfsiz3)

  write(LUPRI,*) 'Generating the crystal far-field interaction matrix Tlattice'
  write(LUPRI,'(A,I4)') 'for the near-field size ',nfsiz

  ! get tensors representing the geometrical structure of
  ! the near-field and the near-field of the first supercell.
  call pbc_WNF_matrix(W_NF,lmax,siz,latvec,lupri)
  !!!call pbc_restore_W_matrix(lmax,W_NF)
#ifdef DEBUGPBC
  normsq = 0
  do i =1,siz
   do j=1,siz
    normsq=normsq+W_nf(i,j)**2
   enddo
  enddo
  write(*,*) 'Norm for WNF matrix = ', normsq
  write(lupri,*) 'Norm for WNF matrix = ', normsq
#endif
     
    

  call pbc_TsupNF_matrix(T_supNF,nfsiz,lmax,siz,square_intermed,latvec,lupri)

#ifdef DEBUGPBC
  normsq = 0
  do i =1,siz
   do j=1,siz
    normsq=normsq+T_supNF(i,j)**2
   enddo
  enddo
  write(*,*) 'Norm for T_supNF matrix = ', normsq
  write(lupri,*) 'Norm for WNF matrix = ', normsq
#endif

  !call pbc_Tmatrix_print(T_supNF,siz,'T_supNF',lupri)
!  debugsumT=0.0_realk
!  do i=1,siz
!   do l=1,siz
!      debugsumT=debugsumT+abs(T_supNF(i,l))
!   enddo
!  enddo
!
!  if(debugsumT .gt. 0.0_realk) then
!    write(*,*) 'debugsumT not equal to 0.0_realk',siz,lmax
!    write(*,*)  debugsumT
!    !call write_matrix(T_supNF,siz,siz)
!    stop
!  endif
  ! iterate Kudin & Scuseria's equation

  Tlat_aux(:,:) = T_supNF(:,:)
  layers_in_T = 3.0_realk * nfsiz + 1.0_realk
  do i = 1, num_its
     write(LUPRI,*) 'pbc_form_Tlattice(): iteration ',i
     if (mod(i,2) .eq. 1) then
        call pbc_iterate_Tlattice(Tlat,Tlat_aux,W_NF,T_supNF,lmax,siz,square_intermed,lupri)
     else
        call pbc_iterate_Tlattice(Tlat_aux,Tlat,W_NF,T_supNF,lmax,siz,square_intermed,lupri)
     end if
     ! count the number of cell layers included in the T tensor
     layers_in_T = 3.0_realk * layers_in_T + 1.0_realk
     ! some debug code
     if (.false.) then
        if (i .eq. 1) then
           debug_nf = 3*nfsiz + 1
           call pbc_TsupNF_matrix(Tdebug,debug_nf,debug_lmax,debug_siz,&
		   & square_intermed,latvec,lupri)
           Tdebug(:,:) = Tdebug(:,:) + T_supNF(1:debug_siz,1:debug_siz)

!           call pbc_matdiff_db(Tdebug,debug_siz,debug_siz,Tlat,siz,siz,15,15,'T1_iter','T1_direct',.not. square_intermed)
        else if (i .eq. 2) then
           debug_nf = 3*(3*nfsiz + 1) + 1
           call pbc_TsupNF_matrix(Tdebug2,debug_nf,debug_lmax,debug_siz,&
		   & square_intermed,latvec,lupri)
           Tdebug2(:,:) = Tdebug2(:,:) + Tdebug(:,:)

!           call pbc_matdiff_db(Tdebug2,debug_siz,debug_siz,Tlat_aux,siz,siz,15,15,'T2_iter','T2_direct',.not. square_intermed)
        end if
     end if
     ! end debug code
  end do
  if (mod(num_its,2) .eq. 0) then
     Tlat = Tlat_aux
  end if

!!! I have to put this on
!  if (reset_Tdiv) then
!     if (pbc_3dims_active()) then
!        Tlat(1:4,1:4) = 0.0_realk
!        write(LUPRI,*) 'Resetting charge-charge and charge-dipole interactions in Tlattice.'
!     else
        Tlat(1,1) = 0.0_realk
!        write(LUPRI,*) 'Resetting charge-charge.'
!     end if
!  end if


  write(LUPRI,'(A,I3,A)') 'pbc_form_Tlattice: Tlattice was iterated ',num_its,' times.'
  write(LUPRI,'(A,E10.3,A)') 'pbc_form_Tlattice: Tlattice includes interactions for ',layers_in_T,' cell layers.'

  write(*,'(A,I3,A)') 'pbc_form_Tlattice: Tlattice was iterated ',num_its,' times.'
  write(*,'(A,E10.3,A)') 'pbc_form_Tlattice: Tlattice includes interactions for ',layers_in_T,' cell layers.'

  ! put in data in upper triangle
  if (.not. square_intermed) then
     write(LUPRI,*) 'pbc_form_Tlattice: Restoring some elements in upper triangle of Tlattice.'
     do l = 0, lmax
        qmin = l*(l+1) - l + 1
        qmax = l*(l+1) + l + 1
        do row = qmin, qmax
           do col = row + 1, qmax
              Tlat(row, col) = Tlat(col, row)
           end do
        end do
     end do
  end if

#ifdef DEBUGPBC
  normsq = 0
  do i =1,siz
   do j=1,siz
    normsq=normsq+Tlat(i,j)**2
   enddo
  enddo
  write(*,*) 'Norm for Tlat matrix = ', normsq
  write(lupri,*) 'Norm for Tlat matrix = ', normsq
#endif

  call pbc_T_enforce_invsym(Tlat,lmax,siz,lupri)

#ifdef DEBUGPBC
  normsq = 0
  do i =1,siz
   do j=1,siz
    normsq=normsq+Tlat(i,j)**2
   enddo
  enddo
  write(*,*) 'Norm for Tlat matrix invsym = ', normsq
  write(lupri,*) 'Norm for Tlat matrix invsym = ', normsq
#endif

!  call pbc_matlab_print(Tlat,256,256,'Tlattice before trans',lupri)

  call pbc_restore_T_matrix(lmax,Tlat)

#ifdef DEBUGPBC
  normsq = 0
  do i =1,siz
   do j=1,siz
    normsq=normsq+Tlat(i,j)**2
   enddo
  enddo
  write(*,*) 'Norm for Tlat matrix restore = ', normsq
  write(lupri,*) 'Norm for Tlat matrix invsym = ', normsq
 ! call pbc_Tmatrix_print(Tlat,siz,'Tlattice',lupri)
#endif
  
!  call pbc_matlab_print(Tlat,256,256,'Tlattice after',lupri)

  !call write_matrix(Tlat,siz,siz)
  !write(*,*) 'debug',siz
  !stop
!  debugsumT=0.0_realk
!  do i=1,siz
!   do l=1,siz
!      debugsumT=debugsumT+abs(Tlat(i,l))
!   enddo
!  enddo
!
!  if(debugsumT .gt. 0.0_realk) then
!    write(*,*) 'debugsumT not equal to 0.0_realk after invsym',siz,lmax
!    write(*,*)  debugsumT
!    !call write_matrix(T_supNF,siz,siz)
!!    stop
!  endif

!  call pbc_matrix_info(siz,siz,Tlat,Tlat,.false.,'Tlattice','dummy')
  
  !call pbc_Tmatrix_print(Tlat,siz,'Tlattice',lupri)

  deallocate(W_NF)
  deallocate(T_supNF)
  deallocate(Tlat_aux)

END SUBROUTINE pbc_form_Tlattice


! This subroutine calculates the M2M translation operator T_NF
! that translates the moments of the central cell to cells in
! the near-field. The central cell is included.
!
! Note that the argument siz is there for convenience and must
! always have the value (1+lmax)**2
!
subroutine pbc_WNF_matrix(W_NF,lmax,siz,latvec,lupri)
  implicit none

  integer, intent(in) :: siz
  integer, intent(in) :: lmax
  integer, intent(in) :: lupri
  real(realk), intent(inout) :: W_NF((1+lmax)**2,(1+lmax)**2)
  real(realk),intent(in) :: latvec(3,3)

  integer :: mx, my, mz, p, m_max1, m_max2, m_max3
  real(realk) :: pos_std(3), pos_lat(3)
  real(realk), allocatable :: W_aux(:,:)

  !logical, external :: ifpbc_active

  if (siz .ne. (1+lmax)**2) call lsquit('pbc_WNF_matrix: Inconsistent sizes.'&
  &,lupri)

  allocate(W_aux(siz,siz))  
  if (.not. allocated(W_aux)) then
     call lsquit('pbc_WNF_matrix: Failed to allocate memory for a W tensor.'&
     &,lupri)
  end if

  ! initialize
  W_NF(:,:) = 0.0_realk

  ! if the periodicity is turned off in some directions we need
  ! to constrain the loops accordingly
  m_max1 = 0
  if (ifpbc_active(1,lupri)) m_max1 = 1
  m_max2 = 0
  if (ifpbc_active(2,lupri)) m_max2 = 1
  m_max3 = 0
  if (ifpbc_active(3,lupri)) m_max3 = 1

  ! loop over first super cell (it is the near-field when nfsiz = 1)
  do mx = -m_max1, m_max1
  do my = -m_max2, m_max2
  do mz = -m_max3, m_max3
     if (abs(mx) + abs(my) + abs(mz) .ne. 0) then
        pos_lat(:) = -1.0_realk * (/ mx, my, mz /)
        !call pbc_lat2std_coord(pos_lat, pos_std)
        call latt_2_std_coord(pos_lat, pos_std,latvec)
        W_aux(:,:) = 0.0_realk
        call pbc_get_ltsqr_W_matrix(lmax,pos_std,W_aux)
        W_NF = W_NF + W_aux
     else
        ! The FMM code routines that we call are picky and will
        ! not accept the task of constructing a translation matrix
        ! for a zero vector. Fortunately, this is trivial so we add
        ! the contribution directly here.
        do p = 1,siz
           W_NF(p,p) = W_NF(p,p) + 1.0_realk
        end do
     end if
  end do
  end do
  end do

  deallocate(W_aux)

end subroutine pbc_WNF_matrix

! Regard the NF of the central cell as a supercell. This subroutine
! calculates the M2L translation operator T_supNF that
! transforms the moments of the central cell into local moments of
! the near-field of this supercell.
subroutine pbc_TsupNF_matrix(T_supNF,nfsiz,lmax,siz,square_flag,latvec,lupri)
  implicit none
  
  integer, intent(in) :: nfsiz
  integer, intent(in) :: lmax, siz
  integer, intent(in) :: lupri
  real(realk), intent(inout) :: T_supNF(siz,siz)
  real(realk),intent(in) ::latvec(3,3)
  logical, intent(in) :: square_flag

  integer :: m_min, m_max, mx, my, mz, layer
  integer :: fac1, fac2, fac3
  real(realk) :: pos_std(3), pos_lat(3)
  real(realk), allocatable :: T_aux(:,:)


  !logical, external :: ifpbc_active

  if (siz .ne. (1+lmax)**2) call lsquit('pbc_TsupNF_matrix: Inconsistent sizes.',lupri)

  if (nfsiz .lt. 1 .or. nfsiz .gt. 22) then
     call lsquit('pbc_TsupNF_matrix: Huge NF size. Bailing out.',lupri)
  end if

  allocate(T_aux(siz,siz))  
  if (.not. allocated(T_aux)) then
     call lsquit('pbc_TsupNF_matrix: Failed to allocate memory for a T tensor.'&
     &,lupri)
  end if

  ! initialize
  T_supNF(:,:) = 0.0_realk

  ! if the periodicity is turned off in some directions we need
  ! to constrain the loops accordingly
  fac1 = 0
  if (ifpbc_active(1,lupri)) fac1 = 1
  fac2 = 0
  if (ifpbc_active(2,lupri)) fac2 = 1
  fac3 = 0
  if (ifpbc_active(3,lupri)) fac3 = 1

  ! loop over near-field of super-cell, but accept
  ! only points outside the near-field of the central cell
  m_min = nfsiz + 1
  m_max = 3 * nfsiz + 1
  do layer = m_max, m_min, -1
     write(*,*) 'layer: ', layer
  do mx = -layer*fac1, layer*fac1
  do my = -layer*fac2, layer*fac2
  do mz = -layer*fac3, layer*fac3
     if ((abs(mx) .eq. layer) .or. (abs(my) .eq. layer) .or. (abs(mz) .eq. layer)) then
        pos_lat(1) = 1.0_realk * real(mx)
        pos_lat(2) = 1.0_realk * real(my)
        pos_lat(3) = 1.0_realk * real(mz)
        !call pbc_lat2std_coord(pos_lat, pos_std)
        call latt_2_std_coord(pos_lat, pos_std,latvec)
        T_aux(:,:) = 0.0_realk
        call pbc_get_FLTSQ_T_matrix(lmax,pos_std,T_aux)
        T_supNF = T_supNF + T_aux
     end if
  end do
  end do
  end do
  end do

  if (square_flag) then
     ! The FMM code gives a matrix T_{lm,jk} which is lower-triangular
     ! or, more precisely, T_{lm,jk} = 0 for j > l. Here we restore the
     ! full square matrix.
     call pbc_restore_T_matrix(lmax,T_supNF)
  end if

  deallocate(T_aux)

end subroutine pbc_TsupNF_matrix

! This subroutine does one iteration according to Eq. (24) in
! Kudin & Scuseria, J Chem Phys 121(7):2886-2890.
! That is, it computes
!
!    Wlat = m2l(scale(Wlat_old), T_NF) + W_supNF
!
subroutine pbc_iterate_Tlattice(Tlat,Tlat_old,W_NF,T_supNF,lmax,siz,square_flag&
&,lupri)
  implicit none

  real(realk), parameter :: stmp = 10.0_realk

  integer, intent(in) :: lmax, siz,lupri
  real(realk), intent(in) :: W_NF(siz,siz), T_supNF(siz,siz)
  real(realk), intent(inout) :: Tlat_old(siz,siz)
  real(realk), intent(out) :: Tlat(siz,siz)
  logical, intent(in) :: square_flag

  real(realk) :: scale_factor

  if (siz .ne. (1+lmax)**2) call lsquit('pbc_iterate_Tlattice: Inconsistent sizes.',lupri)

  !call pbc_matrix_info(siz,siz,Tlat_old,Tlat_old,.false.,'Tlat_old','dummy')

  scale_factor = 3.0_realk

  call pbc_scale_T(Tlat_old,lmax,siz,scale_factor,lupri)

  call pbc_do_m2l(Tlat,Tlat_old,W_NF,lmax,siz,square_flag,lupri)
  Tlat = Tlat + T_supNF

  !if (ffdata%reset_T_diverg) then
  !Tlat(1,1) = 0.0_realk
  !Tlat(2:4,2:4) = 0.0_realk
  !end if

end subroutine pbc_iterate_Tlattice


! This subroutine scales an interaction matrix so that the translations
! it represents become larger by a factor s.
subroutine pbc_scale_T(Tmatrix,lmax,siz,s,lupri)
  implicit none

  integer, intent(in) :: lmax, siz,lupri
  real(realk), intent(in) :: s
  real(realk), intent(inout) :: Tmatrix(siz,siz)

  integer :: l, mu, j, kappa, row, col, npow
  
  if (siz .ne. (1+lmax)**2) call lsquit('pbc_scale_T: Inconsistent sizes.',lupri)
  
  row = 0
  do l = 0, lmax
     do mu = -l, l
        ! update row index
        row = row + 1
        ! reset column index
        col = 0
        do j = 0, lmax
           do kappa = -j, j
              ! update column index
              col = col + 1
              ! which power of the scaling factor?
              npow = l + j + 1
              ! scale!
              Tmatrix(row,col) = Tmatrix(row,col) / s**npow
           end do
        end do
     end do
  end do

end subroutine pbc_scale_T

! This subroutine performs the M2L operation of T1 with W2 and
! stores the result in Tres.
subroutine pbc_do_m2l(Tres,T1,W2,lmax,siz,square_flag,lupri)
	implicit none

	logical :: square_flag
	integer, intent(in) :: lmax, siz,lupri
	real(realk), intent(in) :: T1(siz,siz), W2(siz,siz)
	real(realk), intent(out) :: Tres(siz,siz)

	real(realk) :: fac
	integer :: col, row
	integer :: l,m,j,k,p,q,kmax,qmin,lm,jk,pq

	if (siz .ne. (1+lmax)**2) call lsquit('pbc_do_m2l: Inconsistent sizes.',lupri)

	Tres(:,:) = 0.0_realk

	! compute M2L result
	if (square_flag) then
		do col = 1, siz
			do row = 1, siz
				Tres(row,col) = dot_product(T1(row,col:siz),W2(col:siz,col))
			end do
		end do
		return
	end if


	rows_l: do l = 0,lmax
		rows_m: do m = -l,l
			lm = l*(l+1)+1+m
			cols_j: do j = 0,l
				if (j .eq. l) then
					kmax = m
				else
					kmax = j
				end if
				cols_k: do k = -j, kmax
					jk = j*(j+1)+1+k
					! loop over pq >= jk
					sum_p: do p = j, lmax
						if (p .eq. j) then
							qmin = k
						else
							qmin = -p
						end if
						sum_q: do q = qmin,p
							pq = p*(p+1)+1+q
							! determine factor
							fac = 1.0_realk
							if (q .eq. 0) fac = fac / 2.0_realk
							!if (mod(p,2_long) .eq. 1) fac = -fac
							! add contribution
							Tres(lm,jk) = Tres(lm,jk) + T1(max(lm,pq),min(lm,pq)) &
								& * W2(pq,jk) * fac
						end do sum_q
					end do sum_p
					! put in the factor
					fac = 1.0_realk
					if (k .eq. 0) fac = fac / 2.0_realk
					!if (mod(j,2_long) .eq. 1) fac = -fac              
					Tres(lm,jk) = Tres(lm,jk) / fac
				end do cols_k
			end do cols_j
		end do rows_m
	end do rows_l

end subroutine pbc_do_m2l

! Even if the content of a cell lacks any symmetry the lattice
! itself will still have inversion symmetry. Therefore Tlattice
! must also have inversion symmetry. This subroutine resets elements
! in a T tensor that should be zero due inversion symmetry.
subroutine pbc_T_enforce_invsym(T,lmax,siz,lupri)

  integer, intent(in) :: lmax, siz
  real(realk), intent(inout) :: T(siz,siz)
  integer,intent(in) :: lupri

  integer :: l,j,pmin,pmax,qmin,qmax

  if (siz .ne. (1+lmax)**2) call lsquit('pbc_T_enforce_invsym: Inconsistent sizes.',lupri)

  write(LUPRI,*) 'pbc_T_enforce_invsym: Enforcing inversion symmetry in a T tensor.'

  do l = 0, lmax
     pmin = l*(l+1) - l + 1
     pmax = l*(l+1) + l + 1
     do j = 0, lmax
        if (mod(l+j,2) .eq. 1) then
           qmin = j*(j+1) - j + 1
           qmax = j*(j+1) + j + 1

           T(pmin:pmax,qmin:qmax) = 0.0_realk
        end if
     end do
  end do
end subroutine pbc_T_enforce_invsym

subroutine pbc_get_nfsize(n1,n2,n3,layer,lupri)
  implicit none

  integer, intent(inout) :: n1,n2,n3
  integer,intent(in) :: layer,lupri
  integer :: fdim(3)

  call pbcstruct_get_active_dims(fdim)
!
  if (layer .ge. 1 .and. &
      & layer .le. 25) then
     n1 = layer * fdim(1)
     n2 = layer * fdim(2)
     n3 = layer * fdim(3)
  else
     ! put code to determine NF size automatically here;
     ! for now we just invent a value
     n1 = 5 * fdim(1)
     n2 = 5 * fdim(2)
     n3 = 5 * fdim(3)
     call lsquit('Unreasonable value of realspc_range%coulomb_nf.',lupri)
  end if
end subroutine pbc_get_nfsize

! **************************************************************
! **************************************************************
! **                 D E B U G   T O O L S                    **
! **      functions written only for debugging purposes       **
! **************************************************************
! **************************************************************
              
subroutine pbc_Tmatrix_print(T,siz,Tname,lupri)
  implicit none
  
  integer, intent(in) :: siz,lupri
  real(realk), intent(in) :: T(siz,siz)
  character*(*), intent(in) :: Tname
  
  integer :: l, m, row,lunit
  
#ifdef DEBUGPBC
  lunit=-1
  call LSOPEN(lunit,'Tlattice','UNKNOWN','FORMATTED')
  !write(LUPRI,*) 'The first elements of ',Tname,':'
  !write(LUPRI,*) '        j = 0      j = 1                            j = 2'
  write(lunit,*) 'The first elements of ',Tname,':'
  write(lunit,*) '        j = 0      j = 1                            j = 2'
  do l = 0,3
     do m = -l, l
        row = l*(l+1)+1 + m
        if (m .eq. -l) then
           write(lunit,*) T(row,1:siz)
        else
           write(lunit,*) T(row,1:siz)
        end if
     end do
  end do
  call LSCLOSE(lunit,'KEEP')
#else
  write(LUPRI,*) 'The first elements of ',Tname,':'
  write(LUPRI,*) '        j = 0      j = 1                            j = 2'
  do l = 0,3
     do m = -l, l
        row = l*(l+1)+1 + m
        if (m .eq. -l) then
           write(LUPRI,'(A,I2,A,9D11.3)') 'l =',l,': ',T(row,1:9)
        else
           write(LUPRI,'(A,9D11.3)')      '     : ',T(row,1:9)
        end if
     end do
  end do
#endif

end subroutine pbc_Tmatrix_print




SUBROUTINE pbc_ff_fck(Tlmax,tlat,lmax,nbast,ll,nfdensity,nucmom,&
                         & g_2,E_ff,E_nn,lupri)
IMPLICIT NONE
INTEGER,INTENT(IN) :: Tlmax,lmax,nbast,lupri
real(realk), intent(in) :: Tlat((1+lmax)**2,(1+lmax)**2),nucmom((1+lmax)**2)
real(realk),INTENT(INOUT) :: E_ff,E_nn
TYPE(lvec_list_t), INTENT(INOUT) :: ll
TYPE(matrix) :: nfdensity(size(ll%lvec))
TYPE(matrix),intent(inout) :: g_2(size(ll%lvec))
!LOCAL
CHARACTER(len=40) :: filename
CHARACTER(len=10) :: numstr0,numstr1,numstr2,numstr3
!TYPE(lattice_cell_info_t), allocatable :: sphermom(:)
TYPE(lattice_cell_info_t), pointer :: sphermom(:)
INTEGER :: num_latvec,nf,nk,nrlm,j,s,t,st
INTEGER :: y,x2,y2,z2,jk,lm,delta,m,ii,iunit
integer :: cd,il1,il2,il3,indred
INTEGER(short) :: gab1
!TYPE(matrix) :: rhojk((lmax+1)**2)
TYPE(matrix) :: farfieldtmp
REAL(realk) :: phase1,phase2,phase3
REAL(realk),pointer :: rhojk(:)
!REAL(realk) :: mmfck(nfsze,nbast*nbast),kvec(3)!,debug_mat(nbast,nbast)
REAL(realk) :: Coulombf2,Coulombfst,Coulomb2,Coulombst
!REAL(realk) :: PI=3.14159265358979323846_realk
!real(realk) :: tlatlm((lmax+1)**2),tlatlmnu((lmax+1)**2)
!real(realk) :: tlatlm(:),tlatlmnu(:)
real(realk),pointer :: tlatlm(:),tlatlmnu(:)
real(realk) :: debugsumt
complex(complexk) :: phase
character(len=12) :: diis,stiter
!TYPE(matrix) :: debug_tm !FOR debug only
Integer, save :: iter=0
real(realk),external :: ddot
#ifdef DEBUGPBC
real(realk) :: normsq
#endif

iter=iter+1

if(iter .gt. 1) then
  write(stiter,'(I5)') iter-1
  stiter=adjustl(stiter)
  diis='diis_'//trim(stiter)//'_'
endif

call mat_init(farfieldtmp,nbast,nbast)
num_latvec=size(ll%lvec)
!call pbc_get_nfsize(n1,n2,n3,ll%nneighbour,lupri)
nrlm=(lmax+1)**2
call mem_alloc(tlatlm,nrlm)
call mem_alloc(tlatlmnu,nrlm)
call mem_alloc(rhojk,(1+lmax)**2)
tlatlm=0.0_realk
tlatlmnu=0.0_realk

E_ff=0._realk
E_nn=0._realk

call mem_alloc(sphermom,num_latvec)


ii=0
DO nk=1,num_latvec
   x2=int(ll%lvec(nk)%lat_coord(1))
   y2=int(ll%lvec(nk)%lat_coord(2))
   z2=int(ll%lvec(nk)%lat_coord(3))
   gab1=ll%lvec(nk)%maxgab
   if(gab1 .ge. ll%realthr) then
     sphermom(nk)%is_defined = .true.
     call read_multipole_files(ll,lmax,sphermom(nk),nbast,nk,lupri)
     call pbc_multipl_moment_matorder(sphermom(nk)%getmultipole,lmax,nbast)
     call pbc_mat_redefine_q(sphermom(nk)%getmultipole,lmax,nbast)
   endif
enddo


rhojk=0.0_realk
!write(lupri,*) 'Tlattice contracted with just nuclear moments'
DO jk=1,(lmax+1)**2
   DO nk=1,num_latvec
      x2=ll%lvec(nk)%lat_coord(1)
      y2=ll%lvec(nk)%lat_coord(2)
      z2=ll%lvec(nk)%lat_coord(3)

      if(sphermom(nk)%is_defined) then
        if(nfdensity(nk)%init_magic_tag.EQ.mat_init_magic_value) THEN
          rhojk(jk)=rhojk(jk)-&
                 & mat_dotproduct(nfdensity(nk),sphermom(nk)%getmultipole(jk))
        endif
      endif
   ENDDO
ENDDO

!call pbc_redefine_q(rhojk,lmax)
!call pbc_multipl_moment_order(rhojk,lmax)
!#ifdef DEBUGPBC
debugsumT=0.0_realk
!write(lupri,*) 'electronic moments'
do jk=1,(lmax+1)**2
   debugsumt=debugsumt+(rhojk(jk)+nucmom(jk))**2
enddo
write(lupri,*) 'debugsum cell_mom', debugsumt
write(*,*) 'debugsumt cell_mom', debugsumt
write(lupri,*) 'Charge rhojk(1)', rhojk(1)
write(*,*) 'Charge rhojk(1)', rhojk(1)
write(lupri,*) 'Charge nucmom(1)', nucmom(1)
write(*,*) 'Charge nucmom(1)', nucmom(1)
!#endif

DO lm=1,(lmax+1)**2
  tlatlm(lm)=dot_product(tlat(lm,1:nrlm),nucmom+rhojk)
  tlatlmnu(lm)=dot_product(tlat(lm,1:nrlm),nucmom)
ENDDO

#ifdef DEBUGPBC
debugsumT=0.0_realk
do jk=1,(lmax+1)**2
   debugsumt=debugsumt+tlatlm(jk)**2
!   write(lupri,*) rhojk(jk)+nucmom(jk)
enddo
write(lupri,*) 'debugsum tlatlm^2', debugsumt
write(*,*) 'debugsumt tlatlm^2', debugsumt
#endif
write(lupri,*) 'Charge tlatlm(1)', tlatlm(1)
write(*,*) 'Charge tlatlm(1)', tlatlm(1)

!
!write(lupri,*) 'total moments'

#ifdef DEBUGPBC
debugsumT=0.0_realk
do jk=1,(lmax+1)**2
   debugsumt=debugsumt+(rhojk(jk)+nucmom(jk))**2
!   write(lupri,*) rhojk(jk)+nucmom(jk)
enddo
write(lupri,*) 'debugsum rhojk^2', debugsumt
write(*,*) 'debugsumt rhojk^2', debugsumt
#endif
!
#ifdef DEBUGPBC
write(*,*) 'DEBUG 1'
#endif

DO nk=1,num_latvec

   x2=ll%lvec(nk)%lat_coord(1)
   y2=ll%lvec(nk)%lat_coord(2)
   z2=ll%lvec(nk)%lat_coord(3)
   !if(abs(x2) .gt. ll%col1) CYCLE
   !if(abs(y2) .gt. ll%col2) CYCLE
   !if(abs(z2) .gt. ll%col3) CYCLE
   if(.not. ll%lvec(nk)%is_redundant) then
     call find_latt_vectors(nk,x2,y2,z2,ll)
     call find_latt_index(indred,-x2,-y2,-z2,ll,ll%max_layer)
     if(.not.sphermom(nk)%is_defined) CYCLE
#ifdef DEBUGPBC
     write(*,*) 'x2',x2,nk
     call mat_init(ll%lvec(nk)%oper(1),nbast,nbast)
     call mat_zero(ll%lvec(nk)%oper(1))
#endif
     call mat_zero(farfieldtmp)


     DO lm=1,(lmax+1)**2


         call mat_daxpy(-tlatlm(lm),sphermom(nk)%getmultipole(lm),farfieldtmp)
     
              
     
     ENDDO !lm
     !ENDDO !delta
!{{{    
#ifdef DEBUGPBC 
     call mat_copy(1.0_realk,farfieldtmp,ll%lvec(nk)%oper(1))
     write(lupri,*) ''
     write(lupri,*) 'Total Far-field contribution to fock for layer ' ,x2
     write(lupri,*) 'max of rhojk = ', maxval(rhojk,(lmax+1)**2)
     write(lupri,*) 'max of tlatlm = ', maxval(tlatlm,(lmax+1)**2)
     write(lupri,*) 'max of tlatlmnu = ', maxval(tlatlmnu,(lmax+1)**2)
     write(lupri,*) 'max of nucmom = ', maxval(nucmom,(lmax+1)**2)
     call mat_print(farfieldtmp,1,nbast,1,nbast,lupri)
#endif 
!}}}


     if(.not. ll%store_mats) then

        if(ll%lvec(nk)%oper(2)%init_magic_tag.EQ.mat_init_magic_value) then
          !should be with V_z-el as well
          call mat_daxpy(1._realk,farfieldtmp,g_2(nk))
        endif
!{{{
#ifdef DEBUGPBC
          write(lupri,*) 'Near field coul',nk
          call mat_print(ll%lvec(nk)%oper(2),1,nbast,1,nbast,lupri)
          write(lupri,*) 'Total coul',nk
          call mat_print(g_2(nk),1,nbast,1,nbast,lupri)
#endif
!}}}
        !endif
     else
        if(ll%lvec(nk)%oper(2)%init_magic_tag.EQ.mat_init_magic_value) then
          call mat_daxpy(1._realk,farfieldtmp,ll%lvec(nk)%oper(2))
          call pbc_get_file_and_write(ll,nbast,nbast,nk,4,2,'            ')! 4 and 2 Coul J

     !for debugging only
     !call mat_to_full(debug_tm,1.0_realk,debug_mat)
     !write(lupri,*) 'far-field',x2,y2,z2
     !call pbc_matlab_print(debug_mat,nbast,nbast,'far-field',lupri)
         endif
       endif
     !call pbc_get_file_and_write(ll,nbast,nbast,nk,9,1,'            ')! 4 and 2 Coul J
    ! call pbc_get_file_and_write(ll,nbast,nbast,nk,7,3) !7 and 3 fock matrix
!   else
!     call pbc_get_file_and_write(ll,nbast,nbast,nk,4,2,diis)! 4 and 2 Coul J
!    !call pbc_get_file_and_write(ll,nbast,nbast,nk,7,3,diis)!7 and 3 fock matrix
!   endif
!{{{
#ifdef DEBUGPBC
   write(lupri,'(A35)') 'DEBUGPBC farfield and total Coulomb'
   il1=int(ll%lvec(nk)%lat_coord(1))
   il2=int(ll%lvec(nk)%lat_coord(2))
   il3=int(ll%lvec(nk)%lat_coord(3))
   st=0
   Coulombf2=0._realk
   Coulombfst=0._realk
   Coulomb2=0._realk
   Coulombst=0._realk
    DO s=1,nbast
     DO t=1,nbast
      st=st+1
      Coulombf2=Coulombf2+ll%lvec(nk)%oper(1)%elms(st)**2
      Coulombfst=Coulombfst+st*ll%lvec(nk)%oper(1)%elms(st)**2
      Coulomb2=Coulomb2+ll%lvec(nk)%oper(2)%elms(st)**2
      Coulombst=Coulombst+st*ll%lvec(nk)%oper(2)%elms(st)**2
     ENDDO
    ENDDO
    write(lupri,'(A10,X,I3)') 'iteration:',iter
    write(lupri,'(A9,X,I3,x,I3,X,I3)') 't1 t2 t3:',il1,il2,il3
    write(lupri,'(A17)') 'DEBUGPBC farfield'
    write(lupri,'(A4,X,E16.8)') 'Jf^2:',coulombf2
    write(lupri,'(A10,X,E16.8)') 'sum ijJf^2:',coulombfst
    write(lupri,'(A22)') 'DEBUGPBC total coulomb'
    write(lupri,'(A4,X,E16.8)') 'Jt^2:',coulomb2
    write(lupri,'(A10,X,E16.8)') 'sum ijJt^2:',coulombst
#endif
!}}}

 if(ll%lvec(nk)%oper(2)%init_magic_tag.EQ.mat_init_magic_value) then
   call mat_free(ll%lvec(nk)%oper(2))
 endif
!{{{
#ifdef DEBUGPBC
 if(ll%lvec(nk)%oper(1)%init_magic_tag.EQ.mat_init_magic_value) then
   call mat_free(ll%lvec(nk)%oper(1))
 endif
#endif
!}}}
 endif!is_redundant
   !call mat_free(ll%lvec(nk)%oper(3))
ENDDO
!{{{
IF(ll%compare_elmnts ) THEN
  !compare integrals with the old pbc code

!  write(lupri,*) 'comparing multipole moments with old pbc code'
  
  iunit=-1
  DO nf=1,num_latvec
    x2=ll%lvec(nf)%lat_coord(1)
    y2=ll%lvec(nf)%lat_coord(2)
    z2=ll%lvec(nf)%lat_coord(3)
    call find_latt_index(nk,x2,y2,z2,ll,ll%max_layer)
    write(numstr0,'(I5)')  iter
    write(numstr1,'(I5)')  x2
    write(numstr2,'(I5)')  y2
    write(numstr3,'(I5)')  z2
    numstr0=adjustl(numstr0)
    numstr1=adjustl(numstr1)
    numstr2=adjustl(numstr2)
    numstr3=adjustl(numstr3)
    filename='minmom'//trim(numstr0)//trim(numstr1)//trim(numstr2)//trim(numstr3)//'.dat'
    CALL lsOPEN(IUNIT,filename,'unknown','FORMATTED')
!    DO j=1,nbast
!       write(iunit,*) (mmfck(nf,delta+(j-1)*nbast),delta=1,nbast)
!    ENDDO
    call lsclose(iunit,'KEEP')
    if(abs(x2) .gt. ll%nneighbour) CYCLE
    if(abs(y2) .gt. ll%nneighbour) CYCLE
    if(abs(z2) .gt. ll%nneighbour) CYCLE
    filename='Jtot'//trim(numstr0)//trim(numstr1)//trim(numstr2)//trim(numstr3)//'.dat'
    CALL lsOPEN(IUNIT,filename,'unknown','FORMATTED')
    !DO j=1,nbast
    !   write(iunit,*) (ll%lvec(nk)%fck_vec(delta+(j-1)*nbast),delta=1,nbast)
    !ENDDO
    call mat_print(g_2(nk),1,g_2(nk)%nrow,1,g_2(nk)%ncol,iunit)
    call lsclose(iunit,'KEEP')
  
  ENDDO
ENDIF !compare 
!}}}


!E_ff= ddot(nrlm,rhojk,nrlm,tlatlm,nrlm)
E_ff= dot_product(rhojk+nucmom,tlatlm)
E_ff=E_ff/2.0_realk
!E_nn=ddot(nrlm,nucmom,nrlm,tlatlmnu,nrlm)
E_nn= dot_product(nucmom,tlatlmnu)

#ifdef DEBUGPBC
write(lupri,*) 'Farfield energy', E_ff
write(*,*) 'Farfield energy', E_ff
#endif

!write(*,*) 'E_ff from mmit',E_ff,E_nn

call mem_dealloc(tlatlm)
call mem_dealloc(tlatlmnu)
call mem_dealloc(rhojk)
DO ii=1,num_latvec
  if(sphermom(ii)%is_defined) then
    DO y=1,(lmax+1)**2
      call mat_free(sphermom(ii)%getmultipole(y))
    ENDDO
    call mem_dealloc(sphermom(ii)%getmultipole)
  endif
ENDDO
call mem_dealloc(sphermom)
call mat_free(farfieldtmp)

END SUBROUTINE pbc_ff_fck

SUBROUTINE pbc_comp_nucmom(refcell,nucmom,lmax,lupri)

IMPLICIT NONE

  TYPE(moleculeinfo) :: refcell
  integer,intent(in) :: lupri
  integer, intent(in) :: lmax
  real(realk),intent(INOUT) :: nucmom((1+lmax)**2)
  !Local variables
  TYPE(lattice_cell_info_t), pointer :: sphermom(:)
  TYPE(matrix),pointer :: carnucmom(:)
  TYPE(matrix),pointer :: sphnucmom(:)
  logical :: reset_flag
  real(realk) :: q,pos(3)
  real(realk),pointer :: nucmomtmp(:)
  integer :: i,ncarmom,nsphmom
  integer :: iprint

ncarmom = (lmax+1)*(lmax+2)*(lmax+3)/6
!ncarmom = (lmax+1)**2
nsphmom = (lmax+1)**2
!write(*,*) 'debug 1'
!call mem_alloc(sphermom,nfsze)
!write(*,*) 'debug 2'
call mem_alloc(carnucmom,ncarmom)
!write(*,*) 'debug 3'
call mem_alloc(nucmomtmp,ncarmom)
call mem_alloc(sphnucmom,(1+lmax)**2)
!allocate(nucmom(nsphmom))
!write(*,*) 'debug 4'

 
!call Mat_init(mm_dens_tmp,nbast,nbast)
!call Mat_zero(mm_dens_tmp)

DO i=1,ncarmom
   call Mat_init(carnucmom(i),1,1)
   call Mat_zero(carnucmom(i))
enddo

DO i=1,(1+lmax)**2
   call Mat_init(sphnucmom(i),1,1)
   call Mat_zero(sphnucmom(i))
ENDDO


reset_flag = .true.
do i=1,refcell%natoms
   q=refcell%atom(i)%charge
   pos(1)=refcell%atom(i)%center(1)
   pos(2)=refcell%atom(i)%center(2)
   pos(3)=refcell%atom(i)%center(3)
   call calc_pcharge_mom(nucmomtmp,q,pos,lmax,reset_flag,ncarmom,lupri)
   reset_flag = .false.
enddo

!write(*,*) 'debug 3'
iprint=2
DO i=1,ncarmom
   carnucmom(i)%elms(1)=nucmomtmp(i)
ENDDO
!write(lupri,*) 'debug: Charge from Cartesian nuc.mom. is ',&
!                & nucmomtmp(1),carnucmom(1)%elms
!write(*,*) 'debug: Charge from Cartesian nuc.mom. is ',&
!            & nucmomtmp(1),carnucmom(1)%elms

call II_carmom_to_shermom(sphnucmom,carnucmom,(lmax+1)**2,ncarmom,lmax,&
                          & lupri,iprint)

!write(lupri,*) 'debug: Charge from Spherical nuc.mom. is ',&
!                & sphnucmom(1)%elms
!write(*,*) 'debug: Charge from Spherical nuc.mom. is ',&
!            & sphnucmom(1)%elms

!write(*,*) 'DEBUG 1 segmentation fault',ncarmom

DO i=1,ncarmom
   call Mat_free(carnucmom(i))
   !write(*,*) 'DEBUG segmentation fault',i+1,ncarmom
enddo
!write(*,*) 'debug free carnucmom '
call mem_dealloc(carnucmom)
!write(*,*) 'debug dealloc carnucmom '

call mem_dealloc(nucmomtmp)
!write(*,*) 'debug dealloc nucmomtmp '
!write(lupri,*) 'nucmom before transformation'
DO i=1,(lmax+1)**2
  !write(*,*) 'nucmom(i)',i
  nucmom(i) = sphnucmom(i)%elms(1)
  !write(*,*) 'nucmom done '
!  write(lupri,*) nucmom(i)
ENDDO
!write(*,*) 'debug nucmom(i)=sphnucmom '

DO i=1,nsphmom
   call Mat_free(sphnucmom(i))
enddo
!write(*,*) 'debug free sphnucmom '
!deallocate(sphnucmom)
call mem_dealloc(sphnucmom)
!write(*,*) 'debug dealloc sphnucmom '

call pbc_multipl_moment_order(nucmom,lmax)
!write(*,*) 'debug multipl_moment_order'
call pbc_redefine_q(nucmom,lmax)
write(lupri,*) 'charge nuclei',nucmom(1)
write(*,*) 'charge nuclei',nucmom(1)

!call mem_dealloc(sphermom)
!write(lupri,*) 'nucmom after transformation'
!DO i=1,(lmax+1)**2
!  write(lupri,*) nucmom(i)
!ENDDO


END SUBROUTINE pbc_comp_nucmom


SUBROUTINE calc_pcharge_mom(cmm,q,pos,lmax,reset_flag,ncarmom,lupri)
implicit none
      integer,intent(IN) :: lupri,ncarmom
      real(realk),intent(IN) :: q, pos(3)
      real(realk),intent(INOUT) :: cmm(ncarmom)
      integer,intent(IN) :: lmax
      integer, parameter :: max_qlmax_pbc=21
      logical :: reset_flag
      real(realk) :: pos_pow(3,max_qlmax_pbc+1)
      integer :: i,j,k,ltot,index_ijk,icomp

      write(LUPRI,'(A,F10.4,A,3F10.4,A)') 'calc_pcharge_mom: q = ',q,&
     &     ' at pos = (',pos(1),pos(2),pos(3),')'
      !call flshfo(LUPRI)


      if (lmax .lt. 0) then
         call lsquit('calc_pcharge_mom: negative lmax.',lupri)
      else if (lmax .gt. max_qlmax_pbc) then
         call lsquit('calc_pcharge_mom: lmax is too large.',lupri)
      end if

!     Compute x^n, y^n, z^n for n = 0,...,lmax
      pos_pow(1,1) = 1.0_realk
      pos_pow(2,1) = 1.0_realk
      pos_pow(3,1) = 1.0_realk
      do i = 1,lmax
         do icomp = 1,3
            pos_pow(icomp,i+1) = pos_pow(icomp,i) * pos(icomp)
         end do
      end do

!     The components are generated in the order (0 0 0), (1 0 0),
!     (0 1 0), (0 0 1), (2 0 0), (1 1 0), (1 0 1), (0 2 0), ...

!     Loop over Cartesian total "ang. mom."
      index_ijk = 1
      do ltot = 0, lmax, 1
         do i = ltot, 0, -1
            do j = ltot - i, 0, -1
               k = ltot - i - j
               if (reset_flag) then
                  cmm(index_ijk) = q * pos_pow(1,i+1)&
    &                 * pos_pow(2,j+1) * pos_pow(3,k+1)
               else
                  cmm(index_ijk) = cmm(index_ijk) + q * pos_pow(1,i+1)&
     &                 * pos_pow(2,j+1) * pos_pow(3,k+1)
               end if
               index_ijk = index_ijk + 1
            end do
         end do
      end do
      !write(*,*) 'index_ijk = ' ,index_ijk

      return
END SUBROUTINE calc_pcharge_mom


SUBROUTINE pbc_multipl_moment_order(cell_mulmom,lmax)
implicit none
INTEGER,INTENT(IN) :: lmax
real(realk),intent(INOUT) :: cell_mulmom((lmax+1)**2)
!LOCAL variables
INTEGER :: je,M,KHE,KHE0,jmindex
REAL(realk) :: tmp_mulmom((lmax+1)**2)
REAL(realk) :: sphint


DO je=1, (lmax+1)**2
   tmp_mulmom(je)=cell_mulmom(je)
   cell_mulmom(je)=0.0_realk
enddo
KHE0=0
jmindex=1
DO JE=0,LMAX
            DO M=-JE,JE
!              columns are ordered as 0, 1, -1, 2, -2, etc.
               IF(M.GT.0) THEN
                  KHE = KHE0 + 2*M
               ELSE
                  KHE = KHE0 -2*M+1
               END IF
!              The nuclear charges must already be included in the
!              moments, so here we just add up.
               SPHINT = tmp_mulmom(KHE)
!              Add to the cell's moments
               cell_mulmom(jmindex) = SPHINT
               jmindex = jmindex + 1
            END DO
            KHE0 = KHE0 + ( 2*JE+1 ) 
         END DO


END SUBROUTINE pbc_multipl_moment_order

SUBROUTINE pbc_multipl_moment_matorder(sphnucmom,lmax,nbast)
implicit none
INTEGER,INTENT(IN) :: lmax,nbast
!LOCAL variables
INTEGER :: je,M,KHE,KHE0,jmindex,ab
TYPE(matrix),intent(INOUT) :: sphnucmom((1+lmax)**2)
REAL(realk) :: tmp_mulmom(nbast*nbast,(lmax+1)**2)
REAL(realk) :: sphint


DO je=1, (lmax+1)**2
   tmp_mulmom(:,je)= sphnucmom(je)%elms(:)
   call mat_zero(sphnucmom(je))
enddo
KHE0=0
jmindex=1
DO JE=0,LMAX
            DO M=-JE,JE
!              columns are ordered as 0, 1, -1, 2, -2, etc.
               IF(M.GT.0) THEN
                  KHE = KHE0 + 2*M
               ELSE
                  KHE = KHE0 -2*M+1
               END IF
!              The nuclear charges must already be included in the
!              moments, so here we just add up.
               do ab=1,nbast*nbast
                  SPHINT = tmp_mulmom(ab,KHE)
!                 Add to the cell's moments
                  sphnucmom(jmindex)%elms(ab)= SPHINT
               enddo
               jmindex = jmindex + 1
            END DO
            KHE0 = KHE0 + ( 2*JE+1 ) 
         END DO


END SUBROUTINE pbc_multipl_moment_matorder

! This subroutine redefines the multipole moments.
! It is a copy of the FMM code routine mm_renormalise_qlm
! that has been modified slightly.
!        
subroutine pbc_redefine_q(qlm,lmax)
  implicit none     
  integer, intent(in) :: lmax
  real(realk), intent(inout) :: qlm((1+lmax)**2)

  integer :: l,m,p,pp
  real(realk) :: pref, s

  !write(LUPRI,*) 'pbc_redefine_q: Redefining multipole moments, lmax = ',lmax
         
  ! prefactor to symmetrize T-matrix
  do l = 0, lmax
     if (mod(l,2) .eq. 1) then
        s = -1.0_realk
     else
        s = 1.0_realk
     end if
     pp = l*(l+1) +1
     do m = -l, -1
        pref = -1.0_realk / (sqrt(2.0_realk * factorial(l-m)*factorial(l+m)))
        p = pp+m
        qlm(p) = pref*qlm(p) * s
     end do
     pref = 1.0_realk / factorial(l)
     p = pp  ! m = 0
     qlm(p) = pref*qlm(p) * s
     do m = 1, l
        pref = ((-1)**m) / sqrt(2.0_realk * factorial(l-m)*factorial(l+m))
        p = pp+M
        qlm(p) = pref*qlm(p) * s
     end do
  end do

contains

  real(realk) function factorial(n)
    implicit none
    integer, intent(in) :: n
    integer :: i
    factorial = 1
    do i = n, 2, -1
       factorial = factorial*i
    end do
  end function factorial

end subroutine pbc_redefine_q

! This subroutine redefines the multipole moments for a matrix.
! It is a copy of the FMM code routine mm_renormalise_qlm
! that has been modified slightly.
!        
subroutine pbc_mat_redefine_q(qlm,lmax,nbast)
  implicit none     
  integer, intent(in) :: lmax,nbast
  TYPE(matrix), intent(inout) :: qlm((1+lmax)**2)

  integer :: l,m,p,pp,ab
  real(realk) :: pref, s

  !write(LUPRI,*) 'pbc_redefine_q: Redefining multipole moments, lmax = ',lmax
         
  ! prefactor to symmetrize T-matrix
  do l = 0, lmax
     if (mod(l,2) .eq. 1) then
        s = -1.0_realk
     else
        s = 1.0_realk
     end if
     pp = l*(l+1) +1
     do m = -l, -1
        pref = -1.0_realk / (sqrt(2.0_realk * factorial(l-m)*factorial(l+m)))
        p = pp+m
        do ab=1,nbast*nbast
           qlm(p)%elms(ab) = pref*qlm(p)%elms(ab) * s
        enddo
     end do
     pref = 1.0_realk / factorial(l)
     p = pp  ! m = 0
     do ab=1,nbast*nbast
        qlm(p)%elms(ab) = pref*qlm(p)%elms(ab) * s
     enddo
     do m = 1, l
        pref = ((-1)**m) / sqrt(2.0_realk * factorial(l-m)*factorial(l+m))
        p = pp+M
        do ab=1,nbast*nbast
           qlm(p)%elms(ab) = pref*qlm(p)%elms(ab) * s
        enddo
     end do
  end do

contains

  real(realk) function factorial(n)
    implicit none
    integer, intent(in) :: n
    integer :: i
    factorial = 1
    do i = n, 2, -1
       factorial = factorial*i
    end do
  end function factorial

end subroutine pbc_mat_redefine_q


!SUBROUTINE READ_multipole_files(il1,il2,il3,maxmultmom,sphermom,nbast,lupri)
SUBROUTINE READ_multipole_files(lattice,maxmultmom,sphermom,nbast,nk,lupri)
  IMPLICIT NONE
  !TYPE(lvec_list_t), INTENT(IN) :: ll
  !INTEGER, INTENT(IN) :: il1,il2,il3, maxmultmom,nbast,lupri
  INTEGER, INTENT(IN) :: nk, maxmultmom,nbast,lupri
  TYPE(lvec_list_t), INTENT(INOUT) :: lattice
  TYPE(lattice_cell_info_t), INTENT(INOUT) :: sphermom!(:)
  !LOCAL VARIABLES
  INTEGER :: fileunit,totalmom !For file detection
  INTEGER(short) :: gab1
  INTEGER :: i,j,stat,il1,il2,il3
  INTEGER :: sphermoms,ii, debdum
  CHARACTER(len=10) :: numtostring1,numtostring2,numtostring3
  CHARACTER(LEN=20) :: filename
  REAL(realk) :: testread
  !For debugging
  !REAL(realk) :: debugsumT

  fileunit=9950

  
  totalmom=(maxmultmom+1)**2
  !numnn=0
!  DO dummy=1,numvecs

!     call find_latt_vectors(dummy,il1,il2,il3,ll)

    ! if(abs(il1) .gt. n_neighbour) call lsquit('ERROR in READ_multipole_files&
    ! & vector not in NF',lupri)
    ! if(abs(il2) .gt. n_neighbour) call lsquit('ERROR in READ_multipole_files&
    ! & vector not in NF',lupri)
    ! if(abs(il3) .gt. n_neighbour) call lsquit('ERROR in READ_multipole_files&
    ! & vector not in NF',lupri)
!     numnn=numnn+1

!  ENDDO
  !allocate(sphermom(numnn))

  !DO dummy=1,numnn
  !   allocate(sphermom(dummy)%getmultipole(totalmom))
  !   DO ii=1,totalmom
  !   call Mat_init(sphermom(dummy)%getmultipole(ii),nbast,nbast)
  !   call Mat_zero(sphermom(dummy)%getmultipole(ii))
  !   ENDDO
  !ENDDO

    
    il1=int(lattice%lvec(nk)%lat_coord(1))
    il2=int(lattice%lvec(nk)%lat_coord(2))
    il3=int(lattice%lvec(nk)%lat_coord(3))
    call mem_alloc(sphermom%getmultipole,totalmom)
    DO ii=1,totalmom
       call Mat_init(sphermom%getmultipole(ii),nbast,nbast)
       call Mat_zero(sphermom%getmultipole(ii))
    ENDDO

    !call find_latt_vectors(dummy,il1,il2,il3,ll)

    write(numtostring1,'(I5)')  il1
    write(numtostring2,'(I5)')  il2
    write(numtostring3,'(I5)')  il3
    numtostring1=adjustl(numtostring1)
    numtostring2=adjustl(numtostring2)
    numtostring3=adjustl(numtostring3)
    filename='Qcdlm'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
    !write(*,*) 'filnavn for Qcdlm', filename
    !filename2=&
    !&'kopi'//trim(numtostring1)//trim(numtostring2)//trim(numtostring2)//'.dat'

    call LSOPEN(fileunit,filename, 'old','unformatted')
!   read(fileunit,IOSTAT=stat) testread
    !write(*,*) testread,stat
!   BACKSPACE(fileunit)
!   if(stat .ne. 0) then
    rewind(fileunit) 
!   endif
    !OPEN(UNIT=9998,FILE=filename2)
    debdum=0
    !write(lupri,*) 'Debug 1: ', fileunit
    DO sphermoms=1,totalmom
        
       !DO ii=1,nbast*nbast
       debdum=debdum+1

       call mat_read_from_disk(fileunit,sphermom%getmultipole(sphermoms),.true.)
       !READ(fileunit,*) sphermom(numnn)%getmultipole(sphermoms)%elms!(ii)  
       !do j=1,nbast
       !READ(fileunit) (sphermom%getmultipole(sphermoms)%elms(i+(j-1)*nbast),&
       !i=1,nbast)
       !enddo

       !debugsumT=0.0_realk
       !do j=1,nbast*nbast
       !  debugsumT=debugsumT+sphermom%getmultipole(sphermoms)%elms(j)
       !enddo
       !if(debugsumT .gt. 0.0_realk) then
       !  write(*,*) 'debugsumT > 0.0_realk',il1
       !  write(*,*) debugsumT , sphermoms
       !  write(*,*) sphermom%getmultipole(sphermoms)%elms
       !  stop
       !endif
       !write(9998,*) sphermom(dummy)%getmultipole(sphermoms)%elms!(ii)
       !write(*,*) 'Debug 2: ', dummy, ii,nbast*nbast
       !ENDDO
    ENDDO
    !write(*,*) 'Debug 3: ', filename
    !stop
    call LSCLOSE(fileunit,'KEEP') 

  !ENDDO

!  write(*,*) 'Debug 4: ', numnn
!  stop
  !DO ii=1,totalmom
  !   call Mat_free(sphermom%getmultipole(ii))
  !enddo


END SUBROUTINE READ_multipole_files

SUBROUTINE pbc_multipole_expan_k(lupri,luerr,setting,nbast,lattice,refcell,numvecs,maxmultmom)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: lupri,luerr,numvecs,nbast
  INTEGER,intent(in) :: maxmultmom
  TYPE(lvec_list_t),intent(in) ::lattice
  !REAL(realk), INTENT(IN) :: kvec(3)
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(moleculeinfo),intent(inout) :: refcell
  !TYPE(moleculeinfo),intent(in) :: latt_cell(numvecs)

  INTEGER(short)        :: gab1
  Integer :: refindex,idx,il1,il2,il3,nSpherMom,l
  TYPE(matrix),pointer :: spherMom(:)
  CHARACTER(len=20) ::filename
  CHARACTER(len=10) :: numtostring1,numtostring2,numtostring3
  INTEGER :: wunit,I,j
  
  wunit=9910
  nSpherMom=(maxMultmom+1)**2
  call set_lstime_print(.false.)
  
  write(lupri,*) 'PBC_MULTIPOLE, maxMultmom,nspherMom', maxMultmom,nspherMom
  write(lupri,*) 'DEBUG INSIDE PBC_MULTIPOLE'
  
  call mem_Alloc(spherMom,nSpherMom)

  DO l=1,nSpherMom
     call Mat_init(spherMom(l),nbast,nbast)
     call Mat_zero(spherMom(l))
  ENDDO
  call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)

  DO idx=1,numvecs
     

     !if(abs(il1) .gt. abs(lattice%nneighbour)) CYCLE
     !if(abs(il2) .gt. abs(lattice%nneighbour)) CYCLE
     !if(abs(il3) .gt. abs(lattice%nneighbour)) CYCLE
     !input%Basis%regular%atomtype(:)%shell(:)%segment(1)%exponents(:)
     !phase=cmplx(0,k1*il1+k2*il2+k3*il3)

     gab1=lattice%lvec(idx)%maxgab
     if(gab1 .ge. lattice%realthr) then
       !write(*,*) 'is gab1 greater than -12?'

       call find_latt_vectors(idx,il1,il2,il3,lattice)
       write(lupri,*) 'sphermom computed for ', il1
       call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(idx)%molecule,2)
       !call TYPEDEF_setmolecules(setting,refcell,1,refcell,2)

       call II_get_sphmom(LUPRI,LUERR,SETTING,spherMom,nSpherMom,maxMultmom,&
       & 0.0_realk,0.0_realk,0.0_realk)

!      call pbc_readopmat2(il1,il2,il3,matris,nbast,'OVERLAP',.true.,.false.)
!      call write_matrix(matris,nbast,nbast)

       write(numtostring1,'(I5)')  il1
       write(numtostring2,'(I5)')  il2
       write(numtostring3,'(I5)')  il3
                                         !numstring
       numtostring1=adjustl(numtostring1)
       numtostring2=adjustl(numtostring2)
       numtostring3=adjustl(numtostring3)
       filename='Qcdlm'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
       !write(lupri,*) 'filnavn for multipolfiler ',filename,numtostring1 
       call LSOPEN(wunit,filename,'unknown','unformatted')

       !write(lupri,*) 'sphermom 0',il1
       !call mat_print(spherMom(1), 1, nbast, 1,nbast, lupri)

       DO l=1,nSpherMom
         ! DO j=1,nbast
         ! write(wunit) (spherMom(l)%elms(i+(j-1)*nbast),i=1,nbast)
         ! enddo
          call mat_write_to_disk(wunit,spherMom(l),.true.)
          call mat_free(spherMom(l))
          call Mat_init(spherMom(l),nbast,nbast)
          call Mat_zero(spherMom(l))
       ENDDO


       !write(*,*) 'filnavn for multipolfiler skrevet til disk ',trim(numtostring3)

       call LSCLOSE(wunit,'KEEP')
       wunit=wunit+1
     endif

  ENDDO

!  filename='Qcdlm.dat'
!  call LSOPEN(wunit,filename,'unknown','formatted')
  DO l=1,nSpherMom
!     !call mat_print(spherMom(l),1,SpherMom(l)%nrow,1,SpherMom(l)%ncol,wunit)
!     write(wunit,*) spherMomt(l)%elms
     call mat_free(sphermom(l))
!     call mat_free(sphermomt(l))
  ENDDO

!  call LSCLOSE(wunit,'KEEP')
  call mem_dealloc(sphermom)
  call set_lstime_print(.true.)
  !deallocate(sphermomt)
  write(lupri,*) 'Finished PBC_multipole'

END SUBROUTINE pbc_multipole_expan_k


SUBROUTINE pbc_charge_density(lattice)
  IMPLICIT NONE
  TYPE(lvec_list_t),intent(in) ::lattice
  INTEGER :: nlatvecs
  !INTEGER,INTENT(IN) :: lmax
  !lattice%lvec(j)%std_coord(123)

  nlatvecs=size(lattice%lvec)
  

END SUBROUTINE pbc_charge_density

SUBROUTINE pbc_local_expan_k(lupri,luerr,setting,nbast,lattice,refcell,numvecs,maxmultmom)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: lupri,luerr,numvecs,maxmultmom,nbast
  TYPE(lvec_list_t),intent(in) ::lattice
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(moleculeinfo),intent(inout) :: refcell

  Integer :: refindex,idx,il1,il2,il3,nSpherMom,l
  TYPE(matrix),pointer :: spherMom(:)
  CHARACTER(len=20) ::filename
  CHARACTER(len=10) :: numtostring1,numtostring2,numtostring3
  INTEGER :: wunit
  
  wunit=9910
  nSpherMom=0
  call set_lstime_print(.false.)

  DO l=0,maxmultmom
      nSpherMom=nSpherMom+2*l+1
  ENDDO
  
  write(lupri,*) 'PBC_MULTIPOLE, maxMultmom,nspherMom', maxMultmom,nspherMom
  !write(lupri,*) 'DEBUG INSIDE PBC_MULTIPOLE'
  
  call mem_Alloc(spherMom,nSpherMom)
  DO l=1,nSpherMom
     call Mat_init(spherMom(l),nbast,nbast)
     call Mat_zero(spherMom(l))
  ENDDO
  call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)
  DO idx=1,numvecs
     
     !phase1 = kvec(1)*lattice%lvec(idx)%std_coord(1)
     !phase2 = kvec(2)*lattice%lvec(idx)%std_coord(2)
     !phase3 = kvec(3)*lattice%lvec(idx)%std_coord(3)
     !phase=CMPLX(0.,phase1+phase2+phase3)
     !phase=exp(phase)

     call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(idx)%molecule,2)

     call find_latt_vectors(idx,il1,il2,il3,lattice)
     !if(abs(il1) .le. abs(lattice%nneighbour)) CYCLE
     !if(abs(il2) .le. abs(lattice%nneighbour)) CYCLE
     !if(abs(il3) .le. abs(lattice%nneighbour)) CYCLE
     !input%Basis%regular%atomtype(:)%shell(:)%segment(1)%exponents(:)


     call II_get_sphmom(LUPRI,LUERR,SETTING,spherMom,nSpherMom,maxMultmom,0.0_realk,0.0_realk,0.0_realk)

!       call pbc_readopmat2(il1,il2,il3,matris,nbast,'OVERLAP',.true.,.false.)
!       call write_matrix(matris,nbast,nbast)
!
    write(numtostring1,'(I5)')  il1
    write(numtostring2,'(I5)')  il2
    write(numtostring3,'(I5)')  il3
                                      !numstring
    numtostring1=adjustl(numtostring1)
    numtostring2=adjustl(numtostring2)
    numtostring3=adjustl(numtostring3)
    filename='Qcdlm'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
    write(lupri,*) filename
    write(*,*) 'filnavn debug',filename
     call LSOPEN(wunit,filename,'unknown','formatted')

     write(wunit,*) idx
     DO l=1,nSpherMom
     !call mat_print(spherMom(l),1,SpherMom(l)%nrow,1,SpherMom(l)%ncol,wunit)
     write(wunit,*) spherMom(l)%elms
     ENDDO


     call LSCLOSE(wunit,'KEEP')
     !wunit=wunit+1

  END DO

  DO l=1,nSpherMom
     call mat_free(sphermom(l))
  ENDDO
  call mem_dealloc(sphermom)
  call set_lstime_print(.true.)

 

END SUBROUTINE pbc_local_expan_k



END MODULE pbc_ff_contrib
#endif
