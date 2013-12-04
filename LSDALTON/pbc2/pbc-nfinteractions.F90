#ifdef MOD_UNRELEASED
MODULE pbc_interactions
  USE PRECISION
  USE fundamental
  USE TYPEDEF
  USE memory_handling
  USE integralinterfaceMOD
  USE matrix_module
  USE matrix_operations
  USE pbc_MSC
  USE pbc_matrix_operations
  USE lattice_type
  USE lattice_vectors
  USE lstiming
  CONTAINS

!> \author Johannes Rekkedal
!> \date 2013
!> \brief Finds the cutoff distance to the ref. cell for 1-pt. operators.
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param nbast 			Number of basis func.
!> \param lattice 		Information about the lattice.
!> \param latt_cell 		Molecule info lattice cell.
!> \param refcell 		Molecule info reference cell.
SUBROUTINE find_cutoff_onep(lupri,luerr,setting,nbast,latt_cell,refcell,lattice)
  IMPLICIT NONE
  ! input and output arguments
  INTEGER, INTENT(IN) :: lupri,luerr,nbast
  TYPE(lvec_list_t),INTENT(INOUT) :: lattice
  TYPE(moleculeinfo),INTENT(INOUT),TARGET :: latt_cell(size(lattice%lvec))
  TYPE(moleculeinfo),INTENT(INOUT),TARGET :: refcell
  TYPE(lssetting),INTENT(INOUT) :: setting 
  ! local variables
  INTEGER ::  idx

  call set_lstime_print(.false.)

  do idx=1,size(lattice%lvec) 
     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell, &
		  & 3,latt_cell(idx),4)
     call II_get_maxGabelm_ScreenMat(lupri,luerr,setting,nbast, &
		  & lattice%lvec(idx)%maxGab)
  enddo

  call set_lstime_print(.true.)

END SUBROUTINE find_cutoff_onep


!> \author Johannes Rekkedal
!> \date 2013
!> \brief Finds the cutoff distance to the ref. cell for 2-pt. operators.
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param nbast 			Number of basis func.
!> \param lattice 		Information about the lattice.
!> \param latt_cell 		Molecule info lattice cell.
!> \param refcell 		Molecule info reference cell.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix.
SUBROUTINE find_cutoff_twop(lupri,luerr,setting,nbast,lattice, &
		latt_cell,refcell,numvecs,nfdensity)
	IMPLICIT NONE
	! input and output arguments
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs ! nlayer 
	TYPE(moleculeinfo),INTENT(IN) :: latt_cell(numvecs)
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(matrix),INTENT(IN) :: nfdensity(numvecs)
	TYPE(lvec_list_t),INTENT(INOUT) ::lattice
	! local variables
	TYPE(LSSETTING),INTENT(INOUT) :: setting 
	TYPE(MATRIX),POINTER :: Kx(:), K_tmp(:)
	REAL(realk) :: exact_xch
	INTEGER :: j,k,checknf
	INTEGER :: idx,index2,index3,newcell
	INTEGER :: l1,l2,l3
	INTEGER :: l11,l12,l13,il21,il22,il23,il31,il32,il33
	LOGICAL :: ntimes

	call set_lstime_print(.false.)
	write(lupri,*) 'Finding nlayer for Kop',numvecs, lattice%ldef%is_active(1),&
		lattice%max_layer

	call mem_alloc(K_tmp,1)
	call mem_alloc(Kx,1)
	call mat_init(Kx(1),nbast,nbast)
	call mat_zero(Kx(1))
	call mat_init(K_tmp(1),nbast,nbast)
	call mat_zero(K_tmp(1))

	ntimes = .false.

	do l3=0,lattice%max_layer*lattice%fdim(3)
		if(ntimes) exit
		do l2=0,lattice%max_layer*lattice%fdim(2)
			if(ntimes) exit
			do l1=0,lattice%max_layer*lattice%fdim(1)
				if(ntimes) exit
				call find_latt_index(idx,l1,l2,l3,lattice,lattice%max_layer)
				call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),3)
				setting%samemol(1,3)=.false.
				setting%samemol(3,1)=.false.
				do index2=1,numvecs
					call find_latt_vectors(index2,il21,il22,il23,lattice)
					do index3=1,numvecs

						call find_latt_vectors(index3,il31,il32,il33,lattice)
						if(.not. lattice%lvec(index3)%dm_computed) CYCLE
						! ToDo Remove: base on CS screening at some point
						l11=il21+il31
						if( abs(l11) .gt.lattice%max_layer) CYCLE
						! todo l1 blir større enn max slik at newcell blir negativ
						! sjekk dette, er ikke sikker på om l1 skal være slik.
						! Likedan man l2 og l3
						l12=il22+il32
						if( abs(l12) .gt.lattice%max_layer) CYCLE
						l13=il23+il33
						if( abs(l13) .gt.lattice%max_layer) CYCLE
						call find_latt_index(newcell,l11,l12,l13,lattice,lattice%max_layer)
						call find_latt_index(checknf,il31,il32,il33,&
							lattice,lattice%max_layer)
					
						setting%samefrag=.false.
						setting%samemol=.false.

						! ToDo Should be turned on at some point
						call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx), &
							& 3,latt_cell(index2),2,latt_cell(newcell),4)
						call mat_zero(K_tmp(1))
						call II_get_exchange_mat(lupri,luerr,setting, &
							& nfdensity(checknf:checknf),1,.false.,K_tmp)
						call mat_daxpy(0.5_realk,K_tmp(1),Kx(1))
					enddo
				enddo

				exact_xch=0.0_realk
				do j=1,nbast*nbast
					exact_xch=exact_xch+abs(kx(1)%elms(j))
				enddo

				if(exact_xch .lt. 1E-10 .and. .not. ntimes) then
					k = max(l3,l2)
					k = max(k,l1)
					lattice%kx1=k
					lattice%kx2=k
					lattice%kx3=k
					ntimes=.true.
					write(lupri,*) 'nlayer Kx =', k
					write(lupri,*) 'lattice%kx1 =', lattice%kx1
					write(lupri,*) 'lattice%kx2 =', lattice%kx3
					write(lupri,*) 'lattice%kx2 =', lattice%kx3
				endif 

				CALL mat_zero(Kx(1))

				if(ntimes) exit
			end do
		end do
	enddo

	call mat_free(K_tmp(1))
	call mat_free(Kx(1))
	call mem_dealloc(K_tmp)
	call mem_dealloc(Kx)
	call set_lstime_print(.true.)
	write(lupri,*) 'Finished with find_cutoff_two_p'

END SUBROUTINE find_cutoff_twop

!> \author Johannes Rekkedal 
!> \date 2013
!> \brief Finds distance between two real(realk) 3D vectors.
!> \param distance 	Output distance between the two 3D vectors.
!> \param tvec 		input vec 1.
!> \param latstdvec 	input vec 2.
SUBROUTINE calc_distance(distance,tvec,latstdvec)
  IMPLICIT NONE
  REAL(realk), INTENT(IN) :: tvec(3),latstdvec(3)
  REAL(realk), INTENT(OUT) :: distance
  
  distance = (tvec(1)-latstdvec(1))**2 + (tvec(2)-latstdvec(2))**2 &
	  & +(tvec(3)-latstdvec(3)**2)
  distance=sqrt(distance)

END SUBROUTINE calc_distance

!> \author Johannes Rekkedal
!> \date 2013
!> \brief For computing overlap integrals over different cells including for 
!> \brief the reference cell.
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param natoms 			Number of atoms in the refcell.
!> \param nbast 			Number of basis func.
!> \param lattice 		Information about the lattice.
!> \param latt_cell 		Molecule info. Threat lattice cell as molecule.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param ovl 		 		Overlap (array of matrices for each lattice cell.
!todo is it nec to calc entire ovl or can we make use of symmetries?
SUBROUTINE pbc_overlap_k(lupri,luerr,setting,natoms,nbast,lattice, &
		& latt_cell,refcell,numvecs,ovl)
	IMPLICIT NONE
	! input and output arguments
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(lvec_list_t),INTENT(INOUT) :: lattice
	TYPE(lssetting),INTENT(INOUT) :: setting 
	TYPE(moleculeinfo),INTENT(IN) :: latt_cell(numvecs)
	TYPE(matrix),TARGET :: ovl(numvecs)
	! local variables
	INTEGER :: i,j,il1,il2,il3
	INTEGER :: idx,refindex
	INTEGER :: maxl1,maxl2,maxl3
	INTEGER(short) :: gab1
	INTEGER,SAVE :: iter=0
	iter=iter+1

	call set_lstime_print(.false.)
	write(lupri,*) 'Number of lattice vectors ', numvecs

	i=1
	call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)
	maxl1=0
	maxl2=0
	maxl3=0
	DO idx=1,numvecs
		lattice%lvec(idx)%ovl_computed=.false.

		!Doing translations
		call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),3)
		setting%samemol(1,3)=.false.
		setting%samemol(3,1)=.false.

		call find_latt_vectors(idx,il1,il2,il3,lattice)
		!So that we do not consider negligible integrals
		!if(abs(il1) .gt. lattice%nneighbour) CYCLE
		!if(abs(il2) .gt. lattice%nneighbour) CYCLE
		!if(abs(il3) .gt. lattice%nneighbour) CYCLE
		gab1=lattice%lvec(idx)%maxgab
		if(gab1 .ge. -12) then
			maxl1=max(abs(il1),maxl1)
			maxl2=max(abs(il2),maxl2)
			maxl3=max(abs(il3),maxl3)
			lattice%lvec(idx)%ovl_computed=.true.

			if(.not. lattice%store_mats) then
				call mat_init(ovl(idx),nbast,nbast)
			endif
			call mat_init(lattice%lvec(idx)%oper(1),nbast,nbast)
			call mat_zero(lattice%lvec(idx)%oper(1))
			! get the overlap matrix for cell between reference and cell l
			call II_get_overlap(lupri,luerr,setting,lattice%lvec(idx)%oper(1))
			if(lattice%store_mats) then !todo necc? are the matrices ever stored??
				call pbc_get_file_and_write(lattice,nbast,nbast,idx,1,1, &
					& '            ')!1 refers to overlap
			else
				call mat_copy(1.0_realk,lattice%lvec(idx)%oper(1),ovl(idx))
			end if

			call mat_free(lattice%lvec(idx)%oper(1))
		endif
	end do
	lattice%oneop1=maxl1
	lattice%oneop2=maxl2
	lattice%oneop3=maxl3

	call set_lstime_print(.true.)
	write(lupri,*) 'finished pbc_overlap_k'

END SUBROUTINE pbc_overlap_k

!> \author Johannes Rekkedal
!> \date 2013
!> \brief For computing the kinetic part of fock matrix.
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param natoms 			Number of atoms in the refcell.
!> \param nbast 			Number of basis func.
!> \param lattice 		Information about the lattice.
!> \param latt_cell 		Molecule info. Threat lattice cell as molecule.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix
!> \param f1 		 		Fock matrix. Initialized with kinetic part here.
!> \param e_kin 			Kinetic energy (optional, calculated if present)
SUBROUTINE pbc_kinetic_k(lupri,luerr,setting,natoms,nbast,lattice, &
		& latt_cell,refcell,numvecs,nfdensity,f_1,e_kin)
	IMPLICIT NONE
	! input and output variables
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(moleculeinfo),DIMENSION(numvecs),INTENT(IN) :: latt_cell
	TYPE(LSSETTING),INTENT(INOUT) :: setting 
	TYPE(lvec_list_t),INTENT(INOUT) ::lattice
	TYPE(matrix),INTENT(IN) :: nfdensity(numvecs)
	TYPE(matrix),INTENT(INOUT),TARGET :: f_1(numvecs)
	REAL(realk),INTENT(OUT),OPTIONAL :: e_kin
	! local variables
	LOGICAL :: calc_e_kin
	INTEGER :: il1,il2,il3
	INTEGER :: idx,refindex,indred
	INTEGER(short) :: gab1
	INTEGER,SAVE :: iter=0

	write(lupri,*) 'Starting pbc_kinetic_k'
	call set_lstime_print(.false.)
	iter=iter+1
	! print latt index
	call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)
	write(lupri,*) 'reference cell index=', refindex,lattice%max_layer

	calc_e_kin=.false.
	if (PRESENT(e_kin)) then
		e_kin=0._realk
		calc_e_kin=.true.
	endif

	do idx=1,numvecs

		lattice%lvec(idx)%f1_computed=.false.
		if(.not. lattice%lvec(idx)%is_redundant) then
			call find_latt_vectors(idx,il1,il2,il3,lattice)
			call find_latt_index(indred,-il1,-il2,-il3,lattice,lattice%max_layer)
			if(il1**2+il2**2+il3**2 .gt. 0) lattice%lvec(indred)%is_redundant =.true.
			call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),3)
			setting%samemol(1,3)=.false.
			setting%samemol(3,1)=.false.
			gab1=lattice%lvec(idx)%maxgab

			if(gab1 .ge. -12) then
				lattice%lvec(idx)%f1_computed=.true.
				lattice%lvec(indred)%f1_computed=.true.
				call mat_init(lattice%lvec(idx)%oper(1),nbast,nbast)
				call mat_zero(lattice%lvec(idx)%oper(1))
				call mat_init(f_1(idx),nbast,nbast)
				call mat_zero(f_1(idx))
				call II_get_kinetic(lupri,luerr,setting,lattice%lvec(idx)%oper(1))
				call mat_copy(1._realk,lattice%lvec(idx)%oper(1),f_1(idx))
				call pbc_get_file_and_write(lattice,nbast,nbast,idx,2,1 &
					& ,'            ')!2 refers to kin

				if (calc_e_kin) then
					if(nfdensity(idx)%init_magic_tag.eq.mat_init_magic_value)then
						if(il1**2+il2**2+il3**2 .gt. 0) then
							E_kin=E_kin+2.*mat_dotproduct( &
								& lattice%lvec(idx)%oper(1),nfdensity(idx))
						else
							E_kin=E_kin+mat_dotproduct( &
								& lattice%lvec(idx)%oper(1),nfdensity(idx))
						endif
					end if
				end if

				write(lupri,*) 'Kinetic mat finished for', il1,il2,il3
				call mat_free(lattice%lvec(idx)%oper(1))
			endif ! maxgab
		endif !is_redundant
	end do

	call set_lstime_print(.true.)
	write(lupri,*) 'Finished pbc_kinetic_k'

END SUBROUTINE pbc_kinetic_k

!> \author Johannes Rekkedal
!> \date 2013
!> \brief For computing nuclear attraction part of the fock matrix 
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param natoms 			Number of atoms in the refcell.
!> \param nbast 			Number of basis func.
!> \param lattice 		Information about the lattice.
!> \param latt_cell 		Molecule info. Threat lattice cell as molecule.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix
!> \param f1 		 		Fock matrix. Nuc part added here.
!> \param e_en 			Kinetic energy (optional, calculated if present)
SUBROUTINE pbc_nucattrc_k(lupri,luerr,setting,natoms,nbast,lattice, &
		& latt_cell,refcell,numvecs,nfdensity,f_1,e_en)
	IMPLICIT NONE
	! input and output variables
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(lssetting),INTENT(INOUT) :: setting 
	TYPE(lvec_list_t),intent(INOUT) :: lattice
	TYPE(moleculeinfo), DIMENSION(numvecs),INTENT(IN) :: latt_cell
	TYPE(matrix),INTENT(IN) :: Nfdensity(numvecs)
	TYPE(matrix),INTENT(INOUT),TARGET :: f_1(numvecs)
	TYPE(matrix) :: H,H1
	REAL(realk),INTENT(OUT),OPTIONAL :: E_en
	! local variables
	INTEGER :: j,k,il11,il12,il13,il21,il22,il23
	INTEGER :: index1,index2,refindex
	INTEGER :: iunit,indred
	INTEGER(short) :: gab1
	INTEGER,SAVE :: iter=0
	CHARACTER(LEN=10) :: numstr1,numstr2,numstr3
	CHARACTER(LEN=40) :: filename
	LOGICAL :: calc_e_en

	call set_lstime_print(.false.)

	iter=iter+1

	!call mem_alloc
	call mat_init(H,nbast,nbast)
	call mat_zero(H)
	call mat_init(H1,nbast,nbast)
	call mat_zero(H1)

	call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)
	write(lupri,*) 'reference cell index=', refindex

	calc_e_en=.false.
	if (PRESENT(E_en)) then
		E_en=0.0_realk
		calc_e_en=.true.
	end if

	do index1=1,numvecs

		call find_latt_vectors(index1,il11,il12,il13,lattice)
		!if(abs(il11) .gt. lattice%nneighbour) CYCLE
		!if(abs(il12) .gt. lattice%nneighbour) CYCLE
		!if(abs(il13) .gt. lattice%nneighbour) CYCLE

		lattice%lvec(index1)%Vz_computed=.false.
		if(.not. lattice%lvec(index1)%is_redundant) then
			call find_latt_index(indred,-il11,-il12,-il13,lattice, &
				& lattice%max_layer)
			if(il11**2+il12**2+il13**2 .gt. 0) then
				lattice%lvec(indred)%is_redundant =.true.
			end if
			gab1=lattice%lvec(index1)%maxgab

			if(gab1 .ge. -12) then	
				lattice%lvec(index1)%Vz_computed=.true.
				lattice%lvec(indred)%Vz_computed=.true.
				call mat_init(lattice%lvec(index1)%oper(1),nbast,nbast)
				call mat_zero(lattice%lvec(index1)%oper(1))
				do index2=1,numvecs
					call find_latt_vectors(index2,il21,il22,il23,lattice)
					if(abs(il21) .gt. lattice%nf) CYCLE 
					if(abs(il22) .gt. lattice%nf) CYCLE
					if(abs(il23) .gt. lattice%nf) CYCLE

					call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index1), &
						& 2,latt_cell(index2),3)

					!setting%samemol(1,3)=.false.
					!setting%samemol(3,1)=.false.
					!setting%samemol(1,2)=.false.
					!setting%samemol(2,1)=.false.
					!setting%samemol(3,2)=.false.
					!setting%samemol(2,3)=.false.
					!These two screenings should be turned on at some point

					call II_get_nucel_mat(lupri,luerr,setting,H1)
					call mat_daxpy(1._realk,H1,lattice%lvec(index1)%oper(1))

				enddo

				call mat_daxpy(1._realk,lattice%lvec(index1)%oper(1),f_1(index1))

				if (calc_e_en) then
					if(nfdensity(index1)%init_magic_tag.eq.mat_init_magic_value)then
						if(il11**2+il12**2+il13**2 .gt. 0) then
							E_en=E_en+2.*mat_dotproduct(lattice%lvec(index1)%oper(1), &
								& nfdensity(index1))
						else
							E_en=E_en+mat_dotproduct(lattice%lvec(index1)%oper(1), &
								& nfdensity(index1))
						endif
					endif
				endif

				call pbc_get_file_and_write(lattice,nbast,nbast, &
					& index1,3,1,'            ')!3 refers to nuc-el
				call mat_free(lattice%lvec(index1)%oper(1))
				call mat_zero(H)
			endif
		endif ! is_redundant
	enddo

	call mat_free(h)
	call mat_free(h1)

	call set_lstime_print(.true.)
	write(lupri,*) 'finished with pbc_nucattrc_k'

END SUBROUTINE pbc_nucattrc_k

!> \author Johannes Rekkedal
!> \date 2013
!> \brief For computing the electron repulsion part of the fock matrix 
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param natoms 			Number of atoms in the refcell.
!> \param nbast 			Number of basis func.
!> \param lattice 		Information about the lattice.
!> \param latt_cell 		Molecule info. Threat lattice cell as molecule.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix
!> \param g2 		 		Fock matrix. e-e rep part added.
!> \param e_J 				El. repulsion (optional, calculated if present)
SUBROUTINE pbc_electron_rep_k(lupri,luerr,setting,natoms,nbast, &
		& lattice,latt_cell,refcell,numvecs,nfdensity,g_2,E_J)
	IMPLICIT NONE
	! input and output
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms 
	TYPE(matrix),intent(in),DIMENSION(numvecs) :: nfdensity
	TYPE(matrix),target,intent(inout) :: g_2(numvecs)
	TYPE(lvec_list_t),intent(inout) ::lattice
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(moleculeinfo),INTENT(IN),DIMENSION(numvecs) :: latt_cell
	TYPE(lssetting),INTENT(INOUT) :: setting 
	TYPE(MATRIX),pointer :: F_tmp(:)
	REAL(realk),INTENT(INOUT),OPTIONAL :: E_J
	! local variables
	REAL(realk) :: valmax
	INTEGER :: i,j,il21,il22,il23,iunit
	INTEGER :: il31,il32,il33, newcell,il1,il2,il3
	INTEGER :: index1,index2,index3,maxl1,maxl2,maxl3
	INTEGER :: l1,l2,l3,indred
	INTEGER(short) :: gab1,gab2,gabmaxsum
	INTEGER,SAVE :: iter=0
	INTEGER(short) :: valm1
	LOGICAL :: calc_e_j

	write(lupri,*) 'Starting pbc_coul_k'
	iter=iter+1
	call set_lstime_print(.false.)

	call mem_alloc(F_tmp,1)
	call mat_init(F_tmp(1),nbast,nbast)
	call mat_zero(F_tmp(1))

	maxl1=0
	maxl2=0
	maxl3=0

	calc_e_j=.false.
	if (PRESENT(e_j)) then
		e_j=0._realk
		calc_e_j=.true.
	end if

	do index1=1,numvecs

		!if(.not. lattice%lvec(index1)%is_redundant) then
		il1=int(lattice%lvec(index1)%lat_coord(1))
		il2=int(lattice%lvec(index1)%lat_coord(2))
		il3=int(lattice%lvec(index1)%lat_coord(3))

		lattice%lvec(index1)%g2_computed=.false.
		lattice%lvec(index1)%J_computed=.false.
		gab1=lattice%lvec(index1)%maxgab

		if(gab1  .lt. -18) CYCLE
		if(lattice%lvec(index1)%is_redundant) then
			call find_latt_index(indred,-il1,-il2,-il3,lattice,lattice%max_layer)
			if(lattice%lvec(indred)%g2_computed) then
				lattice%lvec(index1)%g2_computed=.true.
				lattice%lvec(index1)%J_computed=.true.
			endif
		endif

		if(.not. lattice%lvec(index1)%is_redundant) then
			call find_latt_index(indred,-il1,-il2,-il3,lattice,lattice%max_layer)
			if(il1**2+il2**2+il3**2 .gt. 0) lattice%lvec(indred)%is_redundant =.true.

			do index2=1,numvecs

				il21=int(lattice%lvec(index2)%lat_coord(1))
				il22=int(lattice%lvec(index2)%lat_coord(2))
				il23=int(lattice%lvec(index2)%lat_coord(3))
				gab2=lattice%lvec(index2)%maxgab

				if(gab2  .lt. -18) CYCLE

				!if(abs(il21) .gt. lattice%ndmat) CYCLE !then 
				!if( abs(il22) .gt.lattice%ndmat) CYCLE
				!if(abs(il23) .gt. lattice%ndmat) CYCLE !then 
				if(.not. lattice%lvec(index2)%dm_computed) CYCLE! then

				do index3=1,numvecs

					il31=int(lattice%lvec(index3)%lat_coord(1))
					il32=int(lattice%lvec(index3)%lat_coord(2))
					il33=int(lattice%lvec(index3)%lat_coord(3))

					if( (abs(il31) .le. lattice%nf) &
						& .and. (abs(il32) .le. lattice%nf) &
						& .and. (abs(il33) .le. lattice%nf) ) then

						!if(abs(il31) .gt. lattice%ndmat) CYCLE !then 
						!if( abs(il32) .gt.lattice%ndmat) CYCLE
						!if(abs(il33) .gt. lattice%ndmat) CYCLE !then 
						!if(nfdensity(index3)%init_magic_tag.NE.mat_init_magic_value) CYCLE
						!     ToDo Remove: base on CS screening at some point

						l1=il21+il31
						if( abs(l1) .gt.lattice%max_layer) CYCLE
						l2=il22+il32
						if( abs(l2) .gt.lattice%max_layer) CYCLE
						l3=il23+il33
						if( abs(l3) .gt.lattice%max_layer) CYCLE

						call find_latt_index(newcell,l1,l2,l3,lattice, &
							& lattice%max_layer)

						! ToDo Should be turned on at some point
						!nf1=lattice%nf
						!nf2=lattice%nf
						!nf3=lattice%nf
						!if(il1 .lt. 0)nf1=-lattice%nf
						!if(il2 .lt. 0)nf2=-lattice%nf
						!if(il3 .lt. 0)nf3=-lattice%nf
						!if(il21+il1 .gt. lattice%nf .and. il1 .gt. 0 ) cycle
						!if(il21 .lt. nf1 .and. il1 .lt. 0) CYCLE
						!if(il21 .lt. nf1

						gabmaxsum=gab1+gab2
						call mat_abs_max_elm(nfdensity(index2),valmax)
						if(valmax .gt. 0._realk) valm1=int(log10(valmax),kind=short)

						if(valm1+gabmaxsum  .ge. -12) then
							!if(gabmaxsum  .ge. -12) then
							lattice%lvec(index1)%g2_computed=.true.
							lattice%lvec(index1)%J_computed=.true.
							if(lattice%lvec(index1)%oper(2)%init_magic_tag &
								& .ne. mat_init_magic_value) then
								call mat_init(lattice%lvec(index1)%oper(2),nbast,nbast)
								call mat_zero(lattice%lvec(index1)%oper(2))
							endif

							maxl1=max(il1,maxl1)
							maxl2=max(il2,maxl2)
							maxl3=max(il3,maxl3)

							call TYPEDEF_setmolecules(setting,refcell,1, &
								& latt_cell(index1), 2,latt_cell(index3),3, &
								& latt_cell(newcell),4)
							call mat_zero(F_tmp(1))
							call II_get_coulomb_mat(lupri,luerr,setting, &
								& nfdensity(index2:index2),F_tmp,1)

							call mat_daxpy(0.5_realk,F_tmp(1), &
								& lattice%lvec(index1)%oper(2))
						endif
					endif
				enddo
			enddo

			if(lattice%lvec(index1)%J_computed) then
				if(.not. lattice%store_mats) then
					call mat_init(g_2(index1),nbast,nbast)
					call mat_copy(1._realk,lattice%lvec(index1)%oper(2),g_2(index1))
				endif

				if (calc_e_j) then
					if((abs(il1).le.lattice%ndmat .and. abs(il2).le.lattice%ndmat) &
						& .and. abs(il3).le.lattice%ndmat) then
						if(nfdensity(index1)%init_magic_tag.EQ.mat_init_magic_value) THEN
							if(il1**2+il2**2+il3**2 .gt. 0) then
								E_j=E_j+ &
									& 1.0_realk*mat_dotproduct(lattice%lvec(index1)%oper(2), &
									& nfdensity(index1))
							else
								E_j=E_j+ &
									& 0.5_realk*mat_dotproduct(lattice%lvec(index1)%oper(2), &
									& nfdensity(index1))
							endif
						endif
					endif
				endif
			end if
		endif !is_redundant
		lattice%lvec(index1)%is_redundant=.false.
	enddo
	lattice%col1=maxl1
	lattice%col2=maxl2
	lattice%col3=maxl3

	call mat_free(F_tmp(1))
	call mem_dealloc(F_tmp)

	write(lupri,*) 'Check, nr. lattice vectors',numvecs
	call set_lstime_print(.true.)
	write(lupri,*) 'finished with pbc_elrep_k'

END SUBROUTINE pbc_electron_rep_k

!> \author Johannes Rekkedal
!> \date 2013
!> \brief For computing the excange part of the fock matrix 
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param natoms 			Number of atoms in the refcell.
!> \param nbast 			Number of basis func.
!> \param lattice 		Information about the lattice.
!> \param latt_cell 		Molecule info. Threat lattice cell as molecule.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix
!> \param g2 		 		Fock matrix. e-e rep part added.
!> \param e_k 				Exchange (optional, calculated if present)
SUBROUTINE pbc_exact_xc_k(lupri,luerr,setting,natoms,nbast, &
		& lattice,latt_cell,refcell,numvecs,nfdensity,g_2,E_K)
	IMPLICIT NONE
	! input and output
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(matrix),INTENT(IN) :: nfdensity(numvecs)
	TYPE(matrix),INTENT(INOUT) :: g_2(numvecs)
	TYPE(moleculeinfo),INTENT(IN),DIMENSION(numvecs) :: latt_cell
	TYPE(LSSETTING),INTENT(INOUT) :: SETTING 
	TYPE(lvec_list_t),INTENT(INOUT) ::lattice
	REAL(realk),INTENT(INOUT),OPTIONAL :: E_K
	! local variables
	TYPE(MATRIX),pointer :: K_tmp(:)
	TYPE(MATRIX) :: Kx
	INTEGER :: il21,il22,il23,gabl1,gabl2,gabl3
	INTEGER :: il31,il32,il33,newcell,il1,il2,il3
	INTEGER :: index1,index2,index3,gabind
	INTEGER :: l1,l2,l3,maxl1,maxl2,maxl3,indred
	INTEGER(short) :: gab1,gab2,maxgabsum
	INTEGER, SAVE :: iter=0
	INTEGER(short) :: valm1
	REAL(realk) :: valmax
	LOGICAL :: calc_e_k

	iter=iter+1

	calc_e_k=.false.
	if (PRESENT(e_k)) then
		calc_e_k=.true.
		e_k=0.0_realk
	end if

	call set_lstime_print(.false.)

	call mem_alloc(K_tmp,1)
	call mat_init(Kx,nbast,nbast)
	call mat_zero(Kx)
	call mat_init(K_tmp(1),nbast,nbast)
	call mat_zero(K_tmp(1))

	maxl1=0
	maxl2=0
	maxl3=0

	do index1=1,numvecs

		!call find_latt_vectors(index1,il1,il2,il3,lattice)
		il1=int(lattice%lvec(index1)%lat_coord(1))
		il2=int(lattice%lvec(index1)%lat_coord(2))
		il3=int(lattice%lvec(index1)%lat_coord(3))
		!if(abs(il1) .gt. lattice%kx1) CYCLE
		!if(abs(il2) .gt. lattice%kx2) CYCLE
		!if(abs(il3) .gt. lattice%kx3) CYCLE
		lattice%lvec(index1)%Kx_computed=.false.

		if(lattice%lvec(index1)%is_redundant) then
			call find_latt_index(indred,-il1,-il2,-il3,lattice,lattice%max_layer)
			if(lattice%lvec(indred)%Kx_computed) then
				lattice%lvec(index1)%g2_computed=.true.
				lattice%lvec(index1)%Kx_computed=.true.
			endif
		endif

		if(.not. lattice%lvec(index1)%is_redundant) then
			call find_latt_index(indred,-il1,-il2,-il3,lattice,lattice%max_layer)

			do index2=1,numvecs

				!call find_latt_vectors(index2,il21,il22,il23,lattice)
				il21=int(lattice%lvec(index2)%lat_coord(1))
				il22=int(lattice%lvec(index2)%lat_coord(2))
				il23=int(lattice%lvec(index2)%lat_coord(3))
				gab1=lattice%lvec(index2)%maxgab
				if(gab1 .lt. -18) CYCLE

				do index3=1,numvecs

					!call find_latt_vectors(index3,il31,il32,il33,lattice)
					il31=int(lattice%lvec(index3)%lat_coord(1))
					il32=int(lattice%lvec(index3)%lat_coord(2))
					il33=int(lattice%lvec(index3)%lat_coord(3))

					!if(abs(il31) .gt. lattice%ndmat) CYCLE! then
					!if(abs(il32) .gt. lattice%ndmat) CYCLE !then 
					!if(abs(il33) .gt. lattice%ndmat) CYCLE !then 
					if(.not. lattice%lvec(index3)%dm_computed) CYCLE! then

					!ToDo Remove: base on CS screening at some point
					l1=il21+il31
					if( abs(l1) .gt.lattice%max_layer) CYCLE
					!todo l1 blir større enn max slik at newcell blir negativ
					!sjekk dette, er ikke sikker på om l1 skal være slik.
					!Likedan man l2 og l3
					l2=il22+il32
					if( abs(l2) .gt.lattice%max_layer) CYCLE
					l3=il23+il33
					if( abs(l3) .gt.lattice%max_layer) CYCLE
					gabl1=l1-il1
					gabl2=l2-il2
					gabl3=l3-il3
					if( abs(gabl1) .gt.lattice%max_layer) CYCLE
					if( abs(gabl2) .gt.lattice%max_layer) CYCLE
					if( abs(gabl3) .gt.lattice%max_layer) CYCLE

					call find_latt_index(newcell,l1,l2,l3,lattice, &
						& lattice%max_layer)
					call find_latt_index(gabind,gabl1,gabl2,gabl3,lattice, &
						& lattice%max_layer)

					gab2=lattice%lvec(gabind)%maxGab
					if(gab2 .lt. -18) CYCLE
					maxgabsum=gab1+gab2
					call mat_abs_max_elm(nfdensity(index3),valmax)
					if(valmax.gt.0._realk) valm1=int(log10(valmax),kind=short)

					if(maxgabsum+valm1 .ge. -10) then
						if(lattice%lvec(index1)%oper(1)%init_magic_tag &
							& .ne.mat_init_magic_value) then
							call mat_init(lattice%lvec(index1)%oper(1),nbast,nbast)
							call mat_zero(lattice%lvec(index1)%oper(1))
						endif
						lattice%lvec(index1)%g2_computed=.true.
						lattice%lvec(index1)%Kx_computed=.true.
						maxl1=max(abs(il1),maxl1)
						maxl2=max(abs(il2),maxl2)
						maxl3=max(abs(il3),maxl3)
						lattice%kx1=maxl1
						lattice%kx2=maxl2
						lattice%kx3=maxl3

						call typedef_setmolecules(setting,refcell,1, &
							& latt_cell(index1),3,latt_cell(index2),2, &
							& latt_cell(newcell),4)
						call mat_zero(K_tmp(1))
						call II_get_exchange_mat(lupri,luerr,setting, &
							& nfdensity(index3:index3),1,.false.,K_tmp)
						call mat_daxpy(1._realk,K_tmp(1),Kx)
						call mat_daxpy(0.5_realk,K_tmp(1), &
							& lattice%lvec(index1)%oper(1))

						if(lattice%lvec(indred)%oper(1)%init_magic_tag &
							& .ne.mat_init_magic_value) then
							call mat_init(lattice%lvec(indred)%oper(1),nbast,nbast)
							call mat_zero(lattice%lvec(indred)%oper(1))
						endif
					endif
				enddo
			enddo

			if((il1**2 + il2**2 +il3**2 .gt. 0) &
				& .and. lattice%lvec(index1)%Kx_computed) then
				call mat_trans(lattice%lvec(index1)%oper(1), &
					& lattice%lvec(indred)%oper(1))
				lattice%lvec(indred)%is_redundant=.true.
			endif
			if(lattice%lvec(index1)%Kx_computed) then
				write(lupri,*) 'Kx computed for', il1,il2,il3
			endif
		endif !is_redundant

		if(lattice%lvec(index1)%Kx_computed) then

			if(lattice%store_mats) then
				call pbc_get_file_and_write(lattice,nbast,nbast,index1,5,1, &
					& '            ')!5 refers to
			else
				if(g_2(index1)%init_magic_tag.ne.mat_init_magic_value) then
					call mat_init(g_2(index1),nbast,nbast)
					call mat_zero(g_2(index1))
				endif
				call mat_daxpy(1._realk,lattice%lvec(index1)%oper(1),g_2(index1))
			endif

			if (calc_e_k) then 
				if((abs(il1) .le. lattice%ndmat) &
					& .and. (abs(il2) .le. lattice%ndmat) &
					& .and. (abs(il3) .le. lattice%ndmat)) then
				if(nfdensity(index1)%init_magic_tag.EQ.mat_init_magic_value) then
					E_K=E_K+0.5_realk*mat_dotproduct(lattice%lvec(index1)%oper(1), &
						& nfdensity(index1))
				endif
			endif
		endif
	endif

	if(lattice%lvec(index1)%oper(1)%init_magic_tag.eq.mat_init_magic_value) then
		call mat_free(lattice%lvec(index1)%oper(1))
	endif
	call mat_zero(Kx)
enddo
lattice%Kx1=maxl1
lattice%Kx2=maxl2
lattice%Kx3=maxl3

call mat_free(Kx)
call mat_free(K_tmp(1))
call mem_dealloc(K_tmp)

write(lupri,*) 'Check, nr. lattice vectors',numvecs
call set_lstime_print(.true.)
write(lupri,*) 'finished with pbc_exact_xc_k'

END SUBROUTINE pbc_exact_xc_k


!> \author Johannes Rekkedal
!> \date 2013
!> \brief For computing overlap integrals 
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param nbast 			Number of basis func.
!> \param lattice 		Information about the lattice.
!> \param latt_cell 		Molecule info. Threat lattice cell as molecule.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
SUBROUTINE pbc_overlap_int(lupri,luerr,setting,nbast,lattice,latt_cell,refcell,numvecs)
  IMPLICIT NONE
  ! input and output variables
  INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs
  TYPE(lssetting),INTENT(INOUT) :: SETTING 
  TYPE(moleculeinfo),INTENT(INOUT) :: refcell
  TYPE(moleculeinfo),INTENT(IN) :: latt_cell(numvecs)
  TYPE(lvec_list_t),INTENT(INOUT) :: lattice
  ! local variables
  TYPE(matrix)  :: S
  INTEGER :: il1,il2,il3
  INTEGER ::  indx

  call mat_init(S,nbast,nbast)
  call mat_zero(S)
  write(lupri,*) 'Number of lattice vectors ', numvecs

  do indx=1,numvecs
	  call typedef_setmolecules(setting,refcell,1,latt_cell(indx),3)
     setting%samemol(1,3)=.false.
     setting%samemol(3,1)=.false.

     call find_latt_vectors(indx,il1,il2,il3,lattice)
     if(abs(il1) .gt. lattice%nneighbour) CYCLE
     if(abs(il2) .gt. lattice%nneighbour) CYCLE
     if(abs(il3) .gt. lattice%nneighbour) CYCLE
     call II_get_overlap(lupri,luerr,setting,S)
  end do

  call mat_free(S)
  write(lupri,*) 'finished pbc_overlap_int'

END SUBROUTINE pbc_overlap_int

!
SUBROUTINE pbc_kinetic_int(lupri,luerr,setting,molecule,nbast,fock_mtx,sizef,lattice,latt_cell,refcell,numvecs)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri, luerr,nbast,sizef,numvecs
  COMPLEX(complexk), INTENT(INOUT) :: fock_mtx(sizef*sizef)
  TYPE(moleculeinfo),INTENT(IN) :: molecule
  TYPE(LSSETTING),intent(inout)   :: SETTING 

  ! local variables
  REAL(realk), DIMENSION(12) :: Tvec
  TYPE(moleculeinfo)  :: refcell
  TYPE(moleculeinfo), DIMENSION(numvecs),INTENT(IN) :: latt_cell
  TYPE(MATRIX)  :: kin
  REAL(realk) :: std_vec_length
  INTEGER :: i,j,il1,il2,il3
  INTEGER ::  index, refindex
  INTEGER :: num_latvectors, natoms
  REAL(realk) :: latt_vec_std(3),origin(3)
  REAL(realk) :: fock_tmp(nbast*nbast)
  TYPE(lvec_list_t),intent(INOUT) ::lattice

  write(lupri,*) 'Starting pbc_kinetic_int'
  natoms=molecule%natoms


  origin(1:3)=1.0_realk
  latt_vec_std(1:3)=0.0_realk

  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)

  call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)
  write(lupri,*) 'reference cell index=', refindex,lattice%max_layer
  fock_tmp(:)= 0.0_realk
  tvec(1:12)=0.0_realk
  i=1

  call mat_init(kin,nbast,nbast)
  call mat_zero(kin)
  write(lupri,*) 'nbast',nbast,kin%ncol,kin%nrow

  DO index=1,num_latvectors
     


     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index),3)

     !call calc_distance(distance,lattice%lvec(index)%std_coord,latt_vec_std)
     call find_latt_vectors(index,il1,il2,il3,lattice)
     if(abs(il1) .gt. lattice%nneighbour) CYCLE
     if(abs(il2) .gt. lattice%nneighbour) CYCLE
     if(abs(il3) .gt. lattice%nneighbour) CYCLE
     !phase=cmplx(0,k1*il1+k2*il2+k3*il3)
    

     setting%samemol(1,3)=.false.
     setting%samemol(3,1)=.false.
     call II_get_kinetic(lupri,luerr,setting,kin)

!     if(lattice%compare_elmnts) then
!       !compare integrals with the old pbc code
!       write(lupri,*) 'comparing kinetic elements with old pbc code'
!       matris(:,:) =0.0
!       write(lupri,*) il1,il2,il3
!       call pbc_readopmat2(il1,il2,il3,matris,nbast,'KINETIC',.true.,.false.)
!       call write_matrix(matris,nbast,nbast)
!!       call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
!       write(lupri,*) ''
!       write(lupri,*) 'To compare with the matrix below'
!       CALL mat_print(kin1,1,kin1%nrow,1,kin1%ncol,6)
!     endif

     i=1
     j=0
!     call mat_print(kin1,1,kin1%nrow,1,kin1%ncol,6)
!     write(lupri,*) kin1%elms
!     STOP
  !   DO k=1,nbast*nbast
  !      j=j+1
  !      if(j .gt. nbast) Then
  !        j=1
  !        i=i+1
  !      ENDIF
  !      fock_tmp(j,i)=kin1%elms(k)
    ! ENDDO
    fock_tmp=kin%elms!*exp(phase)
     write(lupri,*) fock_tmp
     i=0
   !  DO k=refindex*nbast-nbast+1,refindex*nbast 
   !    i=i+1
   !    j=
   !  DO m=index*nbast-nbast+1,index*nbast 
   !    j=j+1
   !    fock_mtx(k,m)=fock_tmp(i,j)!*coeff(k,m)
   !    fock_mtx(m,k)=fock_mtx(k,m)!dagger if complex
   !  ENDDO
   !  ENDDO
        fock_mtx=fock_tmp

  END DO

  call mat_free(kin)
  write(lupri,*) 'Finished pbc_kinetic_int'

!!
END SUBROUTINE pbc_kinetic_int
!
!
SUBROUTINE pbc_nucattrc_int(lupri,luerr,setting,molecule,nbast,fock_mtx,sizef,lattice,latt_cell,refcell,numvecs)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri, luerr,nbast,sizef,numvecs
  TYPE(moleculeinfo),INTENT(IN) :: molecule,refcell
  COMPLEX(complexk), INTENT(INOUT) :: fock_mtx(sizef*sizef)
  TYPE(LSSETTING),intent(inout)   :: SETTING 
!  TYPE(DALTONINPUT) :: INPUT
  TYPE(MATRIX)  :: H,H1

  ! local variables
  REAL(realk), DIMENSION(12) :: Tvec
!  TYPE(moleculeinfo), pointer :: lattice_cell(:)
  TYPE(moleculeinfo), DIMENSION(numvecs),INTENT(IN) :: latt_cell
  REAL(realk) :: std_vec_length
  INTEGER :: il11,il12,il13,il21,il22,il23
  INTEGER ::  index1,index2, refindex
  INTEGER ::  num_latvectors, natoms
  REAL(realk)::  latt_vec_std(3),origin(3)
  REAL(realk) :: fock_tmp(nbast*nbast)
  TYPE(lvec_list_t),intent(INOUT) ::lattice

  natoms=molecule%natoms


  call mat_init(H,nbast,nbast)
  call mat_zero(H)
  call mat_init(H1,nbast,nbast)
  call mat_zero(H1)

  origin(1:3)=1.0
  latt_vec_std(:)=0.0_realk


  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)
  
  call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)
  write(lupri,*) 'reference cell index=', refindex

  tvec(1:12)=0.0_realk
 

  DO index1=1,num_latvectors



!  call calc_distance(distance,lattice%lvec(index1)%std_coord,latt_vec_std)
  call find_latt_vectors(index1,il11,il12,il13,lattice)
  if(abs(il11) .gt. lattice%nneighbour) CYCLE
  if(abs(il12) .gt. lattice%nneighbour) CYCLE
  if(abs(il13) .gt. lattice%nneighbour) CYCLE
  !phase=(0,k1*il1+k2*il2+k3*il3)


  DO index2=1,num_latvectors


     call find_latt_vectors(index2,il21,il22,il23,lattice)
     if(abs(il21) .gt. lattice%nneighbour) CYCLE !uncomment this
     if(abs(il22) .gt. lattice%nneighbour) CYCLE
     if(abs(il23) .gt. lattice%nneighbour) CYCLE

    ! call calc_distance(distance,lattice%lvec(index2)%std_coord,latt_vec_std)

     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index1),2,latt_cell(index2),3)

         setting%samemol(1,3)=.false.
         setting%samemol(3,1)=.false.
         setting%samemol(1,2)=.false.
         setting%samemol(2,1)=.false.
         setting%samemol(3,2)=.false.
         setting%samemol(2,3)=.false.

         call II_get_nucel_mat(lupri,luerr,setting,H1)
         call mat_add(1E0_realk,H,1E0_realk,H1,H)
         
!         write(lupri,*) 'H%elms', H%elms


  ENDDO

!     i=1
!     j=0
!     DO k=1,nbast*nbast
!        j=j+1
!        if(j .gt. nbast) Then
!          j=1
!          i=i+1
!        ENDIF
!        fock_tmp(j,i)=H%elms(k)
!     ENDDO
     !call mat_print(H,1,H%nrow,1,H%ncol,6)
     !call write_matrix(fock_tmp,2,2)
     !STOP
     fock_tmp=H1%elms
!     i=0
!     DO k=refindex*nbast-nbast+1,refindex*nbast 
!       i=i+1
!       j=0
!     DO m=index1*nbast-nbast+1,index1*nbast 
!       j=j+1
!       fock_mtx(k,m)=fock_mtx(k,m)+fock_tmp(i,j)!*coeff(k,m)
!       fock_mtx(m,k)=fock_mtx(k,m)!dagger if complex
!     ENDDO
!     ENDDO

     fock_mtx=fock_mtx+fock_tmp!*coeff(k,m)
     call mat_print(H1,1,H1%nrow,1,H1%ncol,6)
!    write(lupri,*) 'H%elms', H%elms
!  if(lattice%compare_elmnts) then
!    !comparing integrals with the old pbc code
!    write(lupri,*) 'comparing nucleus-electron elements with old pbc code'
!    matris(:,:) =0.0
!
!    write(lupri,*) il11,il12,il13
!
!    call pbc_readopmat2(il11,il12,il13,matris,nbast,'NUCVNF2',.true.,.false.)
!    call write_matrix(matris,nbast,nbast)
!    write(lupri,*) ''
!    write(lupri,*) 'To compare with the matrix below'
!    CALL mat_print(H,1,H%nrow,1,H%ncol,6)
!  endif

    call mat_zero(H1)

  ENDDO

  call mat_free(H)
  call mat_free(H1)
  write(lupri,*) 'finished with pbc_nucattrc_int'

END SUBROUTINE pbc_nucattrc_int
!
!
!SUBROUTINE for computing the nuclear repulsions
SUBROUTINE pbc_nucpot(lupri,luerr,setting,molecule,lattice,&
 latt_cell,refcell,numvecs,E_nn)

  implicit none
  ! input and output arguments
  INTEGER, INTENT(IN) :: lupri, luerr, numvecs
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule
  TYPE(moleculeinfo),intent(inout) :: refcell
  TYPE(lvec_list_t),intent(INOUT) ::lattice
  Real(realk),INTENT(INOUT) :: E_nn

  ! local variables
  REAL(realk) :: kvec(3)
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(moleculeinfo),intent(in) :: latt_cell(numvecs)
  INTEGER :: i,j,m,il1,il2,il3
  INTEGER :: index,refindex
  INTEGER :: natoms,iunit,maxl1,maxl2,maxl3
  REAL(realk) :: latt_vec_std(3)
  REAL(realk) :: nucpot



  natoms=molecule%natoms!number of atoms

  call set_lstime_print(.false.)

  latt_vec_std(:)=0.0_realk

  write(lupri,*) 'Number of lattice vectors ', numvecs

  !finds the lattice index for the reference cell.
  !i.e for l1=l2=l3=0
  call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)

  E_nn=0._realk
  !loop over lattice cells

  DO index=1,numvecs
     

  !Doing translations
     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index),3)
     setting%samemol(1,3)=.false.
     setting%samemol(3,1)=.false.

     call find_latt_vectors(index,il1,il2,il3,lattice)
     !So that we do not consider negligible integrals
     if(abs(il1) .gt. lattice%nf) CYCLE
     if(abs(il2) .gt. lattice%nf) CYCLE
     if(abs(il3) .gt. lattice%nf) CYCLE

     if(index .eq. refindex) then

       call II_get_nucpot(lupri,luerr,setting,nucpot)
     else
       call pbc_get_nucpot(lupri,luerr,setting,nucpot)
     endif
     E_nn=E_nn+nucpot

  END DO
  !E_nn=E_nn/2.

 ! write(*,*) 'Debug 2', E_nn

END SUBROUTINE pbc_nucpot

!> \brief Calculates the nuclear repulsion energy contribution
!> \author S. Reine and T. Kjaergaard
!> \modified by J. Rekkedal
!> \date 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nucpot the nuclear repulsion energy contribution
SUBROUTINE pbc_get_nucpot(LUPRI,LUERR,SETTING,NUCPOT)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
integer             :: usemat
INTEGER               :: LUPRI,LUERR
REAL(realk)           :: nucpot
Integer               :: I,J
real(realk)           :: pq(3),distance
logical               :: NOBQBQ
call time_II_operations1()
NOBQBQ = SETTING%SCHEME%NOBQBQ
NUCPOT=0.0E0_realk
DO I=1,SETTING%MOLECULE(1)%p%Natoms
 IF(SETTING%MOLECULE(1)%p%ATOM(I)%phantom)CYCLE
 DO J=1,SETTING%MOLECULE(1)%p%Natoms
  IF(SETTING%MOLECULE(1)%p%ATOM(J)%phantom)CYCLE
  
  if(setting%molecule(1)%p%ATOM(I)%Pointcharge .and. &
    &setting%molecule(1)%p%ATOM(J)%Pointcharge .and. NOBQBQ) cycle

  pq(1) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1)-SETTING%MOLECULE(3)%p%ATOM(J)%CENTER(1)
  pq(2) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2)-SETTING%MOLECULE(3)%p%ATOM(J)%CENTER(2)
  pq(3) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3)-SETTING%MOLECULE(3)%p%ATOM(J)%CENTER(3)
  Distance = sqrt(pq(1)*pq(1)+pq(2)*pq(2)+pq(3)*pq(3))
  NUCPOT=NUCPOT+SETTING%MOLECULE(1)%p%ATOM(I)%Charge*SETTING%MOLECULE(3)%p%ATOM(J)%Charge/Distance
 ENDDO
ENDDO
NUCPOT=NUCPOT/2.
call time_II_operations2(JOB_II_get_nucpot)
END SUBROUTINE pbc_get_nucpot

SUBROUTINE pbc_electron_rep(lupri,luerr,setting,molecule,&
 nbast,lattice,latt_cell,refcell,numvecs,nfdensity,nfsze)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri, luerr ,nbast,numvecs,nfsze
  TYPE(moleculeinfo),INTENT(IN) :: molecule,refcell
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(matrix),intent(inout),DIMENSION(nfsze) :: nfdensity
  TYPE(MATRIX)  :: F(1),F_tmp(1)

  ! local variables
  REAL(realk), DIMENSION(12) :: Tvec
  !TYPE(moleculeinfo), pointer :: lattice_cell(:)
  TYPE(moleculeinfo), INTENT(IN), DIMENSION(numvecs) :: latt_cell
  !TYPE(lattice_cell_info_t),pointer :: lat_cells(:)
  REAL(realk) :: std_vec_length
  INTEGER :: il21,il22,il23
  INTEGER :: il31,il32,il33, newcell,il1,il2,il3
  INTEGER ::  index1,index2,index3
  INTEGER :: l1,l2,l3, num_latvectors, natoms
  INTEGER :: checknf,checknf1,checknf2,checknf3
  REAL(realk):: latt_vec_std(3),origin(3)
  TYPE(lvec_list_t),intent(inout) ::lattice

  natoms=molecule%natoms

!  call mat_init(D(1),nbast,nbast)
!  call mat_zero(D(1))
  call mat_init(F(1),nbast,nbast)
  call mat_zero(F(1))
  call mat_init(F_tmp(1),nbast,nbast)
  call mat_zero(F_tmp(1))

  origin(1:3)=1.0_realk
  latt_vec_std(:)=0.0_realk
 

  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)


  tvec(1:12)=0.0_realk
 

  DO index1=1,num_latvectors

  !call calc_distance(distance,lattice%lvec(index1)%std_coord,latt_vec_std)
  call find_latt_vectors(index1,il1,il2,il3,lattice)
!  if(abs(il1) .gt. lattice%nneighbour) CYCLE
!  if(abs(il2) .gt. lattice%nneighbour) CYCLE
!  if(abs(il3) .gt. lattice%nneighbour) CYCLE

  DO index2=1,num_latvectors

     call find_latt_vectors(index2,il21,il22,il23,lattice)
  if(abs(il21) .gt. lattice%nneighbour) CYCLE
  if(abs(il22) .gt. lattice%nneighbour) CYCLE
  if(abs(il23) .gt. lattice%nneighbour) CYCLE

     !Skal kun ha avstand til origo ikke latt_std_vec her
!     call calc_distance(distance,lattice%lvec(index2)%std_coord,latt_vec_std)


   DO index3=1,num_latvectors

    !write(lupri,*) 'electron rep get'
    !CALL mat_print(F(1),1,F(1)%nrow,1,F(1)%ncol,6)
    !CALL mat_print(D1(1),1,F(1)%nrow,1,F(1)%ncol,6)

    call find_latt_vectors(index3,il31,il32,il33,lattice)
!  if(abs(il31) .gt. lattice%nneighbour) CYCLE
!  if(abs(il32) .gt. lattice%nneighbour) CYCLE
!  if(abs(il33) .gt. lattice%nneighbour) CYCLE

    l1=il21+il31
    if(abs(l1) .gt. lattice%max_layer) CYCLE! then
    
    l2=il22+il32
    if(abs(l2) .gt. lattice%max_layer) CYCLE !then 
    
    l3=il23+il33
    if(abs(l3) .gt. lattice%max_layer) CYCLE !then 
    
    checknf1=l1-il21
    if(abs(checknf1) .gt. n_neighbour) CYCLE! then
    
    checknf2=l2-il22
    if(abs(checknf2) .gt. n_neighbour) CYCLE !then 
    
    checknf3=l3-il23
    if(abs(checknf3) .gt. n_neighbour) CYCLE !then 

    
    call find_latt_index(newcell,l1,l2,l3,lattice,lattice%max_layer)

    call find_latt_index(checknf,checknf1,checknf2,checknf3,lattice,&
     lattice%nneighbour)

    !Changes setting to point at different lattice cells
    call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index1),2,&
     latt_cell(index2),3,latt_cell(newcell),4)

    setting%samemol(1,2)=.false.
    setting%samemol(2,1)=.false.
    setting%samemol(1,3)=.false.
    setting%samemol(3,1)=.false.
    setting%samemol(1,4)=.false.
    setting%samemol(4,1)=.false.
    setting%samemol(2,3)=.false.
    setting%samemol(2,3)=.false.
    setting%samemol(3,4)=.false.
    setting%samemol(4,2)=.false.
    setting%samemol(3,4)=.false.
    setting%samemol(4,3)=.false.

 !   setting%samefrag=.false.
    
    call II_get_coulomb_mat(lupri,luerr,setting,nfdensity(checknf:checknf),F_tmp,1)
    !call II_get_exchange_mat(lupri,luerr,setting,D,F_tmp,1)
    
    call mat_add(1E0_realk,F(1),1E0_realk,F_tmp(1),F(1))
    
    !call exchange()
!    write(lupri,*) 'electron rep fin rep',nbast,F(1)%nrow,F(1)%ncol
!    CALL mat_print(F(1),1,F_tmp(1)%nrow,1,F_tmp(1)%ncol,6)
!     STOP
   ENDDO
  ENDDO
  !store the matrix here

!    if(lattice%compare_elmnts) then
!      
!      !compare integrals with the old pbc code
!      write(lupri,*) 'comparing electron-electron elements with old pbc code'
!      matris(:,:) =0.0
!
!      write(lupri,*) il1,il2,il3
!
!      call pbc_readopmat2(il1,il2,il3,matris,nbast,'COULOMB',.true.,.false.)
!      call write_matrix(matris,nbast,nbast)
!
!      write(lupri,*) ''
!      write(lupri,*) 'To compare with the matrix below'
!      CALL mat_print(F(1),1,F(1)%nrow,1,F(1)%ncol,6)
!
!    endif

    CALL mat_print(F(1),1,F_tmp(1)%nrow,1,F_tmp(1)%ncol,6)
    call mat_zero(F(1))

    write(lupri,*) 'index1',index1, num_latvectors
  ENDDO
  call mat_free(F(1))
  call mat_free(F_tmp(1))

  write(lupri,*) 'finished with pbc_elrep_int'
  call LSQuit('finished with pbc_elrep_int',lupri)

END SUBROUTINE pbc_electron_rep


SUBROUTINE pbc_complete_Fock_mtx(lupri,nbast,fock_mtx,sizef,cut,lattice)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri,nbast,sizef
  REAL(realk), INTENT(INOUT) :: fock_mtx(sizef,sizef)
  REAL(realk), INTENT(IN) :: cut
  ! Local variables
  REAL(realk) :: distance
  TYPE(lvec_list_t), intent(in) ::lattice
  INTEGER :: num_latvectors
  INTEGER :: cellij,i,j
  INTEGER :: diff1,diff2,diff3
  INTEGER :: il1,il2,il3,jl1,jl2,jl3
  REAL(realk) ::  latt_vec_std(3), origin(3)
  INTEGER :: refcell

!  call build_lvec_list(lattice)
  !call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)
  
!  write(lupri,*) 'debug: ',lattice%lvec(2)%std_coord!, latstdvec 

  num_latvectors=size(lattice%lvec)

  origin(1:3)=1.0_realk
  latt_vec_std(:)=0.0_realk


  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)
  !std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  call find_latt_index(refcell,0,0,0,lattice,lattice%max_layer)
  write(lupri,*) 'refcelle', refcell

  DO i=1,num_latvectors
   IF(i .eq. refcell) CYCLE
   fock_mtx(nbast*i+1-nbast:i*nbast,nbast*i+1-nbast:i*nbast)=&
   fock_mtx(refcell*nbast-nbast+1:refcell*nbast,refcell*nbast-nbast+1:refcell*nbast)
   call find_latt_vectors(i,il1,il2,il3,lattice)
   DO j=1,i
      IF(j .eq. refcell) CYCLE
      IF(j .eq. i) CYCLE
      call find_latt_vectors(j,jl1,jl2,jl3,lattice)
      !write(lupri,*) 'SEGMENTATION FAULT',j
!      write(lupri,*) 'refcelle', refcell,j
      diff1=il1-jl1
      if(abs(diff1) .gt. lattice%max_layer) CYCLE !should be max_layer
      if(abs(diff1) .gt. lattice%nneighbour) CYCLE !should be threshold
      diff2=il2-jl2
      if(abs(diff2) .gt.lattice%max_layer) CYCLE !should be max_layer
      if(abs(diff2) .gt.lattice%nneighbour) CYCLE !should be threshold
      diff3=il3-jl3
      if(abs(diff3) .gt. lattice%max_layer) CYCLE !should be max_layer
      if(abs(diff3) .gt. lattice%nneighbour) CYCLE !should be threshold

      call find_latt_index(cellij,diff1,diff2,diff3,lattice,lattice%max_layer)

      call calc_distance(distance,lattice%lvec(cellij)%std_coord,latt_vec_std)

      IF(distance .ge. cut) CYCLE
      fock_mtx(nbast*i+1-nbast:i*nbast,nbast*j+1-nbast:j*nbast)=&
      fock_mtx(refcell*nbast-nbast+1:refcell*nbast,cellij*nbast-nbast+1:cellij*nbast)
      fock_mtx(nbast*j+1-nbast:j*nbast,nbast*i+1-nbast:i*nbast)=&
      fock_mtx(nbast*i+1-nbast:i*nbast,nbast*j+1-nbast:j*nbast)
 !     IF(i /= j) THEN
 !       diff= i-j
 !       fock_mtx(nbast*i+1-nbast:i*nbast,nbast*j+1-nbast:j*nbast)=&
 !       fock_mtx(refcell*nbast-nbast+1:refcell*nbast,diff*nbast-nbast+1:diff*nbast)
 !       fock_mtx(nbast*j+1-nbast:j*nbast,nbast*i+1-nbast:i*nbast)=&
 !       fock_mtx(nbast*i+1-nbast:i*nbast,nbast*j+1-nbast:j*nbast)
 !     ENDIF
 !   ! DO k=refindex*nbast-nbast+1,refindex*nbast 
 !   !   i=i+1
 !   !   j=0
 !   !  DO m=index1*nbast-nbast+1,index1*nbast 
 !   !    fock_mtx(k,m)=fock_mtx(k,m)+fock_tmp(i,j)!*coeff(k,m)
 !   !    fock_mtx(m,k)=fock_mtx(k,m)
 !   !  ENDDO
 !   ! ENDDO
   ENDDO
  ENDDO
  write(lupri,*) 'FOCK MATRIX'
  call write_matrix(fock_mtx,sizef,sizef)

END SUBROUTINE pbc_complete_Fock_mtx

END MODULE pbc_interactions

#endif

!!todo not in use
!SUBROUTINE screening_ovl(basinfo)
!IMPLICIT NONE
!TYPE(basissetinfo) :: basinfo
!INTEGER :: natomtypes
!INTEGER :: i,j,k
!!REAL(realk) :: PI=3.14159265
!REAL(realk), pointer :: minexp(:)
!REAL(realk) :: distance
!
!  natomtypes=basinfo%natomtypes
!  
!  k=0
!  DO i=1,natomtypes
!  DO j=1,basinfo%atomtype(i)%nangmom
!  k=k+1
!  ENDDO
!  ENDDO
!  call mem_alloc(minexp,k)
!  k=0
!  DO i=1,natomtypes
!  DO j=1,basinfo%atomtype(i)%nangmom
!  k=k+1
!  minexp(k)=minval(basinfo%atomtype(i)%shell(j)%segment(2)%exponents)
!  ENDDO
!  ENDDO
!  
!  !write(lupri,*) 'expon ', minval(minexp)
!  distance=sqrt(1./minval(minexp)*log((pi/(2.*minval(minexp)))**3E0_realk*10E0_realk**20))
!  !write(lupri,*) 'distance ',distance
!
!!  distance=sqrt(1./0.2979640)*sqrt(log((pi/(2.*0.2979640))**(3)*10**20))
!!  write(lupri,*) 'distance ',distance
!   call mem_dealloc(minexp)
!
!
!END SUBROUTINE screening_ovl
!
