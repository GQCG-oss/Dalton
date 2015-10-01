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
!> \param refcell 		Molecule info reference cell.
SUBROUTINE find_cutoff_onep(lupri,luerr,setting,nbast,refcell,lattice)
  IMPLICIT NONE
  ! input and output arguments
  INTEGER, INTENT(IN) :: lupri,luerr,nbast
  TYPE(lvec_list_t),INTENT(INOUT) :: lattice
  TYPE(moleculeinfo),INTENT(INOUT),TARGET :: refcell
  TYPE(lssetting),INTENT(INOUT) :: setting 
  ! local variables
  INTEGER ::  idx

  call set_lstime_print(.false.)

  do idx=1,size(lattice%lvec) 
     call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(idx)%molecule,2,refcell, &
		  & 3,lattice%lvec(idx)%molecule,4)
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
!> \param refcell 		Molecule info reference cell.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix.
SUBROUTINE find_cutoff_twop(lupri,luerr,setting,nbast,lattice, &
		& refcell,numvecs,nfdensity)
	IMPLICIT NONE
	! input and output arguments
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs ! nlayer 
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
		& lattice%max_layer

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
				call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(idx)%molecule,3)
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
							& lattice,lattice%max_layer)
					
						setting%samefrag=.false.
						setting%samemol=.false.

						! ToDo Should be turned on at some point
						call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(idx)%molecule, &
							& 3,lattice%lvec(index2)%molecule,2,lattice%lvec(newcell)%molecule,4)
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
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param ovl 		 		Overlap (array of matrices for each lattice cell.
!todo is it nec to calc entire ovl or can we make use of symmetries?
SUBROUTINE pbc_overlap_k(lupri,luerr,setting,natoms,nbast,lattice, &
		& refcell,numvecs,ovl)
	IMPLICIT NONE
	! input and output arguments
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(lvec_list_t),INTENT(INOUT) :: lattice
	TYPE(lssetting),INTENT(INOUT) :: setting 
	TYPE(matrix),TARGET :: ovl(numvecs)
	! local variables
        REAL(realk) :: max_elm,max_elm_layer
	INTEGER :: i,j,il1,il2,il3
	INTEGER :: idx,refindex,layer,layer_old
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

        layer_old = 0
        max_elm_layer = 0.0_realk
	DO idx=1,numvecs
		lattice%lvec(idx)%ovl_computed=.false.

		!Doing translations
		call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(idx)%molecule,3)
		setting%samemol(1,3)=.false.
		setting%samemol(3,1)=.false.

		call find_latt_vectors(idx,il1,il2,il3,lattice)
                layer = max(abs(il1),abs(il2))
                layer = max(layer,abs(il3))

                if ( layer .ne. layer_old) then
                  !Checks if there is need 
                  !in calculating integrals for
                  !this layer
                  if(max_elm_layer .le. lattice%intthr) exit
                  layer_old = layer
                  max_elm_layer = 0.0_realk
                endif
		gab1=lattice%lvec(idx)%maxgab
		if(gab1 .ge. lattice%realthr) then
			maxl1=max(abs(il1),maxl1)
			maxl2=max(abs(il2),maxl2)
			maxl3=max(abs(il3),maxl3)
			lattice%lvec(idx)%ovl_computed=.true.
                        write(*,*) 'idx in overlap',idx,il1,il2,il3

			if(.not. lattice%store_mats) then
				call mat_init(ovl(idx),nbast,nbast)
			endif
			call mat_init(lattice%lvec(idx)%oper(1),nbast,nbast)
			call mat_zero(lattice%lvec(idx)%oper(1))
			! get the overlap matrix for cell between reference and cell l
			call II_get_overlap(lupri,luerr,setting,lattice%lvec(idx)%oper(1))

                        call mat_abs_max_elm(lattice%lvec(idx)%oper(1),max_elm)
                        max_elm_layer=max(max_elm_layer,max_elm)

                        !write(*,*) 'Overlap computed for',il1,il2,il3
			if(lattice%store_mats) then !todo necc? are the matrices ever stored??
				call pbc_get_file_and_write(lattice,nbast,nbast,idx,1,1, &
					& '            ')!1 refers to overlap
			else
				call mat_copy(1.0_realk,lattice%lvec(idx)%oper(1),ovl(idx))
			end if

			call mat_free(lattice%lvec(idx)%oper(1))
		endif
                write(*,*) 'overlap',idx,lattice%lvec(idx)%ovl_computed
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
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix
!> \param f1 		 		Fock matrix. Initialized with kinetic part here.
!> \param e_kin 			Kinetic energy (optional, calculated if present)
SUBROUTINE pbc_kinetic_k(lupri,luerr,setting,natoms,nbast,lattice, &
		& refcell,numvecs,nfdensity,f_1,e_kin)
	IMPLICIT NONE
	! input and output variables
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(LSSETTING),INTENT(INOUT) :: setting 
	TYPE(lvec_list_t),INTENT(INOUT) ::lattice
	TYPE(matrix),INTENT(IN) :: nfdensity(numvecs)
	TYPE(matrix),INTENT(INOUT),TARGET :: f_1(numvecs)
	REAL(realk),INTENT(OUT),OPTIONAL :: e_kin
	! local variables
        REAL(realk) :: max_elm,max_elm_layer
	LOGICAL :: calc_e_kin
	INTEGER :: il1,il2,il3,layer,layer_old
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

        layer_old = 0
        max_elm_layer = 0.0_realk
	do idx=1,numvecs

		lattice%lvec(idx)%f1_computed=.false.
		if(.not. lattice%lvec(idx)%is_redundant) then
			call find_latt_vectors(idx,il1,il2,il3,lattice)
			call find_latt_index(indred,-il1,-il2,-il3,lattice,lattice%max_layer)
			if(il1**2+il2**2+il3**2 .gt. 0) lattice%lvec(indred)%is_redundant =.true.
			call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(idx)%molecule,3)
			setting%samemol(1,3)=.false.
			setting%samemol(3,1)=.false.
			gab1=lattice%lvec(idx)%maxgab

                        layer = max(abs(il1),abs(il2))
                        layer = max(layer,abs(il3))

                        if ( layer .ne. layer_old) then
                          !Checks if there is need 
                          !in calculating integrals for
                          !this layer
                          if(max_elm_layer .le. lattice%intthr) exit
                          layer_old = layer
                          max_elm_layer = 0.0_realk
                        endif

			if(gab1 .ge. lattice%realthr) then
				lattice%lvec(idx)%f1_computed=.true.
				lattice%lvec(indred)%f1_computed=.true.
				call mat_init(lattice%lvec(idx)%oper(1),nbast,nbast)
				call mat_zero(lattice%lvec(idx)%oper(1))
				call mat_init(f_1(idx),nbast,nbast)
				call mat_zero(f_1(idx))
				call II_get_kinetic(lupri,luerr,setting,lattice%lvec(idx)%oper(1))

                                call mat_abs_max_elm(lattice%lvec(idx)%oper(1),max_elm)
                                max_elm_layer=max(max_elm_layer,max_elm)

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
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix
!> \param f1 		 		Fock matrix. Nuc part added here.
!> \param e_en 			Kinetic energy (optional, calculated if present)
SUBROUTINE pbc_nucattrc_k(lupri,luerr,setting,natoms,nbast,lattice, &
		& refcell,numvecs,nfdensity,f_1,e_en)
	IMPLICIT NONE
	! input and output variables
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(lssetting),INTENT(INOUT) :: setting 
	TYPE(lvec_list_t),intent(INOUT) :: lattice
	TYPE(matrix),INTENT(IN) :: Nfdensity(numvecs)
	TYPE(matrix),INTENT(INOUT),TARGET :: f_1(numvecs)
	TYPE(matrix) :: H,H1
	REAL(realk),INTENT(OUT),OPTIONAL :: E_en
	! local variables
        REAL(realk) :: max_elm,max_elm_layer
	INTEGER :: j,k,il11,il12,il13,il21,il22,il23
	INTEGER :: index1,index2,refindex
	INTEGER :: iunit,indred,layer,layer_old
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

        layer_old = 0
        max_elm_layer = 0.0_realk
        do index1=1,numvecs

		call find_latt_vectors(index1,il11,il12,il13,lattice)
		!if(abs(il11) .gt. lattice%nneighbour) CYCLE
		!if(abs(il12) .gt. lattice%nneighbour) CYCLE
		!if(abs(il13) .gt. lattice%nneighbour) CYCLE
                layer = max(abs(il11),abs(il12))
                layer = max(layer,abs(il13))

                if ( layer .ne. layer_old) then
                  !Checks if there is need 
                  !in calculating integrals for
                  !this layer
                  if(max_elm_layer .le. lattice%intthr) exit
                  layer_old = layer
                  max_elm_layer = 0.0_realk
                endif

                lattice%lvec(index1)%Vz_computed=.false.
                if(.not. lattice%lvec(index1)%is_redundant) then
			call find_latt_index(indred,-il11,-il12,-il13,lattice, &
				& lattice%max_layer)
			if(il11**2+il12**2+il13**2 .gt. 0) then
				lattice%lvec(indred)%is_redundant =.true.
			end if
                        gab1=lattice%lvec(index1)%maxgab

			if(gab1 .ge. lattice%realthr) then	
				lattice%lvec(index1)%Vz_computed=.true.
				lattice%lvec(indred)%Vz_computed=.true.
				call mat_init(lattice%lvec(index1)%oper(1),nbast,nbast)
				call mat_zero(lattice%lvec(index1)%oper(1))
				do index2=1,numvecs
					call find_latt_vectors(index2,il21,il22,il23,lattice)
                                        if(abs(il21) .gt. lattice%nf) CYCLE 
                                        if(abs(il22) .gt. lattice%nf) CYCLE
                                        if(abs(il23) .gt. lattice%nf) CYCLE

                                        call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(index1)%molecule, &
                                          & 2,lattice%lvec(index2)%molecule,3)

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

                                call mat_abs_max_elm(lattice%lvec(index1)%oper(1),max_elm)
                                max_elm_layer=max(max_elm_layer,max_elm)
                                                
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
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix
!> \param g2 		 		Fock matrix. e-e rep part added.
!> \param e_J 				El. repulsion (optional, calculated if present)
SUBROUTINE pbc_electron_rep_k(lupri,luerr,setting,natoms,nbast, &
		& lattice,refcell,numvecs,nfdensity,g_2,E_J)
	IMPLICIT NONE
	! input and output
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms 
	TYPE(matrix),intent(in),DIMENSION(numvecs) :: nfdensity
	TYPE(matrix),target,intent(inout) :: g_2(numvecs)
	TYPE(lvec_list_t),intent(inout) ::lattice
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(lssetting),INTENT(INOUT) :: setting 
	TYPE(MATRIX),pointer :: F_tmp(:)
	REAL(realk),INTENT(INOUT),OPTIONAL :: E_J
	! local variables
	REAL(realk) :: valmax
        REAL(realk) :: max_elm,max_elm_layer
	INTEGER :: i,j,il21,il22,il23,iunit
	INTEGER :: il31,il32,il33, newcell,il1,il2,il3
	INTEGER :: index1,index2,index3,maxl1,maxl2,maxl3
	INTEGER :: l1,l2,l3,indred,layer,layer_old
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

        layer_old = 0
        max_elm_layer = 0.0_realk
	do index1=1,numvecs

          !if(.not. lattice%lvec(index1)%is_redundant) then
          il1=int(lattice%lvec(index1)%lat_coord(1))
          il2=int(lattice%lvec(index1)%lat_coord(2))
          il3=int(lattice%lvec(index1)%lat_coord(3))

          lattice%lvec(index1)%g2_computed=.false.
          lattice%lvec(index1)%J_computed=.false.
          gab1=lattice%lvec(index1)%maxgab

          layer = max(abs(il1),abs(il2))
          layer = max(layer,abs(il3))

          if ( layer .ne. layer_old) then
            !Checks if there is need 
            !in calculating integrals for
            !this layer
            if(max_elm_layer .le. lattice%intthr) exit
            layer_old = layer
            max_elm_layer = 0.0_realk
          endif


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
            write(*,*) 'J computed for',il1,il2,il3

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
              !Not sure with this one I will not use densities for 
              !lattice cells where J was too small
              if(index2 .gt. 1 .and. .not. lattice%lvec(index2)%J_computed)&
                &CYCLE

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

                  gabmaxsum=gab1+gab2
                  call mat_abs_max_elm(nfdensity(index2),valmax)
                  if(valmax .gt. 0._realk) valm1=int(log10(valmax),kind=short)

                  if(valm1+gabmaxsum  .ge. lattice%realthr) then
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
                      & lattice%lvec(index1)%molecule, 2,&
                      &lattice%lvec(index3)%molecule,3,& 
                      & lattice%lvec(newcell)%molecule,4)
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

              call mat_abs_max_elm(lattice%lvec(index1)%oper(2),max_elm)
              max_elm_layer=max(max_elm_layer,max_elm)

              if(.not. lattice%store_mats) then
                call mat_init(g_2(index1),nbast,nbast)
                call mat_copy(1._realk,lattice%lvec(index1)%oper(2),g_2(index1))
              endif

              if (calc_e_j) then
                if((abs(il1).le.lattice%ndmat .and. abs(il2).le.lattice%ndmat) &
                  & .and. abs(il3).le.lattice%ndmat) then
                if(nfdensity(index1)%init_magic_tag&
                  &.EQ.mat_init_magic_value) THEN

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
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
!> \param nfdensity 		Density matrix
!> \param g2 		 		Fock matrix. e-e rep part added.
!> \param e_k 				Exchange (optional, calculated if present)
SUBROUTINE pbc_exact_xc_k(lupri,luerr,setting,natoms,nbast, &
		& lattice,refcell,numvecs,nfdensity,g_2,E_K)
	IMPLICIT NONE
	! input and output
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,natoms
	TYPE(moleculeinfo),INTENT(INOUT) :: refcell
	TYPE(matrix),INTENT(IN) :: nfdensity(numvecs)
	TYPE(matrix),INTENT(INOUT) :: g_2(numvecs)
	TYPE(LSSETTING),INTENT(INOUT) :: SETTING 
	TYPE(lvec_list_t),INTENT(INOUT) ::lattice
	REAL(realk),INTENT(INOUT),OPTIONAL :: E_K
	! local variables
	TYPE(MATRIX),pointer :: K_tmp(:)
	TYPE(MATRIX) :: Kx
	INTEGER :: il21,il22,il23,gabl1,gabl2,gabl3
	INTEGER :: il31,il32,il33,newcell,il1,il2,il3
	INTEGER :: index1,index2,index3,gabind,layer,layer_old
	INTEGER :: l1,l2,l3,maxl1,maxl2,maxl3,indred
	INTEGER(short) :: gab1,gab2,maxgabsum
	INTEGER, SAVE :: iter=0
	INTEGER(short) :: valm1
	REAL(realk) :: valmax
        REAL(realk) :: max_elm,max_elm_layer
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

        layer_old = 0
        max_elm_layer = 0.0_realk
	do index1=1,numvecs

          !call find_latt_vectors(index1,il1,il2,il3,lattice)
          il1=int(lattice%lvec(index1)%lat_coord(1))
          il2=int(lattice%lvec(index1)%lat_coord(2))
          il3=int(lattice%lvec(index1)%lat_coord(3))
          !if(abs(il1) .gt. lattice%kx1) CYCLE
          !if(abs(il2) .gt. lattice%kx2) CYCLE
          !if(abs(il3) .gt. lattice%kx3) CYCLE
          lattice%lvec(index1)%Kx_computed=.false.

          layer = max(abs(il1),abs(il2))
          layer = max(layer,abs(il3))

          if ( layer .ne. layer_old) then
            !Checks if there is need 
            !in calculating integrals for
            !this layer
            if(max_elm_layer .le. lattice%intthr) exit
            layer_old = layer
            max_elm_layer = 0.0_realk
          endif

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
                !Not sure with this one I will not use densities for 
                !lattice cells where Kx was too small
                if(index3 .gt. 1 .and. .not. lattice%lvec(index3)%Kx_computed)&
                  &CYCLE

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

                if(maxgabsum+valm1 .ge. lattice%realthr) then
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
                    & lattice%lvec(index1)%molecule,3,&
                    &lattice%lvec(index2)%molecule,2, &
                    & lattice%lvec(newcell)%molecule,4)
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

              !finds max element layer for layer
              call mat_abs_max_elm(lattice%lvec(index1)%oper(1),max_elm)
              max_elm_layer=max(max_elm_layer,max_elm)
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
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
SUBROUTINE pbc_overlap_int(lupri,luerr,setting,nbast,lattice,refcell,numvecs)
  IMPLICIT NONE
  ! input and output variables
  INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs
  TYPE(lssetting),INTENT(INOUT) :: SETTING 
  TYPE(moleculeinfo),INTENT(INOUT) :: refcell
  TYPE(lvec_list_t),INTENT(INOUT) :: lattice
  ! local variables
  TYPE(matrix)  :: S
  INTEGER :: il1,il2,il3
  INTEGER ::  indx

  call mat_init(S,nbast,nbast)
  call mat_zero(S)
  write(lupri,*) 'Number of lattice vectors ', numvecs

  do indx=1,numvecs
	  call typedef_setmolecules(setting,refcell,1,lattice%lvec(indx)%molecule,3)
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

!> \author Johannes Rekkedal
!> \date 2013
!> \brief For computing Kin. E. integrals (i|T|j) 
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param natoms 			Number of atoms.
!> \param nbast 			Number of basis func.
!> \param fock_mtx 		Fock matrix.
!> \param sizef 			Num col/row in fock matrix.
!> \param lattice 		Information about the lattice.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
SUBROUTINE pbc_kinetic_int(lupri,luerr,setting,natoms,nbast,fock_mtx,sizef, &
		& lattice,refcell,numvecs)
	IMPLICIT NONE
	! input and output variables
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,sizef,numvecs,natoms
	COMPLEX(complexk), INTENT(INOUT) :: fock_mtx(sizef*sizef)
	TYPE(lssetting),INTENT(INOUT) :: setting 
	TYPE(moleculeinfo),INTENT(IN) :: refcell
	TYPE(lvec_list_t),intent(INOUT) ::lattice
	! local variables
	TYPE(MATRIX) :: kin
	INTEGER :: il1,il2,il3
	INTEGER :: indx,refindex
!	REAL(realk) :: fock_tmp(nbast*nbast)

	write(lupri,*) 'Starting pbc_kinetic_int'

	call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)
	write(lupri,*) 'reference cell index=', refindex,lattice%max_layer

	call mat_init(kin,nbast,nbast)
	call mat_zero(kin)
!	write(lupri,*) 'nbast',nbast,kin%ncol,kin%nrow

	do indx=1,numvecs
		call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(indx)%molecule,3)
		call find_latt_vectors(indx,il1,il2,il3,lattice)
		if(abs(il1) .gt. lattice%nneighbour) CYCLE
		if(abs(il2) .gt. lattice%nneighbour) CYCLE
		if(abs(il3) .gt. lattice%nneighbour) CYCLE

		setting%samemol(1,3)=.false.
		setting%samemol(3,1)=.false.
		call II_get_kinetic(lupri,luerr,setting,kin)

!		fock_tmp=kin%elms
!		write(lupri,*) fock_tmp
!		fock_mtx=fock_tmp

	end do

	call mat_free(kin)
	write(lupri,*) 'Finished pbc_kinetic_int'

END SUBROUTINE pbc_kinetic_int

!> \author Johannes Rekkedal
!> \date 2013
!> \brief For computing Nuc. E. integrals 
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param natoms 			Number of atoms.
!> \param nbast 			Number of basis func.
!> \param fock_mtx 		Fock matrix.
!> \param sizef 			Num col/row in fock matrix.
!> \param lattice 		Information about the lattice.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice
SUBROUTINE pbc_nucattrc_int(lupri,luerr,setting,natoms,nbast,fock_mtx, &
		& sizef,lattice,refcell,numvecs)
	IMPLICIT NONE
	! input and output
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,sizef,numvecs,natoms
	COMPLEX(complexk), INTENT(INOUT) :: fock_mtx(sizef*sizef)
	TYPE(lssetting),INTENT(INOUT)   :: setting 
	TYPE(moleculeinfo),INTENT(IN) :: refcell
	TYPE(lvec_list_t),intent(INOUT) ::lattice
	! local variables
	TYPE(matrix)  :: H,H1
	INTEGER :: il11,il12,il13,il21,il22,il23
	INTEGER ::  index1,index2,refindex
!	REAL(realk) :: fock_tmp(nbast*nbast)

	call mat_init(H,nbast,nbast)
	call mat_zero(H)
	call mat_init(H1,nbast,nbast)
	call mat_zero(H1)

	call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)
	write(lupri,*) 'reference cell index=', refindex

	do index1=1,numvecs

		call find_latt_vectors(index1,il11,il12,il13,lattice)
		if(abs(il11) .gt. lattice%nneighbour) CYCLE
		if(abs(il12) .gt. lattice%nneighbour) CYCLE
		if(abs(il13) .gt. lattice%nneighbour) CYCLE

		do index2=1,numvecs
			call find_latt_vectors(index2,il21,il22,il23,lattice)
			if(abs(il21) .gt. lattice%nneighbour) CYCLE !uncomment this
			if(abs(il22) .gt. lattice%nneighbour) CYCLE
			if(abs(il23) .gt. lattice%nneighbour) CYCLE

			call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(index1)%molecule&
                          & ,2,lattice%lvec(index2)%molecule,3)

			setting%samemol(1,3)=.false.
			setting%samemol(3,1)=.false.
			setting%samemol(1,2)=.false.
			setting%samemol(2,1)=.false.
			setting%samemol(3,2)=.false.
			setting%samemol(2,3)=.false.

			call II_get_nucel_mat(lupri,luerr,setting,H1)
			call mat_add(1E0_realk,H,1E0_realk,H1,H)
		enddo

!    	fock_tmp=H1%elms
		fock_mtx=fock_mtx+H1%elms
		call mat_print(H1,1,H1%nrow,1,H1%ncol,6)
		call mat_zero(H1)

	enddo

	call mat_free(H)
	call mat_free(H1)
	write(lupri,*) 'finished with pbc_nucattrc_int'

END SUBROUTINE pbc_nucattrc_int

!> \author Johannes Rekkedal
!> \date 2013
!> \brief For computing Nuc. repulsions.
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param natoms 			Number of atoms.
!> \param lattice 		Information about the lattice.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice.
!> \param E_nn 			Nuclear repulsion potential energy.
SUBROUTINE pbc_nucpot(lupri,luerr,setting,natoms,lattice, &
		& refcell,numvecs,e_nn)
  IMPLICIT NONE
  ! input and output arguments
  INTEGER, INTENT(IN) :: lupri,luerr,numvecs,natoms
  TYPE(moleculeinfo),INTENT(INOUT) :: refcell
  TYPE(lvec_list_t),INTENT(INOUT) ::lattice
  TYPE(lssetting),INTENT(INOUT) :: setting 
  REAL(realk),INTENT(INOUT) :: E_nn
  ! local variables
  REAL(realk) :: kvec(3)
  INTEGER :: il1,il2,il3
  INTEGER :: indx,refindex
  REAL(realk) :: nucpot
	
  E_nn=0.0_realk

  call set_lstime_print(.false.)
  write(lupri,*) 'Number of lattice vectors ', numvecs
  !finds the lattice index for the reference cell.
  call find_latt_index(refindex,0,0,0,lattice,lattice%max_layer)

  !loop over lattice cells
  do indx=1,numvecs

	  !Doing translations
     call typedef_setmolecules(setting,refcell,1,lattice%lvec(indx)%molecule,3)
     setting%samemol(1,3)=.false.
     setting%samemol(3,1)=.false.

     call find_latt_vectors(indx,il1,il2,il3,lattice)
     !So that we do not consider negligible integrals
     if(abs(il1) .gt. lattice%nf) CYCLE
     if(abs(il2) .gt. lattice%nf) CYCLE
     if(abs(il3) .gt. lattice%nf) CYCLE

     if(indx .eq. refindex) then
       call II_get_nucpot(lupri,luerr,setting,nucpot)
     else
       call pbc_get_nucpot(lupri,luerr,setting,nucpot)
     endif
     E_nn=E_nn+nucpot

  end do

END SUBROUTINE pbc_nucpot

!> \brief Calculates the nuclear repulsion energy contribution
!> \author S. Reine and T. Kjaergaard
!> \modified by J. Rekkedal
!> \date 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nucpot the nuclear repulsion energy contribution
SUBROUTINE pbc_get_nucpot(lupri,luerr,setting,nucpot)
	IMPLICIT NONE
	TYPE(lssetting) :: setting
	INTEGER :: usemat
	INTEGER :: lupri,luerr
	REAL(realk) :: nucpot
	INTEGER :: i,j
	REAL(realk) :: pq(3),distance
	LOGICAL :: nobqbq

	call time_II_operations1()
	nobqbq=setting%scheme%nobqbq
	nucpot=0.0e0_realk

	do i=1,setting%molecule(1)%p%natoms
		if(setting%molecule(1)%p%atom(i)%phantom) CYCLE

		do j=1,setting%molecule(1)%p%natoms

			if(setting%molecule(1)%p%atom(j)%phantom) CYCLE
			if(setting%molecule(1)%p%atom(i)%pointcharge &
			  	& .and. setting%molecule(1)%p%atom(j)%pointcharge &
				& .and. nobqbq) CYCLE

			pq(1)=setting%molecule(1)%p%atom(i)%center(1) &
				& - setting%molecule(3)%p%atom(j)%center(1)
			pq(2)=setting%molecule(1)%p%atom(i)%center(2) & 
				& - setting%molecule(3)%p%atom(j)%center(2)
			pq(3)=setting%molecule(1)%p%atom(i)%center(3) - &
				& setting%molecule(3)%p%atom(j)%center(3)

			distance=sqrt(pq(1)*pq(1)+pq(2)*pq(2)+pq(3)*pq(3))
			nucpot=nucpot+setting%molecule(1)%p%atom(i)%charge &
				& * setting%molecule(3)%p%atom(j)%charge/distance
		enddo
	enddo
	nucpot=nucpot/2.
	call time_ii_operations2(job_ii_get_nucpot)

END SUBROUTINE pbc_get_nucpot

!> \brief Calculates the nuclear repulsion energy contribution
!> \author S. Reine and T. Kjaergaard
!> \modified by J. Rekkedal
!> \date 2012
!> \param lupri 			Logical print unit output file.
!> \param luerr 			Logical print unit error file.
!> \param setting 		Integral settings.
!> \param natoms 			Number of atoms.
!> \param nbast 			Number of basis functions.
!> \param lattice 		Information about the lattice.
!> \param refcell 		Molecule info. Threat reference cell as molecule.
!> \param numvecs 		Number of unitcells in the BvK lattice.
!> \param nfdensity 		The density matrix.
!> \param nfsze 			Num densmat nearfield.
SUBROUTINE pbc_electron_rep(lupri,luerr,setting,natoms, &
		& nbast,lattice,refcell,numvecs,nfdensity,nfsze)
	IMPLICIT NONE
	! inout and output
	INTEGER, INTENT(IN) :: lupri,luerr,nbast,numvecs,nfsze,natoms
	TYPE(moleculeinfo),INTENT(IN) :: refcell
	TYPE(lssetting),INTENT(INOUT) :: setting 
	TYPE(matrix),INTENT(INOUT),DIMENSION(nfsze) :: nfdensity
	TYPE(lvec_list_t),INTENT(INOUT) ::lattice
	! local variables
	TYPE(MATRIX)  :: F(1),F_tmp(1)
	INTEGER :: il21,il22,il23,l1,l2,l3
	INTEGER :: il31,il32,il33,newcell
	INTEGER ::  index1,index2,index3
	INTEGER :: checknf,checknf1,checknf2,checknf3

	call mat_init(F(1),nbast,nbast)
	call mat_zero(F(1))
	call mat_init(F_tmp(1),nbast,nbast)
	call mat_zero(F_tmp(1))

	do index1=1,numvecs
!		call find_latt_vectors(index1,il1,il2,il3,lattice)
		do index2=1,numvecs

			call find_latt_vectors(index2,il21,il22,il23,lattice)
			if(abs(il21) .gt. lattice%nneighbour) CYCLE
			if(abs(il22) .gt. lattice%nneighbour) CYCLE
			if(abs(il23) .gt. lattice%nneighbour) CYCLE

			do index3=1,numvecs
				call find_latt_vectors(index3,il31,il32,il33,lattice)

				l1=il21+il31
				if(abs(l1) .gt. lattice%max_layer) CYCLE
				l2=il22+il32
				if(abs(l2) .gt. lattice%max_layer) CYCLE 
				l3=il23+il33
				if(abs(l3) .gt. lattice%max_layer) CYCLE 
				checknf1=l1-il21
				if(abs(checknf1) .gt. n_neighbour) CYCLE
				checknf2=l2-il22
				if(abs(checknf2) .gt. n_neighbour) CYCLE 
				checknf3=l3-il23
				if(abs(checknf3) .gt. n_neighbour) CYCLE 

				call find_latt_index(newcell,l1,l2,l3,lattice,lattice%max_layer)
				call find_latt_index(checknf,checknf1,checknf2,checknf3,lattice,&
					& lattice%nneighbour)
				!Changes setting to point at different lattice cells
				call TYPEDEF_setmolecules(setting,refcell,1,lattice%lvec(index1)%molecule,2,&
					& lattice%lvec(index2)%molecule,3,lattice%lvec(newcell)%molecule,4)

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
				call mat_add(1E0_realk,F(1),1E0_realk,F_tmp(1),F(1))

			enddo
		enddo

		CALL mat_print(F(1),1,F_tmp(1)%nrow,1,F_tmp(1)%ncol,6)
		call mat_zero(F(1))

		write(lupri,*) 'index1',index1, numvecs
	ENDDO
	call mat_free(F(1))
	call mat_free(F_tmp(1))

	write(lupri,*) 'finished with pbc_elrep_int'
	call LSQuit('finished with pbc_elrep_int',lupri)

END SUBROUTINE pbc_electron_rep

!> \brief ???
!> \author J. Rekkedal
!> \date 2012
!> \param lupri 			Logical print unit output file.
!> \param nbast 			Number of basis functions.
!> \param fock_mtx 		The fock matrix.
!> \param sizef 			Fock matrix num col/row.
!> \param cut 				
!> \param lattice 		Information about the lattice.
SUBROUTINE pbc_complete_Fock_mtx(lupri,nbast,fock_mtx,sizef,cut,lattice)
	IMPLICIT NONE
	! input and output variables
	INTEGER, INTENT(IN) :: lupri,nbast,sizef
	REAL(realk),INTENT(INOUT) :: fock_mtx(sizef,sizef)
	REAL(realk),INTENT(IN) :: cut
	TYPE(lvec_list_t),INTENT(IN) ::lattice
	! Local variables
	REAL(realk) :: distance
	INTEGER :: num_latvectors
	INTEGER :: cellij,i,j
	INTEGER :: diff1,diff2,diff3
	INTEGER :: il1,il2,il3,jl1,jl2,jl3
	REAL(realk) ::  latt_vec_std(3),tmp_vec(3)
	INTEGER :: refcell

	tmp_vec(:)=1.0_realk 
	latt_vec_std(:)=0.0_realk
	call latt_2_std_coord(tmp_vec,latt_vec_std,lattice%ldef%avec)
	
	call find_latt_index(refcell,0,0,0,lattice,lattice%max_layer)
	write(lupri,*) 'refcelle', refcell

        num_latvectors=size(lattice%lvec)

	do i=1,num_latvectors
		if(i .eq. refcell) CYCLE
		fock_mtx(nbast*i+1-nbast:i*nbast,nbast*i+1-nbast:i*nbast)=&
			& fock_mtx(refcell*nbast-nbast+1:refcell*nbast, refcell*nbast-nbast+1:refcell*nbast)
		call find_latt_vectors(i,il1,il2,il3,lattice)
		do j=1,i
			if(j .eq. refcell) CYCLE
			if(j .eq. i) CYCLE
			call find_latt_vectors(j,jl1,jl2,jl3,lattice)
			diff1=il1-jl1
			if(abs(diff1) .gt. lattice%max_layer) CYCLE !fixme should be max_layer
			if(abs(diff1) .gt. lattice%nneighbour) CYCLE !fixme should be threshold
			diff2=il2-jl2
			if(abs(diff2) .gt.lattice%max_layer) CYCLE !fixme should be max_layer
			if(abs(diff2) .gt.lattice%nneighbour) CYCLE !fixme should be threshold
			diff3=il3-jl3
			if(abs(diff3) .gt. lattice%max_layer) CYCLE !fixme should be max_layer
			if(abs(diff3) .gt. lattice%nneighbour) CYCLE !fixme should be threshold

			call find_latt_index(cellij,diff1,diff2,diff3,lattice,lattice%max_layer)
			call calc_distance(distance,lattice%lvec(cellij)%std_coord,latt_vec_std)

			if(distance .ge. cut) CYCLE
			fock_mtx(nbast*i+1-nbast:i*nbast,nbast*j+1-nbast:j*nbast)=&
				& fock_mtx(refcell*nbast-nbast+1:refcell*nbast, cellij*nbast-nbast+1:cellij*nbast)
			fock_mtx(nbast*j+1-nbast:j*nbast,nbast*i+1-nbast:i*nbast)=&
				& fock_mtx(nbast*i+1-nbast:i*nbast, nbast*j+1-nbast:j*nbast)
		enddo
	enddo
	write(lupri,*) 'FOCK MATRIX'
	call write_matrix(fock_mtx,sizef,sizef)

END SUBROUTINE pbc_complete_Fock_mtx

END MODULE pbc_interactions

